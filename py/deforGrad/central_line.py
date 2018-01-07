#! /usr/bin/env python

#
# Input: obj file from the simulation
# Output: obj file of the central line
#

import argparse
import numpy as np

class FiberBundle:

    def __init__(self, filename):
        self.vertices = []
        self.fibers = []
        edges = []
        with open(filename, 'r') as fin:
            for line in fin:
                tks = line.strip().split()
                if tks[0] == 'v':
                    self.vertices.append( np.array([ float(tks[i]) for i in xrange(1,4)]) )
                elif tks[0] == 'l':
                    edges.append( (int(tks[1])-1, int(tks[2])-1) )

        # now create individual fiber
        for e in edges:
            ll = len(self.fibers)
            found = False
            for ii in xrange(ll-1,-1,-1):
                if self.fibers[ii][-1] == e[0]:
                    self.fibers[ii].append( e[1] )
                    found = True
                    break
            if not found:
                self.fibers.append( [e[0], e[1]] )

        # ----------------------------
        self.nFibers = len( self.fibers ) # number of fibers
        print '%d fibers loaded' % self.nFibers
        self.fiberLen = len( self.fibers[0] )
        for ii in xrange(1,self.nFibers):
            if len( self.fibers[ii] ) != self.fiberLen:
                raise Exception('Inconsistent data.')

    def gen_central_line(self):
        self.line = np.zeros( (3, self.fiberLen) )
        for vi in xrange(self.fiberLen):
            v = np.zeros((3,))
            for jj in xrange(self.nFibers):
                v += self.vertices[ self.fibers[jj][vi] ]
            v /= self.nFibers
            self.line[:,vi] = v

    def output_obj(self, filename):
        (t,L) = self.line.shape
        with open(filename, 'w') as fout:
            for vi in xrange(L):
                print >> fout, 'v %f %f %f' % ( self.line[0,vi], self.line[1,vi], self.line[2,vi] )
            for vi in xrange(1,L):
                print >> fout, 'l %d %d' % (vi, vi+1)
    
    def smooth_curve(self, W, S):
        tt = W/2 + 1
        dw = 1. / float(tt)
        ws = [ (i+1)*dw for i in xrange(W) ]
        for ii in xrange(W):
            if ws[ii] > 1+1E-8:
                ws[ii] = 2 - ws[ii]
        
        sline = np.zeros(self.line.shape)
        for ii in xrange(self.line.shape[1]):
            # print '------'
            ww = 0.
            for jj in xrange(W):
                idx = jj - W/2 + ii
                if idx >= 0 and idx < self.line.shape[1]:
                    # print '(%d, %f)' % (idx, ws[jj])
                    sline[:,ii] += self.line[:,idx]*ws[jj]
                    ww += ws[jj]
            sline[:,ii] /= ww

        if S > 1:
            # ensure the cut from both sides is even
            xx = range(0, sline.shape[1], S)[-1]
            off = (sline.shape[1] - 1 - xx)/2
            #print range(off, sline.shape[1], S)
            self.line = sline[:, off:sline.shape[1]:S]
        else:
            self.line = sline

# ---------------------------------------------------------------------------------

def main(args):
    if args.range is None:
        bundle = FiberBundle(args.infile)
        bundle.gen_central_line()
        if args.smooth is not None:
            bundle.smooth_curve( args.smooth[0], args.smooth[1] )
        bundle.output_obj('curve.obj')
    else:
        for i in xrange(args.range[0], args.range[1]+1):
            fname = args.infile % i
            bundle = FiberBundle(fname)
            bundle.gen_central_line()
            if args.smooth is not None:
                bundle.smooth_curve( args.smooth[0], args.smooth[1] )
            bundle.output_obj('curve-%03d.obj' % i)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Extract the central line of fiber bundle',
            epilog='Changxi Zheng (cxz@cs.columbia.edu)')
    parser.add_argument('infile', metavar='filename', type=str, help='input filename (or filename pattern)')
    parser.add_argument('--range', metavar='N', type=int, nargs=2, help='range of input files')
    parser.add_argument('--smooth', metavar=('W', 'S'), type=int, nargs=2, help='Smooth and downsample the curve. W: smooth window size, S: downsample size. e.g., 13 3')

    args = parser.parse_args()
    main(args)
