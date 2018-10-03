#! /usr/bin/env python

import numpy as np

class FiberBundle:

    def __init__(self, filename):
        self._vertices = []
        self._fibers = []
        edges = []
        with open(filename, 'r') as fin:
            for line in fin:
                tks = line.strip().split()
                if tks[0] == 'v':
                    self._vertices.append( np.array([ float(tks[i]) for i in range(1,4)]) )
                elif tks[0] == 'l':
                    edges.append( (int(tks[1])-1, int(tks[2])-1) )

        # now create individual fiber
        for e in edges:
            ll = len(self._fibers)
            for ii in range(ll-1,-1,-1):
                if self._fibers[ii][-1] == e[0]:
                    self._fibers[ii].append( e[1] )
                    break
            else:
                self._fibers.append( [e[0], e[1]] )

        # ----------------------------
        self._nFibers = len( self._fibers ) # number of fibers
#        print '%d fibers loaded' % self._nFibers
        self._fiberLen = len( self._fibers[0] )
        for ii in range(1,self._nFibers):
            assert len( self._fibers[ii] ) == self._fiberLen

    @property
    def fibers(self):
        return self._fibers

    @property
    def fiberLen(self):
        return self._fiberLen

    @property
    def numFibers(self):
        return self._nFibers

    @property
    def centralLine(self):
        return self._cLine

    @property
    def numCentralParticles(self):
        return self._cLine.shape[1]

    def fiber_vertex(self, fid, vid):
        return  self._vertices[ self._fibers[fid][vid] ]

    def gen_central_line(self):
        self._cLine = np.zeros( (3, self._fiberLen) )
        for vi in range(self._fiberLen):
            v = np.zeros((3,))
            for jj in range(self._nFibers):
                v += self._vertices[ self._fibers[jj][vi] ]
            v /= self._nFibers
            self._cLine[:,vi] = v

    def output_obj(self, filename):
        (t,L) = self._cLine.shape
        with open(filename, 'w') as fout:
            for vi in range(L):
                print >> fout, 'v %f %f %f' % ( self._cLine[0,vi], self._cLine[1,vi], self._cLine[2,vi] )
            for vi in range(1,L):
                print >> fout, 'l %d %d' % (vi, vi+1)
    
    def smooth_curve(self, W, S):
        tt = W/2 + 1
        dw = 1. / float(tt)
        ws = [ (i+1)*dw for i in range(W) ]
        for ii in range(W):
            if ws[ii] > 1+1E-8:
                ws[ii] = 2 - ws[ii]
        
        sline = np.zeros(self._cLine.shape)
        for ii in range(self._cLine.shape[1]):
            # print '------'
            ww = 0.
            for jj in range(W):
                idx = jj - W/2 + ii
                if idx >= 0 and idx < self._cLine.shape[1]:
                    # print '(%d, %f)' % (idx, ws[jj])
                    sline[:,ii] += self._cLine[:,idx]*ws[jj]
                    ww += ws[jj]
            sline[:,ii] /= ww

        if S > 1:
            # ensure the cut from both sides is even
            xx = range(0, sline.shape[1], S)[-1]
            off = (sline.shape[1] - 1 - xx)/2
            #print range(off, sline.shape[1], S)
            self._cLine = sline[:, off:sline.shape[1]:S]
        else:
            self._cLine = sline

    def downsample(self, S):
        # make sure it hasn't been downsampled before
        assert self._fiberLen == self._cLine.shape[1]

        N = self._cLine.shape[1]
        if S > 1:
            # ensure the cut from both sides is even
            xx = range(0, N, S)[-1]
            off = (N - 1 - xx)/2
            #print range(off, sline.shape[1], S)
            #print '%d-%d:   ' % (0, N-1), range(off, N, S)
            self._cLine = self._cLine[:, off:N:S]
            return range(off, N, S)
        else:
            return range(N)

