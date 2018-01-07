#! /usr/bin/env python

import argparse
import numpy as np
from fiber_bundle import FiberBundle
from def_grad_est import DefGradEstimator

def downsample_def_grad(est, ids, fname):
    K = 2 # K must be even number
    assert K % 2 == 0
    W = ids[1] - ids[0] + K
    tt = W/2 + 1
    dw = 1. / float(tt)
    ws = [ (i+1)*dw for i in xrange(W) ]
    for ii in xrange(W):
        if ws[ii] > 1+1E-8:
            ws[ii] = 2 - ws[ii]

    dg = []
    for ii in xrange( len(ids)-1 ):
        totw = 0
        cdg = np.zeros( (3,3) )
        #print '---------------- ', ids[ii], ids[ii+1]
        for jj in xrange(W):
            ci = ids[ii]+jj-(K/2)
            if ci >= 0 and ci < len(est.central_dg):
                cdg += est.central_dg[ci] * ws[jj]
                totw += ws[jj]
                #print ci, ws[jj]
        dg.append( cdg/totw )

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot( range(len(dg)), [d[0,0] for d in dg], '.-' )
    ax.plot( range(len(dg)), [d[1,1] for d in dg], '.-' )
    ax.plot( range(len(dg)), [d[2,2] for d in dg], '.-' )
    #ax.plot( range(len(dg)), [d[0,1] for d in dg], '.-' )
    #ax.plot( range(len(dg)), [d[0,2] for d in dg], '.-' )
    #ax.plot( range(len(dg)), [d[1,2] for d in dg], '.-' )
    plt.show()

    with open(fname, 'w') as fout:
        for d in dg:
            print >> fout, ' '.join( [ '%.10f'%x for x in d.reshape(9) ] )


def main(args):
    # initial fiber bundle
    initFiber = FiberBundle( '%s.obj' % args.infile[0] )
    initFiber.gen_central_line()
    #initFiber.downsample(args.smooth[1])

    # current fiber bundle
    curFiber = FiberBundle( '%s.obj' % args.infile[1] )
    curFiber.gen_central_line()

    # smooth the curves
    initFiber.smooth_curve( args.smooth[0],1 )
    curFiber.smooth_curve( args.smooth[0],1 )

    fname = '%s_CL.obj' % args.infile[1]
    print 'Save original central line to %s ...' % fname
    curFiber.output_obj(fname)

    #est = DefGradEstimator( initFiber, curFiber, '%s.fe' % args.infile[1] )
    #est.average_dg()
    #est.ls_estimate_dg()

    est = DefGradEstimator( initFiber, curFiber, None )
    est.est_affine_dg()

    fname = '%s_CL.fe' % args.infile[1]
    print 'Save original def. grad. to %s ...' % fname
    est.save_dg(fname)

    ids = curFiber.downsample(args.smooth[1])
    fname = '%s_RED.obj' % args.infile[1]
    print 'Save reduced central line to %s (%d particles)...' % (fname, curFiber.numCentralParticles)
    curFiber.output_obj(fname)

    fname = '%s_RED.fe' % args.infile[1]
    downsample_def_grad(est, ids, fname)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Downsample the central line of a yarn and the deformation gradient. The results will stored in two files, one for the downsampled central line and another for the downsampled deformation gradient',
            epilog='Changxi Zheng (cxz@cs.columbia.edu)')
    parser.add_argument('infile', metavar='filename', type=str, nargs=2,
            help='<rest frame file>, <current frame file>. The rest and current frame files are specified without the filename extensions.')
    parser.add_argument('--smooth', metavar=('W', 'S'), type=int, nargs=2, default=(7,7),
            help='Smooth and downsample the curve. W: smooth window size, S: downsample size (e.g., 13 3). If not given, the default values are 7 7.')

    args = parser.parse_args()
    main(args)
