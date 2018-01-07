#! /usr/bin/env python

import argparse
import numpy as np
from fiber_bundle import FiberBundle

def rotationFromVectors(a, b):
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    assert s > 0
    c = np.dot(a, b)
    V = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
    return np.identity(3) + V + np.dot(V, V)/(1.0 + c)


def loadYarnCenter(fname, n):
    vtx = []
    with open(fname) as fin:
        while True:
            line = fin.readline()
            if line[0] == 'v':
                val = np.array([float(x) for x in line[1:].strip().split()])
                assert len(val) == 3 or len(val) == 4
                vtx.append(val[0 : 3])
                if len(vtx) == n + 1:
                    break
    return vtx


def main(args):
    n = args.n

    pts0 = loadYarnCenter("%s.obj" % args.infile[0], n)

    bundle = FiberBundle("%s.obj" % args.infile[1])
    bundle.gen_central_line()
    idx = np.array(bundle.downsample(args.s))
    assert len(idx) == n + 1

    pts1 = loadYarnCenter("%s.obj" % args.infile[2], n)

    dg = []
    with open("%s.fe" % args.infile[3], "r") as fin:
        for line in fin.readlines():
            val = [float(x) for x in line.strip().split()]
            val = np.reshape(val, (3, 3))
            dg.append(val)
    assert len(dg) == n

    fidx = 0.5*(idx[0 : -1] + idx[1 :])

    centers0 = [None]*n
    centers1 = [None]*n
    for i in range(0, n):
        centers0[i] = 0.5*(pts0[i] + pts0[i + 1])
        centers1[i] = 0.5*(pts1[i] + pts1[i + 1])

    R = [None]*n
    t = [None]*n

    for i in range(0, n):
        e = pts0[i + 1] - pts0[i]
        e /= np.linalg.norm(e)
        R[i] = np.dot(dg[i], rotationFromVectors(e, np.array([1.0, 0.0, 0.0])))
        t[i] = centers1[i] - np.dot(R[i], centers0[i])



        
    fname = args.o + (".txt" if args.mitsuba else ".obj")
    with open("physical_150.txt", "w") as fout1:
        with open(fname, "w") as fout:
            if args.mitsuba:
                print >> fout, bundle.numFibers
            for i in range(0, bundle.numFibers):
                if args.mitsuba:
                    print >> fout, bundle.fiberLen
                for j in range(0, bundle.fiberLen):
                    k = np.searchsorted(fidx, j)
                   
                    if k == 0:
                        R0 = R[0]
                        t0 = t[0]
                    elif k >= n:
                        R0 = R[n - 1]
                        t0 = t[n - 1]
                    else:
                        w = (fidx[k] - j)/(fidx[k] - fidx[k - 1])
                        R0 = w*R[k - 1] + (1.0 - w)*R[k]
                        t0 = w*t[k - 1] + (1.0 - w)*t[k]
                    #####
                    if i==0:
                        fout1.writelines('%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f ' % (R0[0,0], R0[0,1], R0[0,2], R0[1,0], R0[1,1], R0[1,2], R0[2,0], R0[2,1], R0[2,2]) )
                        fout1.writelines('%.8f %.8f %.8f \n' % (t0[0]*0.25, t0[1]*0.25, t0[2]*0.25 ) )
                    ######
                    v = bundle.fiber_vertex(i, j)
                    v = np.dot(R0, v) + t0
                    if args.mitsuba:
                        print >> fout, "%.6f %.6f %.6f" % (v[0], v[1], v[2])
                    else:
                        print >> fout, "v %.6f %.6f %.6f" % (v[0], v[1], v[2])
            if not args.mitsuba:
                for i in range(0, bundle.numFibers):
                    for j in range(0, bundle.fiberLen - 1):
                        k = i*bundle.fiberLen + j + 1
                        print >> fout, "l %d %d" % (k, k + 1)
    
        
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Deform a bundle of fibers using simulated yarn-level deformation gradients',
            epilog='Shuang Zhao (shzz@ics.uci.edu)')

    parser.add_argument('infile', metavar='filename', type=str, nargs=4,
            help='<rest yarn shape (obj)>, <rest fiber shape (obj)>, <current yarn shape (obj)>, and <current yarn d.g. (fe)>. All filenames should be specified without extensions.')
    parser.add_argument('n', metavar='N', type=int,
            help='Number of segments (NOT vertices) of the downsampled yarn centerline.')

    parser.add_argument('-o', metavar='filename', type=str, default='output',
            help="Output filename without extension. If not given, the default value is 'output'.")
    parser.add_argument('-s', metavar='S', type=int, default=3,
            help='Downsample size. If not given, the default value is 3.')
    parser.add_argument('--mitsuba', action='store_false',
            help='Output deformed fiber bundles in Mitsuba format (instead of obj).')

    args = parser.parse_args()
    main(args)
