#! /usr/bin/env python

import argparse
import numpy as np
import numpy.linalg as linalg
from fiber_bundle import FiberBundle

from quat_utils import quaternion_from_two_vecs

class DefGradEstimator:

    def __init__(self, initFiber, curFiber, dgFile):
        self._est_type = None
        self._dg = None
        n = 1
        if dgFile is not None:
            self._dg = []
            with open(dgFile, 'r') as fin:
                for line in fin:
                    tks = line.strip().split()
                    assert len(tks) == 9
                    self._dg.append( np.array([ float(x) for x in tks ]).reshape( (3,3) ) )
                    #print linalg.det(self._dg[-1])
                    if len( self._dg ) == n*curFiber.fiberLen - 1:
                        self._dg.append( None )
                        n += 1
            # check the number of DGs
            #tot = reduce( (lambda x, y: x+y), [len(x)-1 for x in curFiber.fibers] )
            assert len(self._dg) == curFiber.numFibers * curFiber.fiberLen

        self._initFiber = initFiber
        self._curFiber  = curFiber
        assert self._initFiber.numFibers == self._curFiber.numFibers
        assert self._initFiber.fiberLen == self._curFiber.fiberLen

    @property
    def central_dg(self):
        return self._cDG

    def average_dg(self):
        assert self._dg is not None
        self._est_type = 0
        self._cDG = []
        for vi in xrange( self._curFiber.fiberLen-1 ):
            dg = np.zeros( (3,3) )
            for jj in xrange( self._curFiber.numFibers ):
                assert self._dg[ self._curFiber.fibers[jj][vi] ] is not None
                dg += self._dg[ self._curFiber.fibers[jj][vi] ]
            self._cDG.append( dg / float( self._curFiber.numFibers ) )

        assert len(self._cDG) == (self._curFiber.fiberLen-1)


    def save_dg(self, filename):
        # --------------------------------------------------
        # import matplotlib.pyplot as plt
        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        # ax.plot( range(len(self._cDG)), [dg[0,0] for dg in self._cDG], '.-' )
        # ax.plot( range(len(self._cDG)), [dg[1,1] for dg in self._cDG], '.-' )
        # ax.plot( range(len(self._cDG)), [dg[2,2] for dg in self._cDG], '.-' )
        # plt.show()
        # --------------------------------------------------
        with open(filename, 'w') as fout:
            for dg in self._cDG:
                print >> fout, ' '.join( [ '%.10f'%x for x in dg.reshape(9) ] )

    # Estimate the deformation gradient by formulating a LS problem
    def ls_estimate_dg(self):
        self._est_type = 1
        self._cDG = [ self._ls_estimate_dg_idx(i) for i in xrange( self._curFiber.fiberLen-1 ) ]
        #self._ls_estimate_dg_idx(75)

    def est_affine_dg(self):
        self._est_type = 2
        self._cDG = [ self._est_affine_dg(i) for i in xrange( self._curFiber.fiberLen-1 ) ]
        #self._est_affine_dg(75)

    def _est_affine_dg(self, vi):
        X0 = (self._initFiber.centralLine[:,vi] + self._initFiber.centralLine[:,vi+1])*0.5  # rest
        D0 = self._initFiber.centralLine[:,vi+1] - self._initFiber.centralLine[:,vi]
        L0 = linalg.norm(D0)
        assert L0 > 1E-6
        D0 /= L0

        # NOTE: d0 is not normalized
        x0 = (self._curFiber.centralLine[:,vi]  + self._curFiber.centralLine[:,vi+1])*0.5   # deformed
        d0 = self._curFiber.centralLine[:,vi+1] - self._curFiber.centralLine[:,vi]
        #l0 = linalg.norm(d0)
        assert linalg.norm(d0) > 1E-6

        # rotational matrix
        xdir = np.array([1., 0., 0.])
        quat = quaternion_from_two_vecs(D0, xdir)   # rotate back to x-dir

        A = np.zeros( (self._curFiber.numFibers * 3, 6) )
        b = np.empty( self._curFiber.numFibers * 3 )
        df1 = d0 / L0

        for i in xrange(self._curFiber.numFibers):
            XC = ( self._initFiber.fiber_vertex(i,vi) + self._initFiber.fiber_vertex(i,vi+1) )*0.5
            DX = quat.rotate( XC - X0 )             # rotate back to the initial frame

            xc = (self._curFiber.fiber_vertex(i,vi)+self._curFiber.fiber_vertex(i,vi+1))*0.5 - x0

            b[i*3:i*3+3] = xc - df1*DX[0]
            A[i*3:i*3+3,0:3] = np.identity(3)*DX[1]
            A[i*3:i*3+3,3:6] = np.identity(3)*DX[2]

        est = linalg.lstsq(A, b)[0]
        ret = np.empty( (3,3) )
        ret[:,0] = df1
        ret[:,1:3]= est.reshape( (3,2), order='F')
        return ret


    def _ls_estimate_dg_idx(self, vi):
        assert self._dg is not None
        X0 = (self._initFiber.centralLine[:,vi] + self._initFiber.centralLine[:,vi+1])*0.5  # rest
        x0 = (self._curFiber.centralLine[:,vi]  + self._curFiber.centralLine[:,vi+1])*0.5   # deformed
        XS = np.empty( (3,self._initFiber.numFibers) )
        LS = np.empty( self._initFiber.numFibers )
        n = 0
        for fi in xrange( self._initFiber.numFibers ):
            XS[:,fi] = ( self._initFiber.fiber_vertex(fi,vi) + self._initFiber.fiber_vertex(fi,vi+1) )*0.5
            LS[fi] = linalg.norm( XS[:,fi] - X0 )
            if LS[fi] > 1E-6: n += 1

        # construct the least-squares problem
        S = 9*self._curFiber.numFibers
        A = np.zeros( (S+n*3, 9) )
        b = np.empty( S+n*3 )
        idoff = S
        #np.set_printoptions(threshold=np.nan, suppress=True, precision=3)
        for i in xrange(self._curFiber.numFibers):
            A[i*9:i*9+9,:] = np.identity(9)
            b[i*9:i*9+9] = self._dg[ self._curFiber.fibers[i][vi] ].reshape(9)
            if LS[i] > 1E-6:
                DX = (XS[:,i] - X0).transpose() / LS[i] # <--------------------
                #assert np.abs(linalg.norm(DX)-1.) < 1E-10
                A[idoff,0:3]   = DX
                A[idoff+1,3:6] = DX
                A[idoff+2,6:9] = DX

                b[idoff:idoff+3] = ((self._curFiber.fiber_vertex(i,vi)+self._curFiber.fiber_vertex(i,vi+1))*0.5 - x0) / LS[i] # <----------------
                idoff += 3
        assert idoff == S+n*3

        #est = linalg.lstsq(A[:S,:], b[:S])[0]
        est = linalg.lstsq(A, b)[0]

        return est.reshape( (3,3) )

        # from mpl_toolkits.mplot3d import Axes3D
        # import matplotlib.pyplot as plt
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')

        # xx = b[S:].reshape( (3,-1), order='F' )

        # cc = A[S:,:].dot(est)
        # xxs = cc.reshape( (3,-1), order='F' )

        # ax.scatter(xxs[0,:], xxs[1,:], xxs[2,:], c='r', marker='o')
        # ax.scatter(xx[0,:], xx[1,:], xx[2,:], c='b', marker='^')
        # #for i in xrange(0,xxs.shape[1],2):
        # #    ax.plot( [xxs[0,i], xx[0,i]], [xxs[1,i], xx[1,i]], [xxs[2,i], xx[2,i]])
        # plt.show()


    def visualize(self, vi):
        if self._est_type is None:
            print 'Nothing to visualize'
            return

        if self._est_type < 2:
            X0 = (self._initFiber.centralLine[:,vi] + self._initFiber.centralLine[:,vi+1])*0.5
            x0 = (self._curFiber.centralLine[:,vi]  + self._curFiber.centralLine[:,vi+1])*0.5

            XS = np.empty( (3,self._initFiber.numFibers) )
            xs = np.empty( (3,self._curFiber.numFibers) )
            #XS2 = np.empty( (3,self._initFiber.numFibers) )
            for fi in xrange( self._initFiber.numFibers ):
                xt = ( self._initFiber.fiber_vertex(fi,vi) + self._initFiber.fiber_vertex(fi,vi+1) )*0.5 - X0
                LS = linalg.norm( xt )
                #XS[:,fi] = self._cDG[vi].dot( xt/LS )
                XS[:,fi] = self._cDG[vi].dot( xt )
                #XS2[:,fi] = self._cDG2[vi].dot( xt/LS )
                #XS[:,fi] = xt
                xs[:,fi] = (( self._curFiber.fiber_vertex(fi,vi)  + self._curFiber.fiber_vertex(fi,vi+1) )*0.5 - x0)

            from mpl_toolkits.mplot3d import Axes3D
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            ax.scatter(XS[0,:], XS[1,:], XS[2,:], c='r', marker='.')
            ax.scatter(xs[0,:], xs[1,:], xs[2,:], c='b', marker='.')
            #ax.scatter(XS2[0,:], XS2[1,:], XS2[2,:], c='g', marker='.')
            plt.show()

        elif self._est_type == 2:
            X0 = (self._initFiber.centralLine[:,vi] + self._initFiber.centralLine[:,vi+1])*0.5  # rest
            D0 = self._initFiber.centralLine[:,vi+1] - self._initFiber.centralLine[:,vi]
            L0 = linalg.norm(D0)
            D0 /= L0
            x0 = (self._curFiber.centralLine[:,vi]  + self._curFiber.centralLine[:,vi+1])*0.5

            XS = np.empty( (3,self._initFiber.numFibers) )
            xs = np.empty( (3,self._curFiber.numFibers) )

            xdir = np.array([1., 0., 0.])
            quat = quaternion_from_two_vecs(D0, xdir)   # rotate back to x-dir

            for fi in xrange( self._initFiber.numFibers ):
                xt = ( self._initFiber.fiber_vertex(fi,vi) + self._initFiber.fiber_vertex(fi,vi+1) )*0.5 - X0
                XS[:,fi] = self._cDG[vi].dot( quat.rotate(xt) )
                xs[:,fi] = (( self._curFiber.fiber_vertex(fi,vi)  + self._curFiber.fiber_vertex(fi,vi+1) )*0.5 - x0)

            from mpl_toolkits.mplot3d import Axes3D
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            ax.scatter(XS[0,:], XS[1,:], XS[2,:], c='r', marker='.')
            ax.scatter(xs[0,:], xs[1,:], xs[2,:], c='b', marker='.')
            plt.show()


# -------------------------------------------------------------------------------------------------------

def main(args):
    # initial fiber bundle
    initFiber = FiberBundle( '%s.obj' % args.infile[0] )
    initFiber.gen_central_line()
    #initFiber.output_obj('init.obj')

    # current fiber bundle
    curFiber = FiberBundle( '%s.obj' % args.infile[1] )
    curFiber.gen_central_line()

    initFiber.smooth_curve( 7,1 )
    curFiber.smooth_curve( 7,1 )

    est = DefGradEstimator( initFiber, curFiber, '%s.fe' % args.infile[1] )

    #est.average_dg()
    #est.ls_estimate_dg()
    est.est_affine_dg()
    est.visualize(75)
    #est.save_dg(args.infile[2])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Find deformation gradient for the central line of a yarn',
            epilog='Changxi Zheng (cxz@cs.columbia.edu)')
    parser.add_argument('infile', metavar='filename', type=str, nargs=3,
            help='<rest frame file>, <current frame file>, and <output filename>. The rest and current frame files are specified without the filename extensions.')

    args = parser.parse_args()
    main(args)
