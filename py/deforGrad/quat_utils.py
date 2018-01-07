#! /usr/bin/env python

import numpy as np
from pyquaternion import Quaternion

def quaternion_from_two_vecs(u, v):
    sqlength = lambda x: np.inner(x, x)

    norm_u_norm_v = np.sqrt(sqlength(u) * sqlength(v))
    w = np.cross(u, v)
    q = Quaternion(scalar=norm_u_norm_v + np.inner(u, v), vector=w)
    return q.normalised

if __name__ == '__main__':
    import numpy.linalg as linalg
    u = np.array([1,9,2])
    u = u / np.linalg.norm(u)
    v = np.array([0,1,2])
    v = v / np.linalg.norm(v)

    quat = quaternion_from_two_vecs( u, v )
    print quat.rotate(u), v
    print type(quat.rotate(u))
