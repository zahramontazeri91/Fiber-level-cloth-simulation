# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 15:32:40 2018

@author: zahra
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import hsv_to_rgb


dg = []
fname = '../input/pattern/yarn4/spacing1.0x/00011/physical_17000_0.txt'
data = np.loadtxt(fname)
N = data.shape[0]
trim = int(N*0.15)

with open(fname, "r") as fin:
    d =0 
    for line in fin.readlines():
        v = [float(x) for x in line.strip().split()]
        assert len(v) == 9
        v = np.reshape(v, (3, 3))
        if (d>=trim and d<N-trim):
            dg.append(v)
        d=d+1
#assert len(dg) == len(vtxL) - 1

m = len(dg)
n = 100
theta = np.linspace(0, 2.0*np.pi, n + 1)
theta = theta[:-1]
y = np.cos(theta)
z = np.sin(theta)

fig = plt.figure(figsize=(10,10))

for r in range(0, m):
    fig.clf()
    ax = fig.gca(projection='3d')

    p0 = np.dot(dg[r], [-1.5, 0.0, 0.0])
    p1 = np.dot(dg[r], [1.5, 0.0, 0.0])
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]], color='gray', linewidth=3, alpha=0.5)
    for i in range(0, n):
        clr = hsv_to_rgb([float(i)/(n + 1), 1.0, 1.0])

        p0 = np.dot(dg[r], [-1.0, y[i], z[i]])
        p1 = np.dot(dg[r], [1.0, y[i], z[i]])
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]], color=clr, linewidth=2, alpha=0.5)

    # Create cubic bounding box to simulate equal aspect ratio
    max_range = 4.0
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten()
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten()
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten()
    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
       ax.plot([xb], [yb], [zb], 'w')


    ax.view_init(azim=90, elev=0)
    ax.set_xlim3d(-0.5*max_range, 0.5*max_range)
    ax.set_xlabel('X', fontsize=30)
    ax.set_ylim3d(-0.5*max_range, 0.5*max_range)
    ax.set_ylabel('Y', fontsize=30)
    ax.set_zlim3d(-0.5*max_range, 0.5*max_range)
    ax.set_zlabel('Z', fontsize=30)

#    plt.show()
    plt.savefig("../../data/cylinder/world_%04d.png" % r, dpi=80)