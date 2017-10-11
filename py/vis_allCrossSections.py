# -*- coding: utf-8 -*-
"""
Visualize yarn intersection with given planes
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect('equal')
        
#with open('../allCrossSection.txt', 'r') as fin:
with open('../test_planes.txt', 'r') as fin:
    plane_num = int(fin.readline())
    for i in range(0, plane_num):


        its_num = int(fin.readline())
        center = [float(val) for val in fin.readline().strip().split(' ')]
        ax.scatter(center[0], center[1], center[2], s=50, color='red')
        
        for j in range(0, its_num):
            pos = [float(val) for val in fin.readline().strip().split(' ')]
            ax.scatter(pos[0], pos[1], pos[2], alpha=0.5, color='blue')      
    # rotate the axes and update
    for ii in xrange(0,90,10):
            ax.view_init(azim=10., elev=ii)
            plt.savefig("vis_cross/plane%d_view%d.png" % (i, ii) )
        



plt.tight_layout()
plt.show()


