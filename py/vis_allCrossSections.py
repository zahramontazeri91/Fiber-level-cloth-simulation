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
        
with open('../allCrossSection.txt', 'r') as fin:
    plane_num = int( fin.readline().split()[1] )
    ply_num = int(fin.readline().split()[1])

    for i in range(0, plane_num):
        print ( "Display intersections with plane %d ... " %i)
        whitespace = fin.readline().split()
        center = fin.readline().split()
        ax.scatter(float(center[1]), float(center[2]), float(center[3]), s=50, color='red')
           
        for p in range(0,ply_num):
            its_num = int(fin.readline().split()[1])
            c = "C" + str(p)
            for j in range(0, its_num):
                pos = [float(val) for val in fin.readline().strip().split(' ')]
                ax.scatter(pos[0], pos[1], pos[2], alpha=0.5, color = c )      

    ax.view_init(azim=10., elev=10)
    plt.savefig("../../data/vis_crossSections/ints_plane%d.png" %i)

        



plt.tight_layout()
plt.show()


