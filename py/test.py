# -*- coding: utf-8 -*-
"""
Created on Wed Nov 01 18:39:32 2017

@author: zahra
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})
fig = plt.figure(figsize=(10,20))
ax = fig.gca(projection='3d')
ax.set_aspect('equal', 'datalim')
        
with open('../../data/allCrossSection.txt', 'r') as fin:
    plane_num = int( fin.readline().split()[1] )
    ply_num = int(fin.readline().split()[1])
    
    for i in range(0, plane_num):

#        matplotlib.rcParams.update({'font.size': 16})   
#        fig = plt.figure(figsize=(10,20))
#        ax = fig.gca(projection='3d')
#        ax.set_aspect('equal', 'datalim')

        print ( "Display intersections with plane %d ... " %i)
        whitespace = fin.readline().split()
           
        for p in range(0,ply_num):
            its_num = int(fin.readline().split()[1])
            center = fin.readline().split()
            ax.scatter(float(center[1]), float(center[2]), float(center[3]), s=50, color='red')

            c = "C" + str(p)
            for j in range(0, its_num):
                pos = [float(val) for val in fin.readline().strip().split(' ')]
                ax.scatter(pos[0], pos[1], pos[2], alpha=0.5,  s=10, color = c )      
        ax.set_xlim3d(-0.1,0.1)
        ax.set_ylim3d(-0.1,0.1)
        ax.set_zlim3d(-8,8)
        ax.view_init(azim=0., elev=0)
        # Turn off tick labels
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
#        ax.set_zlabel('Z axis')

#        plt.savefig("../../data/vis_crossSections/ints_plane%d.png" %i)
#        plt.show()
 
    plt.savefig("../../data/vis_crossSections/ints_plane%d.png" %i)
    plt.show()  