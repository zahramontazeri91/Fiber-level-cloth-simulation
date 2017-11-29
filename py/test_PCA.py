# -*- coding: utf-8 -*-
"""
test pca

"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})
import matplotlib.patches as patches
import math

#a = plt.subplot(111, aspect='equal')



        
with open('../../data/pca_test.txt', 'r') as fin:

    N = int(fin.readline())
    for i in range (0,N):
        print (i)
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')        
#        whitespace = fin.readline().split()
        cntr = [0,0]
        param = fin.readline().split()
        old_x = float(param[0])
        old_y = float(param[1]) 
        longP_x = float(param[2])
        longP_y = float(param[3])
        shortP_x = float(param[4])
        shortP_y = float(param[5])
      
#        ax.add_patch(patches.Ellipse(
#                    (cntr[0], cntr[1]),   # (x,y)
#                    2.0*longP,          # width
#                    2.0*shortP,          # height
#                    alpha=0.5,
#                    angle = angle ,
#                    facecolor = 'yellow') )
           
        ax.annotate("", xy=(longP_x, longP_y), xytext=(0,0), 
                    textcoords='data',
                    arrowprops=dict(arrowstyle="->",
                    color='blue', 
                    linewidth = 2, 
                    connectionstyle="arc3"
                    ),
        )
        ax.annotate("", xy=(shortP_x, shortP_y), xytext=(0,0), 
                    textcoords='data',
                    arrowprops=dict(arrowstyle="->",
                    color='red',
                    linewidth = 2, 
                    connectionstyle="arc3"
                    ),
        )
    
        ax.annotate("", xy=(old_x, old_y), xytext=(0,0), 
                    textcoords='data',
                    arrowprops=dict(arrowstyle="->",
                    color='black',
                    linewidth = 2, 
                    connectionstyle="arc3"
                    ),
        )
        
        plt.tick_params(axis='both', which='major', labelsize=8)
        # set axes range
        plt.xlim(-1,1)
        plt.ylim(-1,1)
        plt.savefig("../../data/vis_crossSections/ints_plane%d.png" %i)
        plt.show() 