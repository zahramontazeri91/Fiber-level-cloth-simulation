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

#a = plt.subplot(111, aspect='equal')
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')

ax.add_patch(
    patches.Ellipse(
        (0.,0.),   # (x,y)
        0.08,          # width
        0.05,          # height
        alpha=0.5,
        facecolor = 'yellow'
    )
)
        
with open('../../data/pca_test.txt', 'r') as fin:
    num = int( fin.readline().split()[0] )
    for i in range (0,num):
        pnt = fin.readline().split()
        ax.scatter(float(pnt[0]), float(pnt[1]), alpha=0.5, color = 'blue', s = 100)
    
    num = int( fin.readline().split()[0] )
    for i in range (0,num):
        pnt = fin.readline().split()
        ax.scatter(float(pnt[0]), float(pnt[1]), alpha=0.5, color = 'red', s = 100)

    whitespace = fin.readline().split()
    cntr = fin.readline().split()
    longP = fin.readline().split()
    shortP = fin.readline().split()            
    ax.annotate("", xy=( float(longP[0]), float(longP[1]) ),
                xycoords='data', xytext=( float(cntr[0]),float(cntr[1]) ), 
                textcoords='data',
                arrowprops=dict(arrowstyle="->",
                color='black', 
                linewidth = 2, 
                connectionstyle="arc3"
                ),
    )
    ax.annotate("", xy=( float(shortP[0]), float(shortP[1]) ),
                xycoords='data', xytext=( float(cntr[0]),float(cntr[1]) ), 
                textcoords='data',
                arrowprops=dict(arrowstyle="->",
                color='g',
                linewidth = 2, 
                connectionstyle="arc3"
                ),
    )

    
    
    plt.tick_params(axis='both', which='major', labelsize=8)
    # set axes range
    plt.xlim(-.1,.1)
    plt.ylim(-.1,.1)
    plt.savefig("../../data/vis_crossSections/ints_plane%d.png" %i)
    plt.show() 