# -*- coding: utf-8 -*-
"""
test pca

"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})
  
with open('../../data/pca_test.txt', 'r') as fin:
    num = int( fin.readline().split()[0] )
    for i in range (0,num):
        pnt = fin.readline().split()
        plt.scatter(float(pnt[0]), float(pnt[1]), alpha=0.5, color = 'blue', s = 100)
    
    num = int( fin.readline().split()[0] )
    for i in range (0,num):
        pnt = fin.readline().split()
        plt.scatter(float(pnt[0]), float(pnt[1]), alpha=0.5, color = 'red', s = 100)

    whitespace = fin.readline().split()
    cntr = fin.readline().split()
    longP = fin.readline().split()
    shortP = fin.readline().split()            
    plt.annotate("", xy=( float(longP[0]), float(longP[1]) ),
                xycoords='data', xytext=( float(cntr[0]),float(cntr[1]) ), 
                textcoords='data',
                arrowprops=dict(arrowstyle="->",
                facecolor='red', 
                connectionstyle="arc3"
                ),
    )
    plt.annotate("", xy=( float(shortP[0]), float(shortP[1]) ),
                xycoords='data', xytext=( float(cntr[0]),float(cntr[1]) ), 
                textcoords='data',
                arrowprops=dict(arrowstyle="->",
                facecolor='blue', 
                connectionstyle="arc3"
                ),
    )
        
    plt.tick_params(axis='both', which='major', labelsize=8)
    # set axes range
    plt.xlim(-.1,.1)
    plt.ylim(-.1,.1)
    plt.savefig("../../data/vis_crossSections/ints_plane%d.png" %i)
    plt.show() 