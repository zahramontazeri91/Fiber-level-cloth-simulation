# -*- coding: utf-8 -*-
"""
Visualize yarn intersection with given planes
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})
  
with open('../../data/allCrossSection2D_test.txt', 'r') as fin:
    plane_num = int( fin.readline().split()[1] )
    ply_num = int(fin.readline().split()[1])
    
    with open('../../data/orientation.txt','r') as fv: 
        for i in range(0, 10):
            print ( "Display intersections with plane %d ... " %i)
            whitespace = fin.readline().split()
            #center = fin.readline().split()  #Because center is always in the middle
            X = []
            Y = []
            plt.scatter(0.0,0.0, alpha=1.0, color = 'red', s = 100)
            
            for p in range(0,ply_num):
                its_num = int(fin.readline().split()[1])
                X = []
                Y = []
                c = "C" + str(p)
                for j in range(0, its_num):
                    pos = [float(val) for val in fin.readline().strip().split(' ')]                
                    X.append(pos[0])
                    Y.append(pos[1])
                plt.scatter(X, Y, alpha=0.5, color = c)
            
            # draw the eigen vectors
            cntr = fv.readline().split()
            longP = fv.readline().split()
            shortP = fv.readline().split()
            whitespace = fv.readline().split()
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