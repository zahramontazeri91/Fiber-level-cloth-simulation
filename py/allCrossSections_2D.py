# -*- coding: utf-8 -*-
"""
Visualize yarn intersection with given planes
"""
import numpy as np
import matplotlib

from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import math


    
with open('../../data/allCrossSection2D.txt', 'r') as fin:
    plane_num = int( fin.readline().split()[1] )
    ply_num = int(fin.readline().split()[1])
    
    with open('../../data/orientation.txt','r') as fv: 
        for i in range(0, plane_num):
            print ( "Display intersections with plane %d ... " %i)
            whitespace = fin.readline().split()
            #center = fin.readline().split()  #Because center is always in the middle
            X = []
            Y = []
            plt.scatter(0.0,0.0, alpha=1.0, color = 'red', s = 100)
            
            # draw ellipse
            cntr = fv.readline().split()
            param = fv.readline().split()
            whitespace = fv.readline().split()
            width = 2.0 * float(param[0])
            height = 2.0 * float(param[1])
            angle = math.degrees(float(param[2]))
            ax = plt.subplot(111, aspect='equal')
            ell = Ellipse((float(cntr[0]),float(cntr[1])), 
                          width, height, alpha=0.3, 
                          angle = angle, facecolor = 'yellow' )  
            ax.add_artist(ell)
            
            for p in range(0,ply_num):    
                its_num = int(fin.readline().split()[1])
                X = []
                Y = []
                c = "C" + str(p)
                for j in range(0, its_num):
                    pos = [float(val) for val in fin.readline().strip().split(' ')]                
                    X.append(pos[0])
                    Y.append(pos[1])
                plt.scatter(X, Y, alpha=0.8, color = c)
            
            
            plt.tick_params(axis='both', which='major', labelsize=8)
            # set axes range
            plt.xlim(-.1,.1)
            plt.ylim(-.1,.1)
            plt.savefig("../../data/vis_crossSections/ints_plane%d.png" %i)
            plt.show()  