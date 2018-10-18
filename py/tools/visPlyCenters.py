# -*- coding: utf-8 -*-
"""
Visualize yarn intersection with given planes
and draw the bounding ellipse
"""
import numpy as np
import matplotlib

from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import math

import matplotlib.patches as mpatches

with open('../../data/plyCenters_simul.txt', 'r') as fin:
    plane_num = int( fin.readline().split()[1] )
    ply_num = int(fin.readline().split()[1])    
    
    with open('../../data/plyCenters_proc.txt', 'r') as fin2:
        plane_num = int( fin2.readline().split()[1] )
        ply_num = int(fin2.readline().split()[1])
    
        
        with open('../../data/orientation.txt','r') as fv: 
            for i in range(0, plane_num):
                print ( "Display intersections with plane %d ... " %i)
                whitespace = fin.readline().split()
                whitespace = fin2.readline().split()
        
                X = []
                Y = []
                plt.scatter(0.0,0.0, alpha=1.0, color = 'red', s = 100) #Because center is always in the middle
                
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
                
                plyCenterX = []
                plyCenterY = []
                plyCenterX2 = []
                plyCenterY2 = []
                            
                plyHelixRad = []
                plyHelixTheta = []
                plyHelixRad2 = []
                plyHelixTheta2 = []

                for p in range(0,ply_num):    
    #                its_num = int(fin.readline().split()[1])
                    X = []
                    Y = []
                    
                    c = "C" + str(p)
                    plyCenter = fin.readline().split() 
                    plyCenterX.append(plyCenter[1])
                    plyCenterY.append(plyCenter[2])
                    plt.scatter(plyCenter[1], plyCenter[2], alpha=0.8, color = 'green', zorder=200)
                    
                    plyHelix = fin.readline().split()
                    plyHelixRad.append(plyHelix[1])
                    plyHelixTheta.append(plyHelix[2])
                    
                    plyCenter = fin2.readline().split() 
                    plyCenterX2.append(plyCenter[1])
                    plyCenterY2.append(plyCenter[2])
                    plt.scatter(plyCenter[1], plyCenter[2], alpha=0.8, color = 'red', zorder=200)
                    
                    plyHelix2 = fin2.readline().split()
                    plyHelixRad2.append(plyHelix2[1])
                    plyHelixTheta2.append(plyHelix2[2])
                    
                plt.plot(plyCenterX, plyCenterY, linewidth=2, color='green')               
                plt.plot(plyCenterX2, plyCenterY2, linewidth=2, color='red')
    
        
                    
                red_patch = mpatches.Patch(color='red', label='Procedural ply-centers')
                green_patch = mpatches.Patch(color='green', label='Simulated ply-centers')
                plt.legend(handles=[red_patch, green_patch],prop={'size': 6}) 
                plt.tick_params(axis='both', which='major', labelsize=8)
                # set axes range
                plt.xlim(-.05,.05)
                plt.ylim(-.05,.05)
                plt.savefig("../../data/vis_crossSections/ints_plane%d.png" %i)
                plt.show()  