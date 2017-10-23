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

import matplotlib.patches as mpatches

plyHelixRad = []
plyHelixTheta = []
plyHelixRad2 = []
plyHelixTheta2 = []
                
with open('../../data/plyCenters_simul.txt', 'r') as fin:
    plane_num = int( fin.readline().split()[1] )
    ply_num = int(fin.readline().split()[1])    
    
    with open('../../data/plyCenters_proc.txt', 'r') as fin2:
        plane_num = int( fin2.readline().split()[1] )
        ply_num = int(fin2.readline().split()[1])
    
        
        with open('../../data/orientation.txt','r') as fv: 
            for i in range(0, plane_num):
#                print ( "Display intersections with plane %d ... " %i)
                whitespace = fin.readline().split()
                whitespace = fin2.readline().split()
        
                X = []
                Y = []
               #plt.scatter(0.0,0.0, alpha=1.0, color = 'red', s = 100) #Because center is always in the middle
                
                plyCenterX = []
                plyCenterY = []
                plyCenterX2 = []
                plyCenterY2 = []
                            


                for p in range(0,ply_num):    
    #                its_num = int(fin.readline().split()[1])
                    X = []
                    Y = []
                    
                    c = "C" + str(p)
                    plyCenter = fin.readline().split() 
                    plyCenterX.append(plyCenter[1])
                    plyCenterY.append(plyCenter[2])
                    #plt.scatter(plyCenter[1], plyCenter[2], alpha=0.8, color = 'green', zorder=200)
                    
                    plyHelix = fin.readline().split()
                    if p==0 : #plot R only for one ply                    
                        plyHelixRad.append(float(plyHelix[1]))
                        plyHelixTheta.append(float(plyHelix[2]))
                    
                  
                    # procedural param
                    plyCenter = fin2.readline().split() 
                    plyCenterX2.append(plyCenter[1])
                    plyCenterY2.append(plyCenter[2])
#                    plt.scatter(plyCenter[1], plyCenter[2], alpha=0.8, color = 'red', zorder=200)
                    
                    plyHelix2 = fin2.readline().split()
                    if p==0 :
                        plyHelixRad2.append(float(plyHelix2[1]))
                        plyHelixTheta2.append(float(plyHelix2[2]))
                    

             
#            print(plyHelixRad)
            #plt.plot(plyCenterX, plyCenterY, linewidth=2, color='green')               
            #plt.plot(plyCenterX2, plyCenterY2, linewidth=2, color='red')
            plt.figure(1)
            planeId = np.linspace(0, plane_num, num=plane_num )
            plt.plot(planeId, plyHelixRad, 'g', planeId, plyHelixRad2, 'r')            
            red_patch = mpatches.Patch(color='red', label='Procedural helix-radius')
            green_patch = mpatches.Patch(color='green', label='Simulated helix-radius')
            plt.legend(handles=[red_patch, green_patch],prop={'size': 6}) 
            plt.tick_params(axis='both', which='major', labelsize=8)
#            plt.xlim(-.05,.05)
            plt.ylim(0.0, 0.055)
            plt.savefig("../../data/vis_crossSections/helix_Rad.png", dpi=1000)
            plt.show() 

            
            planeId = np.linspace(0, plane_num, num=plane_num )
            plt.plot(planeId, plyHelixTheta, 'g', planeId, plyHelixTheta2, 'r')
            red_patch = mpatches.Patch(color='red', label='Procedural helix-theta')
            green_patch = mpatches.Patch(color='green', label='Simulated helix-theta')
            plt.legend(handles=[red_patch, green_patch],prop={'size': 6}) 
            plt.tick_params(axis='both', which='major', labelsize=8)
            plt.savefig("../../data/vis_crossSections/helix_theta.png", dpi=1000)
            plt.show() 