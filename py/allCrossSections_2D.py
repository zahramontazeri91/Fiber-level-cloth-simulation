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


def visIntersections(fname = ""):  
    with open(fname, 'r') as fin:
        plane_num = int( fin.readline().split()[1] )
        ply_num = int(fin.readline().split()[1])
    
        X = []
        Y = []
        for i in range(0, plane_num):
#            print ( "Display intersections with plane %d ... " %i)
            X.append([])
            Y.append([])
            fin.readline().split() 
            fin.readline().split() #read the plane index
            for p in range(0,ply_num):   
                X[-1].append([])
                Y[-1].append([])
                its_num = int(fin.readline().split()[1])
                plyCenter = fin.readline().split()
                plyCenter_x = float(plyCenter[1])
                plyCenter_y = float(plyCenter[2])
                X[-1][-1].append(plyCenter_x)
                Y[-1][-1].append(plyCenter_y)
                for j in range(0, its_num):
                    pos = [float(val) for val in fin.readline().strip().split(' ')]                
                    X[-1][-1].append(pos[0])
                    Y[-1][-1].append(pos[1])
    return X,Y
                    

 
                
                
#####################
fname = '../../data/allCrossSection2D_ref.txt'
outPath = "../../data/vis_crossSections/reference"
#print('Reading the first file ... ')
X1, Y1 = visIntersections(fname)

fname = '../../data/allCrossSection2D_deformedRef.txt'
outPath = "../../data/vis_crossSections/deformedRef"
#print('Reading the next file ... ')
X2, Y2 = visIntersections(fname)


fname = '../../data/allCrossSection2D_deformed.txt'
outPath = "../../data/vis_crossSections/deformed"
#print('Reading the next file ... ')
X3, Y3 = visIntersections(fname)


plane_num = 1184
for i in range(0, plane_num):  
    ply_num = 2
    plt.figure(figsize=(5,15))
    for p in range(0,ply_num):  
        its_num = 80
        for j in range(0, its_num):
            if j==10 and p==0:  
                z = 200
                c = 'red'
                s = 100
            elif j==20 and p==1: 
                z = 300
                c = 'black'
                s = 100
            else:
                z = 1
                c = "C" + str(p)
                s = 50

            plt.subplot(311)
            plt.scatter(X1[i][p][j], Y1[i][p][j], alpha=0.8, color = c, zorder = z, s = s)
            plt.tick_params(axis='both', which='major', labelsize=5)
            plt.xlim(-0.2,0.2)
            plt.ylim(-0.2,0.2)
            plt.title('Reference - %d' %i)
            
            plt.subplot(312)
            plt.scatter(X2[i][p][j], Y2[i][p][j], alpha=0.8, color = c, zorder = z, s = s)
            plt.tick_params(axis='both', which='major', labelsize=5)
            plt.xlim(-0.2,0.2)
            plt.ylim(-0.2,0.2)
            plt.title('Deformed reference - %d' %i)
            
            plt.subplot(313)
            plt.scatter(X3[i][p][j], Y3[i][p][j], alpha=0.8, color = c, zorder = z, s = s)
            plt.tick_params(axis='both', which='major', labelsize=5)
            plt.xlim(-0.2,0.2)
            plt.ylim(-0.2,0.2)
            plt.title('Simulated frame 9900 - %d' %i)
        
    
    plt.savefig("../../data/vis_crossSections/2Dcompare/plane%d.png" %i)
    plt.show()