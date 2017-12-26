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

# In[]:
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
                    

 
                
# In[]:
def plotIntersects(X1, Y1, X2, Y2, X3, Y3, scale, plane_num, ply_num, its_num, title1, title2, title3):
    for i in range(0, plane_num):  
        plt.figure(figsize=(5,15))
        for p in range(0,ply_num):  
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
                plt.scatter(X1[i][p][j]*scale, Y1[i][p][j]*scale, alpha=0.8, color = c, zorder = z, s = s)
                plt.tick_params(axis='both', which='major', labelsize=5)
                plt.xlim(-0.1,0.1)
                plt.ylim(-0.1,0.1)
                plt.title(title1 + ' - %d' %i)
                
                plt.subplot(312)
                plt.scatter(X2[i][p][j]*scale, Y2[i][p][j]*scale, alpha=0.8, color = c, zorder = z, s = s)
                plt.tick_params(axis='both', which='major', labelsize=5)
                plt.xlim(-0.1,0.1)
                plt.ylim(-0.1,0.1)
                plt.title(title2 + ' - %d' %i)
                
                plt.subplot(313)
                plt.scatter(X3[i][p][j]*scale, Y3[i][p][j]*scale, alpha=0.8, color = c, zorder = z, s = s)
                plt.tick_params(axis='both', which='major', labelsize=5)
                plt.xlim(-0.1,0.1)
                plt.ylim(-0.1,0.1)
                plt.title(title3 + ' - %d' %i)
            
        
        plt.savefig("../../data/vis_crossSections/2Dcompare/plane%d.png" %i)
        plt.show()               
# In[]:
# plot procedural data
fname = '../../data/allCrossSection2D_simulate.txt'
X1, Y1 = visIntersections(fname)
fname = '../../data/allCrossSection2D_ref.txt'
#fname = '../../data/allCrossSection2D_compress.txt'
X2, Y2 = visIntersections(fname)
fname = '../../data/allCrossSection2D_curve.txt'
X3, Y3 = visIntersections(fname)
scale = 1
title1 = 'Simulate'
title2 = 'Compression'
title3 = 'Curve-mapping'
plane_num = 120  
ply_num = 2
its_num = 80
plotIntersects(X1, Y1, X2, Y2, X3, Y3, scale, plane_num, ply_num, its_num, title1, title2, title3)

# In[]:
## plot simulated data
#fname = '../../data/allCrossSection2D_ref.txt'
#X1, Y1 = visIntersections(fname)
#fname = '../../data/allCrossSection2D_deformedRef.txt'
#X2, Y2 = visIntersections(fname)
#fname = '../../data/allCrossSection2D_deformed.txt'
#X3, Y3 = visIntersections(fname)
#scale = 1
#title1 = 'Reference'
#title2 = 'Deformed-Ref'
#title3 = 'Simulated'
#plane_num = 300  
#ply_num = 2
#its_num = 80
#plotIntersects(X1, Y1, X2, Y2, X3, Y3, scale, plane_num, ply_num, its_num, title1, title2, title3)