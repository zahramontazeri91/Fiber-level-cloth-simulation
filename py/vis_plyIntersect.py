# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 09:05:44 2017

@author: zahra
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
        X = []
        Y = []
        its_num = int(fin.readline().split()[0])
        for j in range(0, its_num):
            pos = [float(val) for val in fin.readline().strip().split(' ')]                
            X.append(pos[0])
            Y.append(pos[1])
    return X,Y

# In[]:
def plotIntersects(X1, Y1, X2, Y2, its_num, title1):
    plt.figure(figsize=(5,5))
    for j in range(0, its_num):
        
        scale = 1
        s = 50
        plt.scatter(X1[j]*scale, Y1[j]*scale, alpha=1.0, color = 'g', zorder = 100, s = s, label ='Reference')
        plt.tick_params(axis='both', which='major', labelsize=5)
        plt.xlim(-0.1,0.1)
        plt.ylim(-0.1,0.1)
        
        plt.scatter(X2[j]*scale, Y2[j]*scale, alpha=0.3, color = 'r', zorder = 500, s = s, label ='Simulated')
        plt.tick_params(axis='both', which='major', labelsize=5)
        plt.xlim(-0.1,0.1)
        plt.ylim(-0.1,0.1)
        plt.title(title1 )
    
    plt.savefig("../../data/vis_crossSections/2Dcompare/ply.png")
    plt.show() 
        
        
# In[]:
# plot one ply
fname = '../../data/plyShapeMatch_simul.txt'
X1, Y1 = visIntersections(fname)
fname = '../../data/plyShapeMatch_proc.txt'
X2, Y2 = visIntersections(fname)

title1 = 'Simulated vs Reference - 100 - perPly'
its_num = 160
plotIntersects(X1, Y1, X2, Y2, its_num, title1)