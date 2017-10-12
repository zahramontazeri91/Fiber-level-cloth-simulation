# -*- coding: utf-8 -*-
"""
Visualize yarn intersection with given planes
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})
      
with open('../../data/allCrossSection.txt', 'r') as fin:
    plane_num = int( fin.readline().split()[1] )
    ply_num = int(fin.readline().split()[1])
    

    for i in range(0, plane_num):
        print ( "Display intersections with plane %d ... " %i)
        whitespace = fin.readline().split()
        center = fin.readline().split()
        X = []
        Y = []
        cx = []
        cy = []
        cx.append(float(center[1]))
        cy.append(float(center[2]))
        plt.scatter(cx,cy, alpha=1.0, color = 'red', s = 100)
        
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
        plt.tick_params(axis='both', which='major', labelsize=8)
        plt.savefig("../../data/vis_crossSections/ints_plane%d.png" %i)
        plt.show()  
        

          

   

