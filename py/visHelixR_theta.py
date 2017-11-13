# -*- coding: utf-8 -*-
"""
Visualize 2D plot of plyCenters rotation and theta
"""
import numpy as np
import matplotlib

from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import math


fig, ax = plt.subplots(figsize=(15,5))

R = []
theta = []
with open('../parameterziePlyCntr.txt', 'r') as fin:
    N = int(fin.readline())
    print(N)
    for i in range (0,N):
         helix = fin.readline().split()
         R.append(float(helix[0]))
         theta.append(float(helix[1]))

ind = np.arange(N)  

ax.plot([200,200], [0, 0.05], color='black')
ax.plot([1200,1200], [0, 0.05], color='black')   
rects = ax.plot(ind, R, color='r')
plt.savefig("../../data/vis_crossSections/helixR.png")

#ax.plot([200,200], [-3,3], color='black')
#ax.plot([1200,1200], [-3,3], color='black')
#rects = ax.plot(ind, theta, color='b')
#plt.savefig("../../data/vis_crossSections/helixTheta.png")

plt.show()