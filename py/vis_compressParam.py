"""
visualize compress parameter

@author: zahra
"""
import numpy as np
import matplotlib

from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import math


fig, ax = plt.subplots(figsize=(15,5))

lng = []
shrt = []
theta = []
with open('../compress.txt', 'r') as fin:
    N = int(fin.readline())
    for i in range (0,N):
         compress = fin.readline().split()
         lng.append(float(compress[0]))
         shrt.append(float(compress[1]))
         theta.append(float(compress[2]))
 
ind = np.arange(N)    
ax.plot([200,200], [0, 0.06], color='black') 
ax.plot([1200,1200], [0, 0.06], color='black') 
rects = ax.plot(ind, lng, color='r')
rects = ax.plot(ind, shrt, color='b')
plt.savefig("../../data/vis_crossSections/ellipseShape.png")

#ax.plot([200,200], [0, 3], color='black') 
#ax.plot([1200,1200], [0, 3], color='black') 
#rects = ax.plot(ind, theta, color='g')
#plt.savefig("../../data/vis_crossSections/ellipseTheta.png")

plt.show()