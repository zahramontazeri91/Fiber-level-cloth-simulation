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


fig, ax = plt.subplots(2, sharex=True, figsize=(15,10))

lng = []
shrt = []
theta = []
rot = []
with open('../compress_new.txt', 'r') as fin:
    N = int(fin.readline())
    for i in range (0,N):
         compress = fin.readline().split()
         lng.append(float(compress[0]))
         shrt.append(float(compress[1]))
         theta.append(float(compress[2]))
         rot.append(float(compress[3]))
 
ind = np.arange(N)  

  
ax[0].plot([200,200], [0, 0.06], color='black') 
ax[0].plot([1200,1200], [0, 0.06], color='black') 
ax[0].set_ylim([-0.2, 2.5])
rects = ax[0].plot(ind, lng, color='r')
rects = ax[0].plot(ind, shrt, color='b')

            
ax[1].plot([200,200], [0, 3], color='black') 
ax[1].plot([1200,1200], [0, 3], color='black') 
rects = ax[1].plot(ind, theta, color='g')
rects = ax[1].plot(ind, rot, color='y')

plt.savefig("../../data/vis_crossSections/shapeMatch.png")

plt.show()