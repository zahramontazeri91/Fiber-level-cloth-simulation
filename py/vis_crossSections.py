# -*- coding: utf-8 -*-
"""
Visualize yarn intersection with given planes
Format: 
number of fiber intersections
x y z
...
"""

import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 16})
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect('equal')

with open('../../data/crossSection.txt', 'r') as fin:
    n = int(fin.readline())
    center = [float(val) for val in fin.readline().strip().split(' ')]
    ax.scatter(center[0], center[1], center[2], color='red')
    #print (center)
    
    for i in range(0, n):
        pos = [float(val) for val in fin.readline().strip().split(' ')]
        ax.scatter(pos[0], pos[1], pos[2], alpha=0.5, color='blue')      

# rotate the axes and update
for ii in xrange(0,90,10):
        ax.view_init(azim=10., elev=ii)
        plt.savefig("../../data/vis_crossSections/view%d.png" % ii)


plt.tight_layout()
plt.show()