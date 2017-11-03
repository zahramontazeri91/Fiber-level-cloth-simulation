"""
Visualize 3d points in this format:
point_number
x y z
...
"""

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 16})

fig = plt.figure(figsize=(10,20))
ax = fig.gca(projection='3d')
#ax.set_aspect('equal')
        
with open('../../data/test_plyCenter.txt', 'r') as fin:
    n = int(fin.readline())  
        
    for j in range(0, n):
        pos = [float(val) for val in fin.readline().strip().split(' ')]
        
        ax.scatter(pos[0], pos[1], pos[2], alpha=0.5, color='blue')      
    # rotate the axes and update
#    for ii in xrange(0,90,10):
#            ax.view_init(azim=10., elev=ii)
#            plt.savefig("../../data/vis_crossSections/view%d.png" %  ii )
        

ax.view_init(azim=0., elev=0)
ax.set_xlim3d(-0.1,0.1)
ax.set_ylim3d(-0.1,0.1)
ax.set_zlim3d(-8,8)
#plt.tight_layout()
plt.show()