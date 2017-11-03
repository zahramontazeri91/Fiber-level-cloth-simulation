import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})

fig = plt.figure(figsize=(10,20))
ax = fig.gca(projection='3d')
#ax.set_aspect('equal')

with open('../genYarn.txt', 'r') as fin:
    fiber_num = int(fin.readline())
    for i in range(0, fiber_num): 

        vrtx_num = int(fin.readline())
        #fiber = np.zeros((199, 3)) #between 700 to 900
        fiber = np.zeros((199, 3)) 
        cnt = 0
        for j in range(0, vrtx_num): 
#        for j in range(vrtx_num/4, vrtx_num*3/4):               
            pos = [float(val) for val in fin.readline().strip().split(' ')]
            
            if i==0 or i==80: #plot only ply-center
                if j>700 and j<900 :
                    fiber[cnt, :] = pos
                    cnt = cnt+1
            if i==0:
                c='red'
            if i==80:
                c='blue'

        ax.plot(fiber[:, 0], fiber[:, 1], fiber[:, 2], alpha=0.5, color =c)


#plt.tight_layout()
plt.xlim(-.05,.05)
plt.ylim(-0.05, 0.05)
#plt.zlim(-4.05, 4.05)
plt.savefig("../../data/vis_crossSections/29_plyCenter.png")
plt.show()
