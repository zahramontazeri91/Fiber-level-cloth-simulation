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




lng = []
shrt = []
theta = []
rot = []
with open('../compressParams.txt', 'r') as fin:
    N = int(fin.readline())
    for i in range (0,N):
         compress = fin.readline().split()
         lng.append(float(compress[0]))
         shrt.append(float(compress[1]))
         theta.append(float(compress[2]))
         rot.append(float(compress[3]))
 
ind = np.arange(N)  

#f, axs = plt.subplots(2,1,figsize=(15,15))
#ax = plt.subplot(1,2,1)
plt.figure(figsize=(15,10))

plt.subplot(311)
plt.plot([200,200], [0, 0.06], color='black') 
plt.plot([1200,1200], [0, 0.06], color='black') 
#plt.ylim(-0.2, 2.5)
plt.ylim(0.0, 2.5)
plt.plot(ind, lng, color='r', label='Sx')
plt.plot(ind, shrt, color='b', label='Sx')
plt.legend()
plt.title('Simulated data frame 29')
          
plt.subplot(312)         
plt.plot([200,200], [0, 3], color='black') 
plt.plot([1200,1200], [0, 3], color='black') 
#plt.plot(ind, rot, color='y', label=r'$\theta_R$')
plt.plot(ind, theta, color='g', label=r'$\theta_S$')
plt.legend()

plt.subplot(313)         
plt.plot([200,200], [0, 3], color='black') 
plt.plot([1200,1200], [0, 3], color='black') 
plt.plot(ind, rot, color='y', label=r'$\theta_R$')
#plt.plot(ind, theta, color='g', label=r'$\theta_S$')
plt.legend()

plt.savefig("../../data/vis_crossSections/shapeMatch.png")
plt.show()