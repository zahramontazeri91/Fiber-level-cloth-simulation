import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})

fig = plt.figure(figsize=(10,20))
ax = fig.gca(projection='3d')
#ax.set_aspect('equal')

with open('../R_theta.txt', 'r') as fin:
    x=[]
    y=[]
    z=[]
    for i in range (0,1526):
#        x = float( fin.readline().split()[0] )
#        y = float( fin.readline().split()[1] )
#        plt.scatter(i, theta, alpha=0.8, color = 'r', s=2)
        pos = [float(val) for val in fin.readline().strip().split(' ')]
        x.append(pos[0])
        y.append(pos[1])
#        z.append(pos[2])
        pos = [float(val) for val in fin.readline().strip().split(' ')]
   
    z=np.linspace(0,1526,num=1526)
    ax.plot(x,y,z, alpha=0.5, color ='r')
        
    plt.show() 
