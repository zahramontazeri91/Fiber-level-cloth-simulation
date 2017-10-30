import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})

fig = plt.figure(figsize=(10,20))
ax = fig.gca(projection='3d')
#ax.set_aspect('equal')

#with open('../../data/plyCenter_proc.txt', 'r') as fin:
with open('../plyCenter.txt', 'r') as fin:
    x1=[]
    y1=[]
    x2=[]
    y2=[]

    cnt = 0
    for i in range (0,1526):
#        x = float( fin.readline().split()[0] )
#        y = float( fin.readline().split()[1] )
#        plt.scatter(i, theta, alpha=0.8, color = 'r', s=2)
        pos = [float(val) for val in fin.readline().strip().split(' ')]        
        if i>700 and i<900 :
            x1.append(pos[0])
            y1.append(pos[1])

            
        pos = [float(val) for val in fin.readline().strip().split(' ')] #skip ply-center2
        if i>700 and i<900 :
            x2.append(pos[0])
            y2.append(pos[1])
            cnt = cnt+1
        whitespace = fin.readline()

   
    z=np.linspace(0,cnt,num=cnt)
    ax.plot(x1,y1,z, alpha=0.5, color ='r')
    ax.plot(x2,y2,z, alpha=0.5, color ='b')
        
    plt.savefig("../../data/vis_crossSections/frame29_plyCenter.png")
    plt.show() 
