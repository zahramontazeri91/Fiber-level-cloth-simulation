# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 17:42:00 2018

@author: zahra
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D





def transform (dg, indx):
    #def transform(dg):
    fig = plt.figure(figsize=(10,10))
    #ax = fig.gca(projection='3d')
    ax = Axes3D(fig)
    
    # Cylinder
    z=np.linspace(-1, 1, 100)
    y=np.sqrt(1-z**2)
    y_=-1.0*np.sqrt(1-z**2)
    for i in range (-50,50):
        x=i/50.0
        # transfer x, y,z 
        x_rot=[]
        y_rot=[]
        z_rot=[]
        x_rot_=[]
        y_rot_=[]
        z_rot_=[]
        for j in range (0,100):
            
            vec = np.array([x,y[j],z[j]])

#            dg = np.identity(3)
#            a = np.array([0,1,0,0,0,1,1,0,0])
#            dg = a.reshape([3,3])
            vec_rot = np.matmul(dg, vec)
#            if j<3:
#                print(vec, vec_rot)
            x_rot.append(vec_rot[0])
            y_rot.append(vec_rot[1])
            z_rot.append(vec_rot[2])
            
            vec_ = np.array([x,y_[j],z[j]])
            vec_rot_ = np.matmul(dg, vec_)
            x_rot_.append(vec_rot_[0])
            y_rot_.append(vec_rot_[1])
            z_rot_.append(vec_rot_[2])

                
        ax.view_init(azim=90., elev=0)
        ax.set_xlim3d(-1.5,1.5)
        ax.set_ylim3d(-1.5,1.5)
        ax.set_zlim3d(-1.5,1.5)        
        ax.scatter(x_rot,y_rot,z_rot, color='b', alpha=0.3)
        ax.scatter(x_rot_,y_rot_,z_rot_, color='red', alpha=0.3)
    
    ax.set_xlabel("X", fontsize=30)
    ax.set_ylabel("Y", fontsize=30)
    ax.set_zlabel("Z", fontsize=30)
    
    plt.savefig('../../data/dg/dg_' + str(indx) + '.jpg')
#    plt.show()
    

fname = '../input/pattern/yarn4/spacing1.0x/00011/physical_17000_0.txt'
    
#fname = '../input/twist/yarn4/damp/physicalParam/physical_49500_0_world.txt'
data = np.loadtxt(fname)
N = data.shape[0]
#N=1
all_dg=[]
for d in range(0,N):
    trim = int(N*0.3)
    if (d>=trim and d<N-trim):
        dg=np.array([data[d][0], data[d][1], data[d][2],data[d][3], data[d][4], data[d][5],data[d][6], data[d][7], data[d][8] ] )
        dg_ = dg.reshape([3,3])
        all_dg.append(dg_)
        transform(dg_, d)