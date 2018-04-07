# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 16:41:42 2018

@author: zahra
"""

import numpy as np
import matplotlib.pyplot as plt
    

plt.figure(figsize=(9,9))
skipfactor = 500
c = 1
for f in range (16000/skipfactor, 16000/skipfactor + 1):
    frame_no = f*skipfactor
    fname = '../input/pattern/yarn4/spacing1.5x/10100/physical_' + str(frame_no) + '_0.txt'
#    fname = '../input/woven/yarn4/spacing1.0x/00011/physical_' + str(frame_no) + '_0.txt'
#    fname = '../input/twist/yarn4/damp/physical_' + str(frame_no) + '_0.txt'
    
    data = np.loadtxt(fname)
    N = data.shape[0]
    dg_size = data.shape[1]
#    dg0=[]
#    dg1=[]
#    dg2=[]
#    dg3=[]
#    dg4=[]
#    dg5=[]
#    dg6=[]
#    dg7=[]
#    dg8=[]
#    for i in range (0,N):
#        dg0.append(data[i,0])
#        dg1.append(data[i,1])
#        dg2.append(data[i,2])
#        
#        dg3.append(data[i,3])
#        dg4.append(data[i,4])
#        dg5.append(data[i,5])
#        
#        dg6.append(data[i,6])
#        dg7.append(data[i,7])
#        dg8.append(data[i,8])
#    
#    
#    
#    plt.subplot(3,3,1)
#    plt.plot(dg0, label='dg0')
#    plt.ylim([-1.5,1.5])
#    
#    plt.subplot(3,3,2)
#    plt.plot(dg1, label='dg1')
#    plt.ylim([-1.5,1.5])
#    
#    plt.subplot(3,3,3)
#    plt.plot(dg2, label='dg2')
#    plt.ylim([-1.5,1.5])
#    
#    plt.subplot(3,3,4)
#    plt.plot(dg3, label='dg3')
#    plt.ylim([-1.5,1.5])
#    
#    plt.subplot(3,3,5)
#    plt.plot(dg4, label='dg4')
#    plt.ylim([-1.5,1.5])
#    
#    plt.subplot(3,3,6)
#    plt.plot(dg5, label='dg5')
#    plt.ylim([-1.5,1.5])
#    
#    plt.subplot(3,3,7)
#    plt.plot(dg6, label='dg6')    
#    plt.ylim([-1.5,1.5])
#    
#    plt.subplot(3,3,8)
#    plt.plot(dg7, label='dg7')
#    plt.ylim([-1.5,1.5])
#    
#    plt.subplot(3,3,9)
#    plt.plot(dg8, label='dg8')
#    plt.ylim([-1.5,1.5])
    

    for i in range (0,N):
        color = 'C' + str(c)
        c = c+1
        dg = (data[i,:])
        p = plt.plot(dg, alpha=0.1)
    
        
plt.savefig('../../data/dg_twist_49500_2.jpg') 
plt.show()