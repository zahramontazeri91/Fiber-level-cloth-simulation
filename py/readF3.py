# -*- coding: utf-8 -*-
"""
read 3x3 force matrix for each fiber and add them to find for the centerline

@author: zahra
"""
import numpy as np

#path = 'D:/sandbox/fiberSimulation/dataSets/train_teeth1220/'
#path = 'D:/sandbox/fiberSimulation/dataSets/train_stretching1222_old/'
path = 'D:/sandbox/fiberSimulation/dataSets/train_teeth1220_old/'
dataset = '1220'

vrtNum = 300-1
fiberNum = 160

for i in range (0,30):
    f = i * 10
    if f < 10 :
        frameNum = '0000'+ str(f) + '00'
    elif f <100 :
       frameNum = '000'+ str(f) + '00' 
    else:
       frameNum = '00'+ str(f) + '00' 
    fname_read = path + 'frame_' + frameNum + 'fiber_00.f3'
    fname_write = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/' + dataset + '/physicalParam/physical_' + str(f) + '_world.txt'
    
    force = np.zeros(vrtNum*9).reshape(vrtNum,9)
    with open(fname_write, 'w') as fout:
        with open(fname_read, 'r') as fin:
            for f in range (0,fiberNum): 
                for v in range (0,vrtNum):
                    line = np.fromstring( fin.readline(), dtype=float, sep=' ' )
                    force[v] = force[v] + line
                    
            
            for v in range (0,vrtNum):
                for i in range (0,9): # write 3x3 force matrix
                    fout.writelines('%.8f ' % (force[v,i]) )
                fout.writelines('\n')
                
        fout.close()
        
    print('done!')