# -*- coding: utf-8 -*-
"""
read 3x3 force matrix for each fiber and add them to find for the centerline

@author: zahra
"""
import numpy as np

path = 'D:/sandbox/fiberSimulation/dataSets/train_teeth1220/'
vrtNum = 299
fiberNum = 160

for i in range (0,59):
    if i*5 < 10 :
        frameNum = '0000'+ str(i*5) + '00'
    elif i*5 <100 :
       frameNum = '000'+ str(i*5) + '00' 
    else:
       frameNum = '00'+ str(i*5) + '00' 
    fname_read = path + 'frame_' + frameNum + 'fiber_00.f3'
    fname_write = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/1220/matrix_f3_' + str(i*5) + '.txt'
    
    force = np.zeros(vrtNum*9).reshape(vrtNum,9)
    with open(fname_write, 'w') as fout:
        with open(fname_read, 'r') as fin:
            for f in range (0,fiberNum): 
                for v in range (0,vrtNum):
                    line = np.fromstring( fin.readline(), dtype=float, sep=' ' )
                    force[v] = force[v] + line
                    
            
            for v in range (0,vrtNum):
                for i in range (0,9):
                    fout.writelines('%.8f ' % (force[v,i]) )
                fout.writelines('\n')
                
        fout.close()
        
    print('done!')