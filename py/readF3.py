# -*- coding: utf-8 -*-
"""
read 3x3 force matrix for each fiber and add them to find for the centerline

@author: zahra
"""
import numpy as np

path = 'D:/sandbox/fiberSimulation/dataSets/fiber_data_fortraining_zahra/'
fname_read = path + 'frame_0006000fiber_00.f3'
fname_write = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/NN/force_frame_0006000.txt'

force = np.zeros(149*9).reshape(149,9)

with open(fname_write, 'w') as fout:
    with open(fname_read, 'r') as fin:
        for f in range (0,160): 
            for v in range (0,149):
                line = np.fromstring( fin.readline(), dtype=float, sep=' ' )
                force[v] = force[v] + line
                
        
        for v in range (0,149):
            for i in range (0,9):
                fout.writelines('%.8f ' % (force[v,i]) )
            fout.writelines('\n')
            
    fout.close()
    
print('done!')