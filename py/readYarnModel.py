# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 18:04:21 2017
Read yarn model output x y z twist and write to matrix_R and matrix_S format
@author: zahra
"""

import numpy as np

#path = 'D:/sandbox/fiberSimulation/dataSets/train_teeth1220/'
path = 'D:/sandbox/fiberSimulation/dataSets/train_stretching1222/'
#path = 'D:/sandbox/fiberSimulation/dataSets/nn_train_woven1224/'
dataset = '1222'

vrtNum = 300

for i in range (0,30):
    f = i * 10
    if f < 10 :
        frameNum = '0000'+ str(f) + '00'
    elif f <100 :
       frameNum = '000'+ str(f) + '00' 
    else:
       frameNum = '00'+ str(f) + '00' 
    fname_read = path + 'frame_' + frameNum + 'fiber_00.f3'
    fname_write_twist = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/' + dataset + '/matrix_R_' + str(f) + '.txt'
    fname_write_curve = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/' + dataset + '/centerYarn_' + str(f) + '.txt'
 
    
    with open(fname_write_twist, 'w') as fout_R:
        with open(fname_write_curve, 'w') as fout_a:
            with open(fname_read, 'r') as fin:
                fout_a.writelines('%d\n ' % (vrtNum) )
                fout_R.writelines('%d\n ' % (vrtNum) )
                R = 0.0
                for f in range (0,vrtNum):                         
                    line = np.fromstring( fin.readline(), dtype=float, sep=' ' )
                    fout_a.writelines('%.6f %.6f %.6\n ' % (line[0], line[1], line[2]) )
                    R = R + line[3]
                    fout_R.writelines('%.6f\n ' % (R) )

            fout_a.close()
        fout_R.close()
    print('done!')

