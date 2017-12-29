# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 09:40:07 2017

@author: zahra
"""

"""
Read obj file and write in the data format of simulated data
"""
vrtNum = 50
fiberNum = 158
yarnNum = 8
#path = 'D:/sandbox/fiberSimulation/dataSets/train_stretching1222/'
#path = 'D:/sandbox/fiberSimulation/dataSets/train_teeth1220/'
path = 'D:/sandbox/fiberSimulation/dataSets/nn_train_woven1224/'
dataset = '1224'

cnt = 0
for i in range (0,77):
    if i*5 < 10 :
        frameNum = '0000'+ str(i*5) + '00'
    elif i*5 <100 :
       frameNum = '000'+ str(i*5) + '00' 
    else:
       frameNum = '00'+ str(i*5) + '00' 
       
    for y in range (0,yarnNum):
        fname_read = path + 'frame_' + frameNum + 'fiber_0' + str(y) + '.obj'
        fname_write = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/data/'+dataset+'/simul_frame_' + str(cnt*5) + '.txt'
        cnt = cnt + 1;
        print(fname_read)
        with open(fname_write, 'w') as fout:
            with open(fname_read, 'r') as fin:
                fout.writelines('%d\n' %fiberNum)
                for f in range (0,fiberNum):
                    fout.writelines('%d\n' %vrtNum) 
                    for v in range (0,vrtNum):
                        line = fin.readline().split()
                        vx = float(line[1]) * 0.25 
                        vy = float(line[2]) * 0.25 
                        vz = float(line[3]) * 0.25 
                        fout.writelines('%.6f %.6f %.6f\n' % (vx, vy, vz) )
        
            fin.close()
        fout.close()    
        print('done!')
                
    