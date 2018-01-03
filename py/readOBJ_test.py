# -*- coding: utf-8 -*-
"""
Created on Wed Jan 03 08:50:28 2018

@author: zahra
"""

# -*- coding: utf-8 -*-
"""
Read obj file and write in the data format of simulated data
"""
vrtNum = 300/3
path = 'D:/sandbox/fiberSimulation/dataSets/test_teeth1231_ready/'
dataset = '1231_test'

for i in range (0,18):
    f = i*10
    if f < 10 :
        frameNum = '0000'+ str(f) + '00'
    elif f <100 :
       frameNum = '000'+ str(f) + '00' 
    else:
       frameNum = '00'+ str(f) + '00' 

    fname_read = path + 'frame_' + frameNum + 'fiber_00.obj'
    print(f)
    fname_write = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/centerYarn_' + str(f) + '.txt'
    with open(fname_write, 'w') as fout:
        with open(fname_read, 'r') as fin:
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
                
    

