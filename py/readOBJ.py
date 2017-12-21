# -*- coding: utf-8 -*-
"""
Read obj file and write in the data format of simulated data
"""
vrtNum = 300
path = 'D:/sandbox/fiberSimulation/dataSets/train_teeth1220/'
for i in range (0,36):
    if i*5 < 10 :
        frameNum = '0000'+ str(i*5) + '00'
    elif i*5 <100 :
       frameNum = '000'+ str(i*5) + '00' 
    else:
       frameNum = '00'+ str(i*5) + '00' 

    fname_read = path + 'frame_' + frameNum + 'fiber_00.obj'
    fname_write = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/data/1220_frame_' + frameNum + 'fiber_00.txt'
    with open(fname_write, 'w') as fout:
        with open(fname_read, 'r') as fin:
            fout.writelines('%d\n' %160)
            for f in range (0,160):
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
                
    

