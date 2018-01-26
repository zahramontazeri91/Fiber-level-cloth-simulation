# -*- coding: utf-8 -*-
"""
Read obj file and write in the data format of simulated data
"""


def readOBJ(vrtNum, fiberNum, totalYarn, lastFrame, dataset, path, skipFactor ):
    for i in range (0,lastFrame/skipFactor + 1):
        f = i*skipFactor
        if f < 10 :
            frameNum = '0000'+ str(f) + '00'
        elif f <100 :
           frameNum = '000'+ str(f) + '00' 
        else:
           frameNum = '00'+ str(f) + '00' 
        for y in range (0,totalYarn):
            if y < 10:
                yarnNum = 'fiber_0' + str(y)
            else:
                yarnNum = 'fiber_' + str(y)
            fname_read = path + dataset+'/yarn/' + 'frame_' + frameNum + yarnNum + '.obj'
            fname_write = path + dataset+'/curves/curve_' + str(f) + '_' + str(y) + '.txt'
            with open(fname_write, 'w') as fout:
                with open(fname_read, 'r') as fin:
#                    fout.writelines('%d\n' %fiberNum)
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
                
    
# In[]:
lastFrame = 980
vrtNum = 121
fiberNum = 1
dataset = 'spacing0.5x_00011'
path = 'D:/sandbox/fiberSimulation/dataSets/fixed_teeth/test/'
skipFactor = 5 
totalYarn = 1
readOBJ(vrtNum, fiberNum, totalYarn, lastFrame, dataset, path, skipFactor )