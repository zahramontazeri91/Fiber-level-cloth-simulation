# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 09:40:07 2017

@author: zahra
"""

"""
Read obj file and write in the data format of simulated data
"""
vrtNum = 150
fiberNum = 1
yarnNum = 1
dataset = 'pattern/yarn100/spacing0.5x/10'
f0 = 0
f1 = 18500
skipFactor = 500

pathObj = '../../dataSets/'
pathMtsb = '../../output/'
cnt = 0
for f in range (f0, f1+1, skipFactor):
    frame = str(f).zfill(7) 
    for y in range (0,yarnNum):
        yarn = str(y).zfill(2)
        
        fname_read = pathObj + dataset + '/yarn/frame_' + frame + 'fiber_' + str(yarn) + '.obj'
        fname_write = pathMtsb + dataset + '/curve_' + frame + '_' + str(yarn) + '.txt'

        cnt = cnt + 1;
        print(fname_read)
        with open(fname_write, 'w') as fout:
            with open(fname_read, 'r') as fin:
                fout.writelines('%d\n' %fiberNum)
                for f in range (0,fiberNum):
                    fout.writelines('%d\n' %vrtNum) 
                    for v in range (0,vrtNum):                     
                        line = fin.readline().split()
#                        print(line)
                        vx = float(line[1]) * 0.25 
                        vy = float(line[2]) * 0.25 
                        vz = float(line[3]) * 0.25 
                        fout.writelines('%.6f %.6f %.6f\n' % (vx, vy, vz) )
#        
            fin.close()
        fout.close()    
        print('done!')
                
    