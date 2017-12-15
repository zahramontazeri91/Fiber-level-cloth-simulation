# -*- coding: utf-8 -*-
"""
Read obj file and write in the data format of simulated data
"""
path = 'D:/sandbox/fiberSimulation/dataSets/1215/'
fname_read = path + 'frame_0000000fiber_00.obj'
fname_write = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/data/1215_frame0000000.txt'
with open(fname_write, 'w') as fout:
    with open(fname_read, 'r') as fin:
        fout.writelines('%d\n' %160)
        for f in range (0,160):
            fout.writelines('%d\n' %150) 
            for v in range (0,150):
                line = fin.readline().split()
                vx = float(line[1])
                vy = float(line[2])
                vz = float(line[3])
                fout.writelines('%.6f %.6f %.6f\n' % (vx, vy, vz) )
                
    fout.close()
    
print('done!')
                
    

