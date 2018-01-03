# -*- coding: utf-8 -*-
"""
read 3x3 force matrix for each fiber and add them to find for the centerline

@author: zahra
"""
import numpy as np




# In[]:
def readExternal(fname_read):
    with open(fname_read, 'r') as fin:
        force = np.zeros((ds_vrtNum-1)*9).reshape((ds_vrtNum-1),9)
        for v in range (0,ds_vrtNum-1):
            line = np.fromstring( fin.readline(), dtype=float, sep=' ' )
            force[v] = line
        # use average to use them for vertx not segment
        force1 = np.zeros((ds_vrtNum)*9).reshape((ds_vrtNum),9)
        force1[0] = force[0]
        for v in range (0,len(force)-1):
            force1[v+1] = (force[v] + force[v+1] )/2.0
        force1[len(force)] = force[len(force)-1]        
        
        
        # interpolate two values between each two
        force2 = np.zeros((vrtNum)*9).reshape((vrtNum),9)
        force2[0] = force1[0]
        i = 1
        for v in range (0,len(force1)-1):
            force2[i] = force1[v]
            force2[i+1] = (force1[v+1] - force1[v])/3.0 + force1[v]
            force2[i+2] = 2.0*(force1[v+1] - force1[v])/3.0 + force1[v]
            i = i+3
        force2[i] = force1[len(force1)-1]
        force2[i+1] = force1[len(force1)-1]

        return force2
    print('read external done!')   

# In[]:
def readInternal(fname_read):
    with open(fname_read, 'r') as fin:
        force = np.zeros((ds_vrtNum)*6).reshape((ds_vrtNum),6)
        for v in range (0,ds_vrtNum):
            line = np.fromstring( fin.readline(), dtype=float, sep=' ' )
            force[v] = line      
        
        # interpolate two values between each two
        force2 = np.zeros((vrtNum)*6).reshape((vrtNum),6)
        force2[0] = force[0]
        i = 1
        for v in range (0,len(force)-1):
            force2[i] = force[v]
            force2[i+1] = (force[v+1] - force[v])/3.0 + force[v]
            force2[i+2] = 2.0*(force[v+1] - force[v])/3.0 + force[v]
            i = i+3
        force2[i] = force[len(force)-1]
        force2[i+1] = force[len(force)-1]

        return force2
    print('read internal done!')      
    
# main
# In[]: 
#path = 'D:/sandbox/fiberSimulation/dataSets/train_teeth1231_ready/'
path = 'D:/sandbox/fiberSimulation/dataSets/train_stretch1233_ready/'
dataset = '1233'
vrtNum = 300
ds_vrtNum = vrtNum/3
wrt_path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/' + dataset 
for i in range (0,18):
    f = i * 10
    if f < 10 :
        frameNum = '0000'+ str(f) + '00'
    elif f <100 :
       frameNum = '000'+ str(f) + '00' 
    else:
       frameNum = '00'+ str(f) + '00' 
    fn_write = wrt_path + '/physical_' + str(f) + '.txt'
    fn_read_ext = path + 'frame_' + frameNum + 'fiber_00_RED.fe'
    fn_read_int = path + 'frame_' + frameNum + 'fiber_00.sforce'
    
    ext_force = np.zeros((vrtNum)*9).reshape((vrtNum),9)
    ext_force =  readExternal(fn_read_ext) 
    
    int_force = np.zeros((vrtNum)*6).reshape((vrtNum),6)
    int_force =  readInternal(fn_read_int)     
    
    with open(fn_write, 'w') as fout:
            for v in range (0,vrtNum):
                for i in range (0,9): # write 3x3 force matrix
                    fout.writelines('%.8f ' % (ext_force[v,i]) )
                for i in range (0,6): # write internal force
                    fout.writelines('%.8f ' % (int_force[v,i]) )
                fout.writelines('\n')