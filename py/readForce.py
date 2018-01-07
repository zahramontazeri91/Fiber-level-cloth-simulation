# -*- coding: utf-8 -*-
"""
read 3x3 force matrix for each fiber and add them to find for the centerline

@author: zahra
"""
import numpy as np
from sklearn.preprocessing import MinMaxScaler

    
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

# In[]:
def readCenterYarn(fname_read):
    with open(fname_read, 'r') as fin:
        cntr = np.zeros((ds_vrtNum)*3).reshape((ds_vrtNum),3)
        twist = np.zeros((ds_vrtNum))
        for v in range (0,ds_vrtNum):
            line = fin.readline().split()
            cntr[v,0] = float(line[1]) * 0.25
            cntr[v,1] = float(line[2]) * 0.25
            cntr[v,2] = float(line[3]) * 0.25
#            twist[v] = float(line[4]) 
        
        # interpolate two values between each two
        twist2 = np.zeros((vrtNum))
        twist2[0] = twist[0]
        i = 1
        for v in range (0,len(twist)-1):
            twist2[i] = twist[v]
            twist2[i+1] = (twist[v+1] - twist[v])/3.0 + twist[v]
            twist2[i+2] = 2.0*(twist[v+1] - twist[v])/3.0 + twist[v]
            i = i+3
        twist2[i] = twist[len(twist)-1]
        twist2[i+1] = twist[len(twist)-1]

        return cntr, twist2
    print('read centerline done!')     
# main
# In[]: 
path = 'D:/sandbox/fiberSimulation/dataSets/teeth_spacing2_ready/'
#path = 'D:/sandbox/fiberSimulation/dataSets/training_tmp/'
#path = 'D:/sandbox/fiberSimulation/dataSets/train_teeth1231_ready/'
#path = 'D:/sandbox/fiberSimulation/dataSets/train_stretch1233_ready/'
dataset = '1231_2' #####CHANGE FILE FORMAT CENTERLINE FOR TEST DATA#####
vrtNum = 300
ds_vrtNum = vrtNum/3
skipFactor = 5 
wrt_path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/' + dataset 
for i in range (0,180/5 + 1):
    f = i * skipFactor
    if f < 10 :
        frameNum = '0000'+ str(f) + '00'
    elif f <100 :
        frameNum = '000'+ str(f) + '00' 
    else:
       frameNum = '00'+ str(f) + '00' 
    fn_write_force = wrt_path + '/physicalParam/physical_' + str(f) + '_world.txt'
    fn_read_ext = path + 'frame_' + frameNum + 'fiber_00.fe'
    fn_read_int = path + 'frame_' + frameNum + 'fiber_00.sforce'
    fn_read_center = path + 'frame_' + frameNum + 'fiber_00_RED.obj'
    
    ext_force = np.zeros((vrtNum)*9).reshape((vrtNum),9)
    ext_force =  readExternal(fn_read_ext) 
    
    int_force = np.zeros((vrtNum)*6).reshape((vrtNum),6)
    int_force =  readInternal(fn_read_int)     
    scaler = MinMaxScaler(feature_range=(0, 1))
    scaler.fit(int_force)
    int_force = scaler.transform(int_force)
    
    with open(fn_write_force, 'w') as fout:
            for v in range (0,vrtNum):
                for i in range (0,9): # write 3x3 force matrix
                    fout.writelines('%.8f ' % (ext_force[v,i]) )
                for i in range (0,6): # write internal force
                    fout.writelines('%.8f ' % (int_force[v,i]) )
                fout.writelines('\n')

    centerYarn = np.zeros((vrtNum)*3).reshape((vrtNum),3)
    twist = np.zeros(vrtNum)          
    centerYarn, twist = readCenterYarn(fn_read_center);

    
    fn_write_cntr = wrt_path + '/centerYarn_' + str(f) + '_ds.txt' 
    with open(fn_write_cntr, 'w') as fout:
        fout.writelines('%d \n' % (ds_vrtNum) )
        for v in range (0,ds_vrtNum):
                fout.writelines('%.8f %.8f %.8f \n' % (centerYarn[v,0], centerYarn[v,1], centerYarn[v,2]) )
            
#    fn_write_twist = wrt_path + '/twist_' + str(f) + '.txt' 
#    with open(fn_write_twist, 'w') as fout:
#        r = 0.0
#        fout.writelines('%d \n' % (vrtNum) )
#        for v in range (0,vrtNum):         
#            r = r + twist[v]
#            fout.writelines('%.8f \n' % (r) )