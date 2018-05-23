# -*- coding: utf-8 -*-
"""
Created on Tue May 22 12:39:42 2018

@author: zahra
"""

import numpy as np
import matplotlib.image as image
from sklearn.preprocessing import normalize

# In[]: 
#load pattern
def loadPattern(fn):
    img = image.imread(fn)
    arr = np.array(img)
    pattern = np.array(arr[:,:,0])  
    return pattern
# In[]: 
def writeOBJ(fn_obj, all_x, all_y, all_z):

    with open(fn_obj, "w") as fout_obj:
        for i in range (1, all_x.shape[0]): #start from line1 because o is trash thanks to numpy initization
            fout_obj.writelines('v %.8f %.8f %.8f\n' % (all_x[i], all_y[i], all_z[i]) )
        for i in range (1, all_x.shape[0]-1 ):
            fout_obj.writelines('l %d %d \n' % (i, i+1) )
            
# In[]: 
def writeFE(fn_fe, all_x, all_y, all_z):
    with open(fn_fe, "w") as fout_fe:
        for i in range (0, all_x.shape[0] - 1):
            v0 = np.array( [all_x[i], all_y[i], all_z[i] ] )
            v1 = np.array( [all_x[i+1], all_y[i+1], all_z[i+1] ] )
            
            tan = ( v1-v0 )
#            tan = normalize(tan[:,np.newaxis], axis=0).ravel()
            
            normal = np.array( [1,0,0 ] )
            if (v0[0] == v1[0]): #same x so normal is x axis
                normal = np.array( [1,0,0 ] )
            if (v0[1] == v1[1]): #same z so normal is z axis
                normal = np.array( [0,0,1 ] )
                
            binorm = np.cross(tan, normal)
    
            fout_fe.writelines('%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f \n' % \
                               (tan[0],normal[0], binorm[0], \
                               tan[1],normal[1], binorm[1], \
                               tan[2],normal[2], binorm[2] ) )            
            
# In[]:     
sample = 7
height = 0.4
segLen = height*2 #because of spacing1.0x and so tanh is simpler
top = height/2.0
bottom = -1.0*top
step = segLen/sample
path = "D:/sandbox/fiberSimulation/dataSets/woven/arbitrary_pattern/100x100/yarn/"

sz = 100
pattern_type = '100x100'
fn_pattern = pattern_type + '_pattern_1.png'
pattern = loadPattern(fn_pattern)

noise = np.random.normal(0,0.001,sz)

# write yarns along x 
for yarn in range(0,sz):
    all_x = np.empty([1])
    all_y = np.empty([1])
    all_z = np.empty([1])
    
    
    for seg in range(0,sz-1):
        q = str(int(pattern[yarn,seg]) ) + str(int(pattern[yarn,seg+1]) )
#        print(q)
        
        start = segLen*seg
        end = start + segLen - step
        x = np.arange(start, end, step)
        dis = segLen*yarn
        z = np.full( (int(x.shape[0])), dis+noise[seg])
        all_x = np.append(all_x, x, axis=0)
        all_z = np.append(all_z, z, axis=0)
        compress = 4.0/segLen #4 because tanh lies between -2 to 2 in order to have amplitude between [-1,1]
        x_tanh = x - start - segLen/2.0
        
        
        y = height/2.0 * np.tanh(compress*x_tanh)
         
#        if (q=='00'):
#            y = np.full( (int(x.shape[0])), bottom)
##            tg_y = 0
#        elif (q=='11'):
#            y = np.full( (int(x.shape[0])), top) 
#            tg_y = 0
#        elif (q=='01'):
#            y = height/2.0 * np.tanh(compress*x_tanh)
##            tg_y =  height/2.0 *compress * (1.0 - np.tanh(compress*x_tanh) )
#        elif (q=='10'):         
#            y = -1.0 * height/2.0 * np.tanh(compress*x_tanh)

    #    print(y)
        
#        all_tg_y = np.append(all_tg_y, tg_y, axis=0)
        all_y = np.append(all_y, y+noise[seg], axis=0)
    
    fn_fe = path + '/frame_0000000fiber_' + str(yarn).zfill(2) + '.fe'
    fn_obj = path + '/frame_0000000fiber_' + str(yarn).zfill(2) + '.obj'
    print(fn_obj)
    writeOBJ(fn_obj, all_x, all_y, all_z)
    writeFE(fn_fe, all_x, all_y, all_z)
    
# write yarns along z
for yarn in range(0,sz):
    all_x = np.empty([1])
    all_y = np.empty([1])
    all_z = np.empty([1])
    for seg in range(0,sz-1):
        q = str(int(1 - pattern[seg, yarn]) ) + str(int(1 - pattern[seg+1, yarn]) ) #inverse the pattern
#        print(q)
        
        start = segLen*seg
        end = start + segLen - step
        z = np.arange(start, end, step)
        dis = segLen*yarn
        x = np.full( (int(z.shape[0])), dis+noise[seg])
        all_x = np.append(all_x, x, axis=0)
        all_z = np.append(all_z, z, axis=0)
        compress = 4.0/segLen #4 because tanh lies between -2 to 2 in order to have amplitude between [-1,1]
        z_tanh = z - start - segLen/2.0
        
        y = height/2.0 * np.tanh(compress*z_tanh)
        
#        if (q=='00'):
#            y = np.full( (int(x.shape[0])), bottom)
#        elif (q=='11'):
#            y = np.full( (int(x.shape[0])), top)  
#        elif (q=='01'):
#            y = height/2.0 * np.tanh(compress*z_tanh)
#        elif (q=='10'):         
#            y = -1.0 * height/2.0 * np.tanh(compress*z_tanh)
#        print(y)
        
        all_y = np.append(all_y, y+noise[seg], axis=0)
    
    fn_fe = path + '/frame_0000000fiber_' + str(yarn+sz).zfill(2) + '.fe'
    fn_obj = path + '/frame_0000000fiber_' + str(yarn+sz).zfill(2) + '.obj'
    print(fn_obj)
    writeOBJ(fn_obj, all_x, all_y, all_z)
    writeFE(fn_fe, all_x, all_y, all_z)
    
 # In[]: GENERATE DG   
# deformation gradient in simulation space:
# Multiplying DG to x-axis would give us the tangent
# So A.[1 0 0] = tan
# The first column of A should be the tan and the other two column are two orthonormal vector.
# We can hard code one of them, since all the curves are in 2D and get the other with cross-product



     
    
    
    