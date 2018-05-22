# -*- coding: utf-8 -*-
"""
Created on Tue May 22 12:39:42 2018

@author: zahra
"""

import numpy as np
import matplotlib.image as image

# In[]: 
#load pattern
def loadPattern(fn):
    img = image.imread(fn)
    arr = np.array(img)
    pattern = np.array(arr[:,:,0])  
    return pattern
# In[]: 
def writefiles(fn_obj, all_x, all_y, all_z):

    with open(fn_obj, "w") as fout_obj:
        for i in range (1, all_x.shape[0]): #start from line1 because o is trash thanks to numpy initization
            fout_obj.writelines('v %.8f %.8f %.8f\n' % (all_x[i], all_y[i], all_z[i]) )
        for i in range (1, all_x.shape[0]-1 ):
            fout_obj.writelines('l %d %d \n' % (i, i+1) )
            
# In[]:     
cur_x = []
cur_y = []
cur_z = []
sample = 7
height = 0.4
segLen = height*2 #because of spacing1.0x and so tanh is simpler
top = height/2.0
bottom = -1.0*top
step = segLen/sample

sz = 100
pattern_type = '100x100'
fn_pattern = pattern_type + '_pattern_1.png'
pattern = loadPattern(fn_pattern)

# write yarns along x 
for yarn in range(0,sz):
    all_x = np.empty([1])
    all_y = np.empty([1])
    all_z = np.empty([1])
    for seg in range(0,sz-1):
        q = str(int(pattern[yarn,seg]) ) + str(int(pattern[yarn,seg+1]) )
#        print(q)
        
        start = segLen*seg
        end = start + segLen
        x = np.arange(start, end, step)
        dis = segLen*yarn
        z = np.full( (int(x.shape[0])), dis)
        all_x = np.append(all_x, x, axis=0)
        all_z = np.append(all_z, z, axis=0)
        compress = 4.0/segLen #4 because tanh lies between -2 to 2 in order to have amplitude between [-1,1]
        
        if (q=='00'):
            y = np.full( (int(x.shape[0])), bottom)
        elif (q=='11'):
            y = np.full( (int(x.shape[0])), top)  
        elif (q=='01'):
            x_tanh = x - start - segLen/2.0
            y = height/2.0 * np.tanh(compress*x_tanh)
        elif (q=='10'):
            x_tanh = x - start - segLen/2.0
            y = -1.0 * height/2.0 * np.tanh(compress*x_tanh)
    #    print(y)
        
        all_y = np.append(all_y, y, axis=0)
    
    
    fn_obj = pattern_type + '/centerline_' + str(yarn).zfill(2) + '.obj'
    print(fn_obj)
    writefiles(fn_obj, all_x, all_y, all_z)
    
# write yarns along z
for yarn in range(0,sz):
    all_x = np.empty([1])
    all_y = np.empty([1])
    all_z = np.empty([1])
    for seg in range(0,sz-1):
        q = str(int(1 - pattern[seg, yarn]) ) + str(int(1 - pattern[seg+1, yarn]) ) #inverse the pattern
#        print(q)
        
        start = segLen*seg
        end = start + segLen
        z = np.arange(start, end, step)
        dis = segLen*yarn
        x = np.full( (int(z.shape[0])), dis)
        all_x = np.append(all_x, x, axis=0)
        all_z = np.append(all_z, z, axis=0)
        compress = 4.0/segLen #4 because tanh lies between -2 to 2 in order to have amplitude between [-1,1]
        z_tanh = z - start - segLen/2.0
        
        if (q=='00'):
            y = np.full( (int(x.shape[0])), bottom)
        elif (q=='11'):
            y = np.full( (int(x.shape[0])), top)  
        elif (q=='01'):
            y = height/2.0 * np.tanh(compress*z_tanh)
        elif (q=='10'):         
            y = -1.0 * height/2.0 * np.tanh(compress*z_tanh)
#        print(y)
        
        all_y = np.append(all_y, y, axis=0)
    
    
    fn_obj = pattern_type + '/centerline_' + str(yarn+sz).zfill(2) + '.obj'
    print(fn_obj)
    writefiles(fn_obj, all_x, all_y, all_z)
    
    
    
    