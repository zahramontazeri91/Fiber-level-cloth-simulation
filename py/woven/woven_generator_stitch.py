# -*- coding: utf-8 -*-
"""
Created on Thu May 17 17:08:20 2018

@author: zahra
"""
import re
import numpy as np
import matplotlib.image as image

# read pattern
sz = 10 #must be divisible by 5
img = image.imread('100x100_pattern_1.png')
arr = np.array(img)
data = np.array(arr[:,:,0])
pattern_digit = 5
yarnType = 'yarn8/'
init = 70
cylinderLen = 6 #vertices covered in each segment of 5-digit pattern
z_step_size = 0.07
patternLen = cylinderLen*pattern_digit*0.07
# In[]:
def cropPattern(p):
    objFile = 'pattern/' + yarnType + p + '_1.0x.obj'
    feFile = 'pattern/' + yarnType + p + '_1.0x.fe'
    objSeg=[]
    feSeg=[] 
    twistSeg=[]

    p_start = init 
    p_end = p_start + 2 * (cylinderLen*pattern_digit) #mult by 2 to cover circles
    
    fin_obj=open(objFile)
    lines=fin_obj.readlines()
    for i in range (p_start, p_end+1):
        particle = ( float(lines[i].split()[1]), float(lines[i].split()[2]), float(lines[i].split()[3]) ) #first split is 'v'
        objSeg.append(np.array(particle) )
        twist = float(lines[i].split()[4])
        twistSeg.append(twist)

    fin_fe=open(feFile)
    lines=fin_fe.readlines()    
    for j in range (p_start, p_end+1):
        dg = np.array([ [ float(lines[j].split()[0]), float(lines[j].split()[1]), float(lines[j].split()[2]) ], \
                        [ float(lines[j].split()[3]), float(lines[j].split()[4]), float(lines[j].split()[5]) ], \
                        [ float(lines[j].split()[6]), float(lines[j].split()[7]), float(lines[j].split()[8]) ] ])
        
        feSeg.append(dg)
        
    
    return objSeg, feSeg, twistSeg    
# In[]:
def preload():
    p = '00011'
    obj_00011, fe_00011, twist_00011 = cropPattern(p)
    p = '10100'
    obj_10100, fe_10100, twist_10100 = cropPattern(p)
    p = '11110'
    obj_11110, fe_11110, twist_11110 = cropPattern(p)
    return obj_00011, fe_00011, twist_00011, obj_10100, fe_10100, twist_10100, obj_11110, fe_11110, twist_11110


# In[]:        
def updateCurve(obj_ref, fe_ref, twist_ref, idx, isFlip):
    offset = idx*cylinderLen
    obj = np.array(obj_ref[offset:offset + pattern_digit*cylinderLen])
    fe = np.array(fe_ref[offset:offset + pattern_digit*cylinderLen])
    twist = twist_ref[offset:offset + pattern_digit*cylinderLen]
    
    if (isFlip):
        rotateZ = np.array([ [-1,0,0],[0,-1,0],[0,0,1] ])
        obj = np.matmul(obj, rotateZ)
        fe = np.matmul(fe, rotateZ)

    return obj, fe, twist
# In[]:
def shiftCurve(obj, translateZ, translateX):
    obj[:,2] = obj[:,2] + translateZ
    obj[:,0] = obj[:,0] + translateX
    return obj

# In[]:
def rotateCurve(obj, fe):
    rotateY = np.array([ [0,0,-1],[0,1,0],[1,0,0] ])
    obj = np.matmul(obj, rotateY)
    fe = np.matmul(fe, rotateY)
    return obj, fe
    
# In[]:
def whichpattern(query):
    patterns = ["00011", "10100", "11110", "11100", "01011", "00001" ]
    isFlip = 0
    for p in patterns:
    #    print(p)
        found = re.search(query, p+p)
        if (found):
            idx = found.start() # founf the pattern and index
            break
            
    pattern_idx = patterns.index(p)
    half = int (len(patterns)/2 )
    if (pattern_idx >= half ): 
        isFlip = 1
        p = patterns[pattern_idx - half]
    return p, idx, isFlip

# In[]:
def writefiles(obj, fe, twist, fn_obj, fn_fe):
    print(obj.shape[0], fe.shape[0])
    assert (obj.shape[0] == twist.shape[0] )
    assert (obj.shape[0] == fe.shape[0] )
    with open(fn_obj, "w") as fout_obj:
        for i in range (1, obj.shape[0]): #start from line1 because o is trash thanks to numpy initization
            fout_obj.writelines('v %.8f %.8f %.8f %.8f\n' % (obj[i][0], obj[i][1], obj[i][2], twist[i]) )
        for i in range (1, obj.shape[0]-1 ):
            fout_obj.writelines('l %d %d \n' % (i, i+1) )

    with open(fn_fe, "w") as fout_fe:
        for i in range (1, fe.shape[0]-1): #start from line1 because 0 is trash thanks to numpy initization
            fout_fe.writelines('%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f \n' % \
                               (fe[i][0][0],fe[i][0][0], fe[i][0][0], \
                               fe[i][0][0],fe[i][0][0], fe[i][0][0], \
                               fe[i][0][0],fe[i][0][0], fe[i][0][0] ) )
            
# In[]:
#This is for cropping both centerline and deformation from reference to generate arbitrary woven
yarn_num = 6

obj_00011, fe_00011, twist_00011, obj_10100, fe_10100, twist_10100, obj_11110, fe_11110, twist_11110 = preload()


for y in range (0,yarn_num):

    fn_obj = 'yarn/woven_flower_' + str(y) + '.obj'
    fn_fe = 'yarn/woven_flower_' + str(y) + '.fe'
    seg = 0
    translateX = 0.1 * y
    
    fe_all = np.empty([1,3,3])
    obj_all = np.empty([1,3])
    twist_all = np.empty([1])   

    for d in range (0, sz-pattern_digit + 1, pattern_digit):
        print ( ' *** pattern ****' )
        q = ""
        for w in range (0,pattern_digit):
            q += str(int(data[y][d+w]) )
         
        print('query: ', q)
        query = q
        p, idx, isFlip = whichpattern(query)    
        print('patern: ', p, idx, isFlip)
        
        if (p=='00011'):
            obj_ref, fe_ref, twist_ref = obj_00011, fe_00011, twist_00011
        elif (p=='10100'):
            obj_ref, fe_ref, twist_ref = obj_10100, fe_10100, twist_10100
        elif (p=='11110'):     
            obj_ref, fe_ref, twist_ref = obj_11110, fe_11110, twist_11110
        
        obj, fe, twist = updateCurve(obj_ref, fe_ref, twist_ref, idx, isFlip)
    #    print(obj.shape, fe.shape)
        translateZ = 2.2 * seg #######################
        seg = seg+1
        print(translateZ, translateX)
        shiftCurve(obj, translateZ, translateX)
        
        
        obj_all = np.append(obj_all, obj, axis=0)
        fe_all = np.append(fe_all, fe, axis=0)
        twist_all = np.append(twist_all, twist, axis=0)
     
    if (y >= int(yarn_num/2) ) :
        obj_all, fe_all = rotateCurve(obj_all, fe_all)
    writefiles(obj_all, fe_all, twist_all, fn_obj, fn_fe)  
    


        
    
    










