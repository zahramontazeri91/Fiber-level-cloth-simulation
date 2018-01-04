# -*- coding: utf-8 -*-
"""
Take the middle part of the data and append them together

@author: zahra
"""
import numpy as np
from shutil import copyfile


# In[]:
def buildTrainY(vrtNum, trimPercent, first_frame, last_frame, not_frame, filename):
    c=first_frame
    for i in range (first_frame,last_frame+1):
        if i==not_frame:
            continue
        j = c*skipFactor
        fname_write = path + filename + str(j - first_frame*skipFactor) + '.txt'
        with open(fname_write, 'w') as fout:
            fname_read = path + 'matrix_S_' + str(j) + '.txt'
            with open(fname_read, 'r') as fin:
                for cs in range (0,vrtNum):
                    line = fin.readline()
                    ignor = int(vrtNum*trimPercent)                   
                    if (cs >= ignor and cs <= vrtNum - ignor):
                        fout.writelines(line)
        fout.close()                
        print(str(c) + '-th TrainY added\n')
        c = c+1
# In[]:

def buildTrainX_conv(vrtNum, trimPercent, first_frame, last_frame, not_frame, sigma, filename):
    c = first_frame
    for i in range (first_frame,last_frame+1):
        if i==not_frame:
            continue
        j = c*skipFactor
        fname_write = path + filename + str(j - first_frame*skipFactor) + '.txt'
        with open(fname_write, 'w') as fout:
            fname_read = path + 'physical_' + str(j) + '.txt'            
            with open(fname_read, 'r') as fin:
                line = fin.read().splitlines() #read the whole file as a list of string
                assert(len(line) == vrtNum)
                ignor = int(vrtNum*trimPercent)  
                for cs in range (0,vrtNum):                 
                    if (cs >= ignor and cs <= vrtNum - ignor):
                        total = ""
                        for w in range (sigma ,0, -1): 
                            total += line[cs - w] + ' ' 
                        for w in range (0,sigma + 1 ):                               
                            total += line[cs + w] + ' '
                            
                        fout.writelines(total+'\n')
        fout.close()                
        print( str(c) + '-th TrainX added\n')
        c = c+1

# In[]:
def appendAll(first_frame, last_frame, filename, test_out):
    fname_write = path + filename + 'all.txt'
    with open(fname_write, 'w') as fout:
        for i in range (first_frame,last_frame+1-test_out): #-1 because one frame is gone for testing
            fname_read = path + filename + str(i*skipFactor - first_frame*skipFactor) + '.txt'
            with open(fname_read, 'r') as fin:
                for j in range (0,vrtNum):
                    line = fin.readline()
                    fout.writelines(line)
    fout.close()                
    print('all files appended!\n')

# In[]:
def appendDatasets():
    path1 = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/1231/NN/'
    path2 = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/1233/NN/'
    X_train_all_1 = np.loadtxt(path1 + "trainX_all.txt",delimiter=None)
    Y_train_all_1 = np.loadtxt(path1 + "trainY_all.txt",delimiter=None)
    X_train_all_2 = np.loadtxt(path2 + "trainX_all.txt",delimiter=None)
    Y_train_all_2 = np.loadtxt(path2 + "trainY_all.txt",delimiter=None)
    
    X_train_all = np.concatenate((X_train_all_1,X_train_all_2), axis=0)
    Y_train_all = np.concatenate((Y_train_all_1,Y_train_all_2), axis=0)
    
    w_path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'
    np.savetxt(w_path+'trainX_all_1231_1233.txt', X_train_all, fmt='%.6f', delimiter=' ')
    np.savetxt(w_path+'trainY_all_1231_1233.txt', Y_train_all, fmt='%.6f', delimiter=' ')
# In[]:
dataset = '1231'
path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/'
vrtNum = 300
skipFactor = 10
trimPercent = 0.2
first_frame = 0
last_frame = 17
test_frame = 3
not_frame = -1
test_out = 0
sigma = 6 #window_size = 2*sigma + 1
filename = 'NN/trainY_'
buildTrainY(vrtNum, trimPercent, first_frame, last_frame, not_frame, filename)
filename = 'NN/trainX_'
buildTrainX_conv(vrtNum, trimPercent, first_frame, last_frame, not_frame, sigma, filename)
appendAll(first_frame, last_frame, 'NN/trainX_', test_out)
appendAll(first_frame, last_frame, 'NN/trainY_', test_out)

#appendDatasets()

# In[]:
#build test data:
#first_frame = 1
#last_frame = 35
#not_frame = 0
#filename = 'NN/testY_'
#buildTrainY(nb_seg, trimPercent, first_frame, last_frame, not_frame, filename)
#filename = 'NN/testX_'
#buildTrainX_conv(nb_seg, trimPercent, first_frame, last_frame, not_frame, sigma, filename)


