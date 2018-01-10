# -*- coding: utf-8 -*-
"""
Take the middle part of the data and append them together

@author: zahra
"""
import numpy as np
from shutil import copyfile
from sklearn.preprocessing import MinMaxScaler

# In[]:
def buildTrainY(vrtNum, trimPercent, first_frame, last_frame, not_frame, stride, filename):
    c=first_frame
    for i in range (first_frame,last_frame+1):
        if i==not_frame:
            continue
        j = c*skipFactor
        fname_write = path + filename + str(j - first_frame*skipFactor) + '.txt'
        with open(fname_write, 'w') as fout:
            fname_read = path + 'matrix_S_' + str(j) + '.txt'
            with open(fname_read, 'r') as fin:
                cs = 0
                while (cs<vrtNum):
                    line = fin.readline()
                    ignor = int(vrtNum*trimPercent)                   
                    if (cs >= ignor and cs <= vrtNum - ignor):
                        fout.writelines(line)
                    cs = cs + stride
        fout.close()                
        c = c+1
    print( 'all TrainY added\n')
# In[]:
def normalizeInternal(vrtNum, first_frame, last_frame):
    c = first_frame
    for i in range (first_frame,last_frame+1):
        j = c*skipFactor
        fname_read = path + 'physical_' + str(j) + '.txt'
        data = np.loadtxt(fname_read)
        internal = data[:,9:]
        scaler = MinMaxScaler(feature_range=(0, 1))
        scaler.fit(internal)
        internal = scaler.transform(internal)
        #append it back to external
        data = np.concatenate((data[:,0:9],internal), axis=1)
        fname_write = path + 'physical_' + str(j) + '_norm.txt'
        np.savetxt(fname_write, data, fmt='%.8f' )
            
        c = c+1
# In[]:
def buildTrainX_conv(vrtNum, trimPercent, first_frame, last_frame, not_frame, sigma, stride, filename):
    c = first_frame
    for i in range (first_frame,last_frame+1):
        if i==not_frame:
            continue
        j = c*skipFactor
        fname_write = path + filename + str(j - first_frame*skipFactor) + '.txt'
        with open(fname_write, 'w') as fout:
            fname_read = path + 'physical_' + str(j) + '_norm.txt'            
            with open(fname_read, 'r') as fin:
                line = fin.read().splitlines() #read the whole file as a list of string
                if (len(line) != vrtNum) :
                    print (len(line), vrtNum)
                assert(len(line) == vrtNum)
                ignor = int(vrtNum*trimPercent)
                cs = 0
                while (cs<vrtNum):                    
                    if (cs >= ignor and cs <= vrtNum - ignor):
                        total = ""
                        for w in range (sigma ,0, -1): 
#                            TO TAKE ONLY PART OF THE PHYSICAL PARAMS
                            tmp = line[cs - w].split()
                            line[cs - w] = ' '.join(tmp[0:12])
                            total += line[cs - w] + ' ' 
                        for w in range (0,sigma + 1 ):  
                            tmp = line[cs + w].split()
                            line[cs + w] = ' '.join(tmp[0:12])                             
                            total += line[cs + w] + ' '
                            
                        fout.writelines(total+'\n')
                    cs = cs + stride
        fout.close()                
        c = c+1
    print( 'all TrainX added\n')
    
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
def appendDatasets(dataset1, dataset2):
    path1 = w_path + dataset1 + '/NN/'
    path2 = w_path + dataset2 + '/NN/'
    X_train_all_1 = np.loadtxt(path1 + "trainX_all.txt",delimiter=None)
    Y_train_all_1 = np.loadtxt(path1 + "trainY_all.txt",delimiter=None)
    X_train_all_2 = np.loadtxt(path2 + "trainX_all.txt",delimiter=None)
    Y_train_all_2 = np.loadtxt(path2 + "trainY_all.txt",delimiter=None)
    
    X_train_all = np.concatenate((X_train_all_1,X_train_all_2), axis=0)
    Y_train_all = np.concatenate((Y_train_all_1,Y_train_all_2), axis=0)
    
    
    np.savetxt(w_path + 'all/NN/' + 'trainX_all.txt', X_train_all, fmt='%.6f', delimiter=' ')
    np.savetxt(w_path + 'all/NN/' + 'trainY_all.txt', Y_train_all, fmt='%.6f', delimiter=' ')
# In[]:
dataset = 'spacing1.5x'
#dataset = 'spacing3.0x_rotate_test'

w_path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'
path = w_path + dataset + '/'
vrtNum = 300
skipFactor = 5
trimPercent = 0.15
first_frame = 0
last_frame = 240/skipFactor 
test_frame = -1
not_frame = -1
test_out = 0 #binary flag  
sigma = 5 #window_size = 2*sigma + 1
stride = 1

normalizeInternal(vrtNum, first_frame, last_frame)

# In[]:
#build train data:

#filename = 'NN/trainY_'
#buildTrainY(vrtNum, trimPercent, first_frame, last_frame, not_frame, stride, filename)
#filename = 'NN/trainX_'
#buildTrainX_conv(vrtNum, trimPercent, first_frame, last_frame, not_frame, sigma, stride, filename)
#appendAll(first_frame, last_frame, 'NN/trainX_', test_out)
#appendAll(first_frame, last_frame, 'NN/trainY_', test_out)
#
#appendDatasets('all', dataset)
# In[]:
#build test data:
filename = 'NN/testX_'
buildTrainX_conv(vrtNum, trimPercent, first_frame, last_frame, not_frame, sigma, stride, filename)


