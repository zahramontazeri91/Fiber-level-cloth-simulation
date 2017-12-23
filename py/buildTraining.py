# -*- coding: utf-8 -*-
"""
Take the middle part of the data and append them together

@author: zahra
"""
import numpy as np
from shutil import copyfile


##################################
def buildTrainY(dataset, nb_seg, trimPercent, first_frame, last_frame, not_frame, filename):
    c=first_frame
    for i in range (first_frame,last_frame+1):
        if i==not_frame:
            continue
        j = c*5
        fname_write = path + filename + str(j - first_frame*5) + '.txt'
        with open(fname_write, 'w') as fout:
            if j < 10 :
                fname_read = path + 'param_' + dataset + '_0000' + str(j) + '00.txt'
            elif j < 100 :
                fname_read = path + 'param_' + dataset + '_000' + str(j) + '00.txt'
            else:
                fname_read = path + 'param_' + dataset + '_00' + str(j) + '00.txt'
            with open(fname_read, 'r') as fin:
                for cs in range (0,nb_seg):
                    line = fin.readline()
                    ignor = int(nb_seg*trimPercent)                   
                    if (cs >= ignor and cs <= nb_seg - ignor):
                        fout.writelines(line)
        fout.close()                
        print(str(c) + '-th TrainY added\n')
        c = c+1
#################

def buildTrainX_conv(dataset, nb_seg, trimPercent, first_frame, last_frame, not_frame, sigma, filename):
    c = first_frame
    for i in range (first_frame,last_frame+1):
        if i==not_frame:
            continue
        j = c*5
        fname_write = path + filename + str(j - first_frame*5) + '.txt'
        with open(fname_write, 'w') as fout:
            if j < 10 :
                fname_read = path + 'force_' + dataset + '_0000' + str(j) + '00.txt'
            elif j < 100 :
                fname_read = path + 'force_' + dataset + '_000' + str(j) + '00.txt'
            else:
                fname_read = path + 'force_' + dataset + '_00' + str(j) + '00.txt'            
            with open(fname_read, 'r') as fin:
                line = fin.read().splitlines()
                assert(len(line) == nb_seg)
                ignor = int(nb_seg*trimPercent)  
                for cs in range (0,nb_seg):                 
                    if (cs >= ignor and cs <= nb_seg - ignor):
                                               
                        total = ""
                        for w in range (sigma ,0, -1): 
                            total += line[cs - w] + ' ' 
                        for w in range (0,sigma + 1 ):                               
                            total += line[cs + w] + ' '
                            
                        fout.writelines(total+'\n')
        fout.close()                
        print( str(c) + '-th TrainX added\n')
        c = c+1
#################
def appendAll(first_frame, last_frame, filename, test_out):
    fname_write = path + filename + 'all.txt'
    with open(fname_write, 'w') as fout:
        for i in range (first_frame,last_frame+1-test_out): #-1 because one frame is gone for testing
            fname_read = path + filename + str(i*5 - first_frame*5) + '.txt'
            with open(fname_read, 'r') as fin:
                for j in range (0,nb_seg):
                    line = fin.readline()
                    fout.writelines(line)
    fout.close()                
    print('all files appended!\n')
################
path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/NN/'
dataset = '1220'
nb_seg = 299
trimPercent = 0.2
first_frame = 3
last_frame = 35
test_frame = 32
not_frame = 32
test_out = 1
sigma = 6 #window_size = 2*sigma + 1
filename = 'trainY_'
buildTrainY(dataset, nb_seg, trimPercent, first_frame, last_frame, not_frame, filename)
filename = 'trainX_'
buildTrainX_conv(dataset, nb_seg, trimPercent, first_frame, last_frame, not_frame, sigma, filename)
appendAll(first_frame, last_frame, 'trainX_', test_out)
appendAll(first_frame, last_frame, 'trainY_', test_out)
#build test data:
first_frame = test_frame
last_frame = test_frame
not_frame = 0
filename = 'testY_'
buildTrainY(dataset, nb_seg, trimPercent, first_frame, last_frame, not_frame, filename)
filename = 'testX_'
buildTrainX_conv(dataset, nb_seg, trimPercent, first_frame, last_frame, not_frame, sigma, filename)