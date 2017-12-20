# -*- coding: utf-8 -*-
"""
Take the middle part of the data and append them together

@author: zahra
"""
import numpy as np
from shutil import copyfile

path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/NN/'


def buildTrainY(nb_seg, trimPercent, first_frame, last_frame):
    for i in range (first_frame,last_frame+1):
        fname_write = path + 'trainY_' + str(i-first_frame) + '.txt'
        with open(fname_write, 'w') as fout:
            if i < 10 :
                fname_read = path + 'param_frame_000' + str(i) + '000.txt'
            else:
                fname_read = path + 'param_frame_00' + str(i) + '000.txt'
            with open(fname_read, 'r') as fin:
                for cs in range (0,nb_seg):
                    line = fin.readline()
                    ignor = int(nb_seg*trimPercent)                   
                    if (cs >= ignor and cs <= nb_seg - ignor):
                        fout.writelines(line)
        fout.close()                
        print('TrainY added\n')
#################

def buildTrainX_conv(nb_seg, trimPercent, first_frame, last_frame, sigma):
    for i in range (first_frame,last_frame+1):
        fname_write = path + 'trainX_' + str(i-first_frame) + '.txt'
        with open(fname_write, 'w') as fout:
            if i < 10 :
                fname_read = path + 'force_frame_000' + str(i) + '000.txt'
            else:
                fname_read = path + 'force_frame_00' + str(i) + '000.txt'
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
        print('TrainX added\n')

def appendAll(first_frame, last_frame, filename):
    fname_write = path + filename + 'all.txt'
    with open(fname_write, 'w') as fout:
        for i in range (first_frame,last_frame+1):
            if i < 10 :
                fname_read = path + filename + str(i-first_frame) + '.txt'
            else:
                fname_read = path + filename + str(i-first_frame) + '.txt'
            with open(fname_read, 'r') as fin:
                for j in range (0,nb_seg):
                    line = fin.readline()
                    fout.writelines(line)
    fout.close()                
    print('all files appended!\n')
################
def buildTestX_conv(nb_seg, trimPercent, first_frame, last_frame, sigma):
    for i in range (first_frame,last_frame+1):
        fname_write = path + 'testX_' + str(i-first_frame) + '.txt'
        with open(fname_write, 'w') as fout:
            if i < 10 :
                fname_read = path + 'force_frame_000' + str(i) + '000.txt'
            else:
                fname_read = path + 'force_frame_00' + str(i) + '000.txt'
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
        print('TrainX added\n')
#################
nb_seg = 149
trimPercent = 0.40
first_frame = 6
last_frame = 14
sigma = 1
buildTrainY(nb_seg, trimPercent, first_frame, last_frame)
buildTrainX_conv(nb_seg, trimPercent, first_frame, last_frame, sigma)
appendAll(first_frame, last_frame, 'trainX_')
appendAll(first_frame, last_frame, 'trainY_')
#build test data:
first_frame = 15
last_frame = 15
buildTestX_conv(nb_seg, trimPercent, first_frame, last_frame, sigma)