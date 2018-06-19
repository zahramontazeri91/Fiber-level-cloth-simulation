''' open this file in python3.5 where tensorflow is installed
change the parameters for each dataset accordingly
Make sure yarnTypes/dataset matches with the dataset here 
NOTE THAT YARN8 IS HARD CODED'''


import sys
sys.path.insert(0, 'deforGrad')
from deform import main
from NN import main_NN
import os


# set paramters for the dataset: (later read these parameters from dataset.txt)
skipFactor = 200
downSample = 1
vrtNum = 264 #397 ###before upsampling
fiberNum = 111
isTrain = 0
yarnType = 'yarn8'



fn = "F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/yarnTypes/" + yarnType + "/datasets.txt"
with open(fn, 'r') as fin:
    info = fin.readline().split()
    dataset = info[0]
    firstFrame = int(info[1])
    lastFrame = int(info[2])
    yarn0 = int(info[3])
    yarn1 = int(info[4])


str1 = "yarnTypes/" + yarnType + "/config_step2.txt"
str2 = "yarnTypes/" + yarnType + "/datasets.txt"
########################## read input
print ("*************** phase0: READ INPUT ***************\n")
path = "F:/sandbox/fiberSimulation/dataSets/" + dataset+"/yarn/"
main (path, dataset, vrtNum, fiberNum, isTrain,  firstFrame, lastFrame, yarn0, yarn1, skipFactor, downSample)

########################## phase1
print ("*************** phase1: GENERATE NN INPUT ***************\n")
os.chdir('F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/x64/Release')
os.system('YarnGeneration 1 %s %s -w 5 -s %d -k %d -v %d -c 0 -rx 10 -ry 10 -rz 10 -rad 0.1' %(str1, str2, downSample, skipFactor, vrtNum))

########################## NN
print ("*************** phase1.5: NN ***************\n")
vrtx_us = vrtNum*downSample
main_NN(yarnType, downSample, dataset, firstFrame, lastFrame, yarn0, yarn1, skipFactor, vrtx_us )

########################## phase2
print ("*************** phase2: APPLY NN OUTPUT ***************\n")
os.system('YarnGeneration 2 %s %s -w 5 -s %d -k %d -v %d -c 1 -rx 10 -ry 10 -rz 10 -rad 0.1' %(str1, str2, downSample, skipFactor, vrtNum))
#os.system('YarnGeneration 2 %s %s -w 5 -s %d -k %d -v %d -c 0 -rx 10 -ry 10 -rz 10 -rad 0.1' %(str1, str2, downSample, skipFactor, vrtNum))



