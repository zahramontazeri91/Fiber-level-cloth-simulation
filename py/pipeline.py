''' open this file in python3.5 where tensorflow is installed
change the parameters for each dataset accordingly
Inputs are listed in yarnTypes/#yarnType/dataset '''

import sys
sys.path.insert(0, 'deforGrad')
from deform import main
from NN import main_NN
#from linearReg import main_NN
import os

# set paramters for the dataset: (later read these parameters from dataset.txt)
downSample = 2
vrtNum = 150 #397 ###before upsampling
isTrain = 0
trimPercent = 0.0 #larger than 0 if isTrain
yarnType = 'yarn4'

if (yarnType=='yarn4'):
    fiberNum = 160
elif (yarnType=='yarn8'):
    fiberNum = 111
elif (yarnType=='yarn11'):
    fiberNum = 120

# example: single_yarn/yarn4/stretch 0 49500 0 1 -k 500
fn = "F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/yarnTypes/" + yarnType + "/datasets.txt"
with open(fn, 'r') as fin:
    info = fin.readline().split()
    dataset = info[0]
    firstFrame = int(info[1])
    lastFrame = int(info[2])
    yarn0 = int(info[3])
    yarn1 = int(info[4])
    if (info[5] == '-v'):
        skipFactor = int(info[6])

str1 = "yarnTypes/" + yarnType + "/config_step2.txt"
str2 = "yarnTypes/" + yarnType + "/datasets.txt"

# In[]
########################## read input
print ("*************** phase0: READ INPUT ***************\n")
path = "F:/sandbox/fiberSimulation/dataSets/" + dataset+"/yarn/"
main (path, dataset, vrtNum, fiberNum, isTrain,  firstFrame, lastFrame, yarn0, yarn1, skipFactor, downSample)

# In[] 
########################## phase1
print ("*************** phase1: GENERATE NN INPUT ***************\n")
os.chdir('F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/x64/Release')
os.system('YarnGeneration 1 %s %s -w 5 -s %d -t %d -x %f -k %d -v %d -c 0 -rx 10 -ry 10 -rz 10 -rad 0.1' %(str1, str2, downSample, isTrain, trimPercent, skipFactor, vrtNum))

# In[]
########################## NN
print ("*************** phase1.5: NN ***************\n")
vrtx_us = vrtNum*downSample
main_NN(yarnType, downSample, dataset, firstFrame, lastFrame, yarn0, yarn1, skipFactor, vrtx_us )

# In[]
########################## phase2
print ("*************** phase2: APPLY NN OUTPUT ***************\n")
os.system('YarnGeneration 2 %s %s -w 5 -s %d -t %d -x %f -k %d -v %d -c 1 -rx 10 -ry 10 -rz 10 -rad 0.1' %(str1, str2, downSample, isTrain, trimPercent, skipFactor, vrtNum)) #deform the yarn
#os.system('YarnGeneration 2 %s %s -w 5 -s %d -t %d -x %f -k %d -v %d -c 0 -rx 10 -ry 10 -rz 10 -rad 0.1' %(str1, str2, downSample, isTrain, trimPercent, skipFactor, vrtNum)) #without deformation

# In[]
########################## phase2
os.chdir('F:/sandbox/fiberSimulation/yarn_generation_project/scene')

for i in range (int(firstFrame/skipFactor), int(lastFrame/skipFactor+1)):
    f = i*skipFactor
    for y in range (yarn0, yarn1):
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
        print(fn)
        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn))
        os.rename("fibers.exr", '../results/' + dataset + '/NN_bcsdf_' + str(f) + '_' + str(y) + '_highFrameRate.exr')
