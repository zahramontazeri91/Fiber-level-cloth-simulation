''' open this file in python3.5 where tensorflow is installed
change the parameters for each dataset accordingly
Make sure yarnTypes/dataset matches with the dataset here '''


import sys
sys.path.insert(0, 'deforGrad')
from deform import main
from NN import main_NN
import os


# set paramters for the dataset: (later read these parameters from dataset.txt)
#skipFactor = 500
#downSample = 2
#vrtNum = 51
#fiberNum = 111 #160
#yarn0 = 0
#yarn1 = 12
#isTrain = 0
#firstFrame = 6000
#lastFrame = 6000
#dataset = 'woven/6x6'
########################
skipFactor = 200
downSample = 1
vrtNum = 397 ###before upsampling
fiberNum = 111
yarn0 = 150
yarn1 = 151
isTrain = 0
dataset = 'woven/arbitrary_pattern/150x100'
firstFrame = 400
lastFrame = 400



########################## read input
path = "F:/sandbox/fiberSimulation/dataSets/" + dataset+"/yarn/"
main (path, dataset, vrtNum, fiberNum, isTrain,  firstFrame, lastFrame, yarn0, yarn1, skipFactor, downSample)

########################## phase1
os.chdir('F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/x64/Release')
os.system('YarnGeneration 1 yarnTypes/yarn8/config_step2.txt yarnTypes/yarn8/datasets.txt -w 5 -s %d -k %d -v %d -c 0 -rx 10 -ry 10 -rz 10 -rad 0.1' %(downSample, skipFactor, vrtNum))

########################## NN
vrtx_us = vrtNum*downSample
main_NN('yarn8',downSample, dataset, firstFrame, lastFrame, yarn0, yarn1, skipFactor, vrtx_us )

########################## phase2
os.system('YarnGeneration 2 yarnTypes/yarn8/config_step2.txt yarnTypes/yarn8/datasets.txt -w 5 -s %d -k %d -v %d -c 0 -rx 10 -ry 10 -rz 10 -rad 0.1' %(downSample, skipFactor, vrtNum))


