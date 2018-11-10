''' open this file in python3.5 where tensorflow is installed
change the parameters for each dataset accordingly
Inputs are listed in yarnTypes/#yarnType/dataset '''
import os
########################## set dataset parameters
# set paramters for the dataset: (later read these parameters from dataset.txt)
yarnType = 'yarn8'
vrtNum = 400 ###before upsampling
#150 #trainign
#397 #woven
#51 #for 6x6
isTrain = 0
trimPercent = 0.1 #larger than 0 if isTrain
upsample = 1
upsampleMore = 1 # 3 for teeth and 9 for stretch
hasCylinder = 0
pos = 0

########################## Read dataset.txt
# example in dataset.txt : single_yarn/yarn4/stretch 0 49500 0 1 -k 500
fn = "F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/yarnTypes/" + yarnType + "/datasets.txt"
with open(fn, 'r') as fin:  
    info = fin.readline().split()
    dataset = info[0]
    firstFrame = int(info[1])
    lastFrame = int(info[2])
    yarn0 = int(info[3])
    yarn1 = int(info[4])
    if (info[5] == '-k'):
        skipFactor = int(info[6])

# In[] 
########################## phase1
print ("*************** phase1: GENERATE NN INPUT ***************\n")
os.chdir('F:/YarnGeneration/x64/Release')
os.system('YarnGeneration 1 %s -w 5 -s %d -t %d -x %f -k %d -v %d' %(yarnType, upsample, isTrain, trimPercent, skipFactor, vrtNum))

# In[]print(type(tf.Session().run(tf.constant([1,2,3]))))
########################## NN
#from nonNN import main_NN
from NN import main_NN
print ("*************** phase1.5: NN ***************\n")
vrtx_us = vrtNum*upsample
main_NN(yarnType, upsample, dataset, firstFrame, lastFrame, yarn0, yarn1, skipFactor, vrtx_us )

# In[]
########################## phase2
print ("*************** phase2: APPLY NN OUTPUT ***************\n")
os.chdir('F:/YarnGeneration/x64/Release')
#os.system('YarnGeneration 2 %s -w 5 -s %d -s2 %d -t %d -x %f -k %d -v %d -c 1 -vol 0 -rx 10 -ry 10 -rz 10 -rad 0.1' %(yarnType, upsample, upsampleMore, isTrain, trimPercent, skipFactor, vrtNum)) #deform the yarn
#os.system('YarnGeneration 2 %s -w 5 -s %d -s2 %d -t %d -x %f -k %d -v %d -c 0 -vol 0 -rx 10 -ry 10 -rz 10 -rad 0.1' %(yarnType, upsample, upsampleMore, isTrain, trimPercent, skipFactor, vrtNum)) #without deformation

# In[]
########################## write mitsuba xml file
import sys
sys.path.insert(0, '../scene')
import generate_single

if hasCylinder:
    xmlfile = 'fibers.xml'
    spp = 32
    generate_single.genScene (xmlfile, spp, yarnType )
#    generate_single.genScene (xmlfile, spp )
    
    ########################## rendering
    def cubic_ease_function( t, t0, t1, L):
        ta=t0+(t1-t0)/3 
        tb=t1-(t1-t0)/3
        yh = (L * 2.0) / (t1 - t0 + tb - ta)
        if(t < t0 or t > t1):
            return 0.0
        else:
            if(t < ta):
                return (yh * (t0 - t) * (t0 - t) * (t0 - 3.0 * ta + 2.0 * t)) / ((t0 - ta) * (t0 - ta) * (t0 - ta))
            elif (t > tb):
                return (yh * (t1 - t) * (t1 - t) * (t1 - 3.0 * tb + 2.0 * t)) / ((t1 - tb) * (t1 - tb) * (t1 - tb))
            else:
                return yh
            
    t0 = 0.0
    t1 = 1.0
    L = 0.8/4
    dt = 0.00004 * skipFactor
    #dt = (t1-t0)/( (lastFrame-firstFrame)/skipFactor )
    pos = L
    t = 0.0
    c = 0
# In[]
############################ 
os.chdir('F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/scene')

for f in range (firstFrame, lastFrame+1, skipFactor):
    
    if hasCylinder:    
        t = c * dt
        c = c + 1 #frame counter 
        yh = cubic_ease_function( t, t0, t1, L)
        pos = pos - yh*dt
        print('frame: ', f, 'time: ', t,'position: ', pos)
    ####
    for y in range (yarn0, yarn1):
        
        h1 = pos
        h2 = -h1
          
#        fn = "../fibersim/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../output/" + dataset + "/genYarn_wo_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ "_us.txt"
#        fn = "../fibersim/slipstitchrib_2.3.txt"
        
        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D h1=%.8f -D h2=%.8f -D fn="%s" fibers.xml' % (h1, h2, fn))
        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba knitted_stretch.xml')
        os.rename("knitted_stretch.exr", '../../results/' + dataset + '/simul_' + str(f) + '_' + str(y) + '.exr')
        
    
