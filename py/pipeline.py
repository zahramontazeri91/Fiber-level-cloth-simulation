''' open this file in python3.5 where tensorflow is installed
change the parameters for each dataset accordingly
Inputs are listed in yarnTypes/#yarnType/dataset '''
import os
########################## set dataset parameters
#info = 'yarn100/train/spacing1.0x/00011 7000 17000 0 1 -k 500 -v 150 -t 1 -s 2'
#info = 'yarn4/single_stretch 30000 30000 0 1 -k 500 -t 1 -v 150 -s 2 -ss 5'
#info = 'yarn8/single_teeth_1.2_00110 500 1000 0 1 -k 500 -t 1 -v 150 -s 2 -ss 5'
#info = 'yarn11/single_teeth_1.6_10 0 45500 0 1 -k 500 -t 1 -v 150 -s 2 -ss 5'

info = 'yarn8/slipstitchrib_push 14000 23000 0 1 -k 1000 -t 0 -v 85832 -s 1 -vol 1'
#info = 'yarn8/slipstitchrib_stretch 0 49500 0 1 -k 500 -t 0 -v 3372 -s 1 -vol 1'
#info = 'yarn8/slipstitchrib_stretch 0 49500 0 1 -k 1000 -t 0 -v 85832 -s 1 -vol 1'

########################## Read info
yarnType = info.split("/")[0]
split = info.split()
dataset = split[0]
firstFrame = int(split[1])
lastFrame = int(split[2])
yarn0 = int(split[3])
yarn1 = int(split[4])
isTrain = 0
skipFactor = 500
vrtNum = 150 #before upsampling
upsample = 2
n = len(split)
for i in range (5, n-1): #5 first terms already reserved
    if (split[i] == '-t'):
        isTrain = int(split[i + 1])
    if (split[i] == '-k'):
        skipFactor = int(split[i + 1])
    if (split[i] == '-v'):
        vrtNum = int(split[i + 1])
    if (split[i] == '-s'):
        upsample = int(split[i + 1])

# In[] 
########################## phase1
print ("*************** phase1: GENERATE NN INPUT ***************\n")
os.chdir('F:/YarnGeneration/x64/Release')
os.system('YarnGeneration 1 %s' %(info))

# In[]
########################## NN\
#from nonNN import main_NN
from NN import main_NN
print ("*************** phase1.5: NN ***************\n")
vrtx_us = vrtNum*upsample
main_NN(yarnType, upsample, dataset, firstFrame, lastFrame, yarn0, yarn1, skipFactor, vrtx_us )

# In[]
########################## phase2
print ("*************** phase2: APPLY NN OUTPUT ***************\n")
os.chdir('F:/YarnGeneration/x64/Release')
os.system('YarnGeneration 2 %s -c 1' %(info)) #prints both compressed and non-compressed resutls

# In[]
# use scene/render_ct2.py for volume rendering (use python2.7)

# In[]
########################## write mitsuba xml file

os.chdir('F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/scene')
import generate_single
import generate_tiled
hasCylinder = 1
xmlfile = 'fibers.xml'
spp = 16
isFlat = 0
#generate_single.genScene (xmlfile, spp)
#generate_single.genScene (xmlfile, spp, yarnType )
#generate_tiled.genScene (xmlfile, spp, isFlat, yarnType)
  
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
pos = L
c = 2 #start simulation at dt (Raymond: "The simulation used t+2dt as the parameter currently.")

    
for f in range (firstFrame, lastFrame+1, skipFactor):
    
    if hasCylinder:    
        t = c * dt
        c = c + 1 #frame counter 
        yh = cubic_ease_function( t, t0, t1, L)
        pos = pos - yh*dt

    for y in range (yarn0, yarn1):
        h1 = pos
        h2 = -h1
        
#        fn = "../output/" + dataset + "/curve_" + str(f).zfill(7) + "_" + str(y).zfill(2) + ".txt"
#        fn = "../fibersim/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../output/" + dataset + "/genYarn_wo_" + str(f) + "_" + str(y)+ ".txt"
        fn = "../output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ "_us.txt"
#        fn = "../fibersim/slipstitchrib_2.3.txt"
        
        print(fn)
        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D h1=%.8f -D h2=%.8f -D fn="%s" fibers.xml' % (h1, h2, fn))
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba knitted_stretch.xml')
        os.rename("fibers.exr", '../../results/' + dataset + '/test_path_tile_' + str(f) + '_' + str(y) + '.exr')
        
  