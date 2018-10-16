import numpy as np
from fiber_bundle import FiberBundle

#activate writing centerlines in mitsuba format to use that for rendering mitsuba-ct2 with large yarn-radius instead of rendering fibers
mitsuba_centerline = 0
hasTwist = 1

def rotationFromVectors(a, b):
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    if (s==0): #if it's already aligned with x-axis
        return np.identity(3)
    # TODO: for rotation 180 deg
    c = np.dot(a, b)
    V = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
    return np.identity(3) + V + np.dot(V, V)/(1.0 + c)

def transform(vrtNum, fiberNum, cntr_n_obj, cntr_n, cntr_n_mitsuba, twist_n, dg_n, internal_n, physicalParam_trans, def_obj, src_obj, isTrain, downSample):
#    pts0 = []
#    with open(cntr_0_obj, "r") as fin:
#        while True:
#            line = fin.readline()
#            if line[0] != 'v':
#                break
#            val = np.array([float(x) for x in line[1:].strip().split()])
#            pts0.append(val[0 : 3])
#            if len(pts0) == vrtNum:
#                break               
                
    pts1 = []
    twist = []
    with open(cntr_n_obj, "r") as fin:
        with open(cntr_n, "w") as fout:
            pre = 0
            fout.writelines('%d \n' % (vrtNum) )
            while True:
                line = fin.readline()
                if line[0] != 'v':
                    break
                val = np.array([float(x) for x in line[1:].strip().split()]) 
                pts1.append(val[0 : 3])
                if (hasTwist):
                    pre = val[3] 
                twist.append(pre)
                fout.writelines('%.8f %.8f %.8f \n' % (val[0]*0.25, val[1]*0.25, val[2]*0.25) )
                if len(pts1) == vrtNum:
                    break 
#    print(pts1[1][0])
############# print centerlines in mitsuba format
    if (mitsuba_centerline):            
        with open(cntr_n_mitsuba, "w") as fout:
            fout.writelines('%d \n' % (1) ) #1 curve
            fout.writelines('%d \n' % (len(pts1)) )
            for i in range (0,len(pts1)):
                fout.writelines('%.8f %.8f %.8f \n' % (pts1[i][0]*0.25, pts1[i][1]*0.25, pts1[i][2]*0.25) )
#################
    dg = []
    with open( dg_n, "r") as fin:
#        for line in fin.readlines():
        for i in range (0, len(pts1)-1): #-1 because DG is one less than num of particles
            line = fin.readline()
            val = [float(x) for x in line.strip().split()]
            val = np.reshape(val, (3, 3))
            dg.append(val)
    n = len(dg)
        
    hasInternal = 0
    if (hasInternal==1):
        stretch = []
        bend = []
        with open( internal_n, "r") as fin:
            for i in range(0,n+1):
                line = fin.readline()
                val = [float(x) for x in line.strip().split()]
                val_s = np.reshape(val[0:3], (3,1))
                val_b = np.reshape(val[3:], (3,1))
                stretch.append(val_s)
                bend.append(val_b)
        n1 = len(stretch)
        n2 = len(bend)  
        assert n1==n2 and n1==n+1
            
    if (isTrain==1):        
        bundle = FiberBundle(src_obj)
        bundle.gen_central_line()
#        idx = np.array(bundle.downsample(downSample))

    idx = np.arange(0,vrtNum*downSample,downSample)
  
    
    if (len(pts1) != n + 1 or len(idx) != n + 1):
        print (len(pts1), len(idx), n + 1 )
    assert len(pts1) == n + 1 and len(idx) == n + 1
    
    fidx = 0.5*(idx[0 : -1] + idx[1 :])
#    print(idx)
#    print(fidx)
    
#    centers0 = [None]*n
#    centers1 = [None]*n
#    for i in range(0, n):
#        centers0[i] = 0.5*(pts0[i] + pts0[i + 1])
#        centers1[i] = 0.5*(pts1[i] + pts1[i + 1])
    
    R = [None]*n
#    t = [None]*n
    
    # np.set_printoptions(formatter={"float" : lambda x : "%.6f" % x})
    for i in range(0, n):
#        e = pts0[i + 1] - pts0[i]
#        e /= np.linalg.norm(e)
#        R[i] = np.dot(dg[i], rotationFromVectors(e, np.array([1.0, 0.0, 0.0])))
        R[i] = dg[i]
#        t[i] = centers1[i] - np.dot(R[i], centers0[i])

    fiber_length = vrtNum *downSample                
    with open(physicalParam_trans, "w") as fout1:
        with open(twist_n, "w") as fout_twist:
            fout_twist.writelines('%d \n' %fiber_length)
            for j in range(0, fiber_length):
                k = np.searchsorted(fidx, j)
                if k == 0:
                    theta = twist[0]
                    R0 = R[0]
#                    t0 = t[0]
                    if hasInternal==1:
                        S0 = stretch[0]
                        B0 = bend[0] 
                elif k >= n:
                    theta = twist[n-1]
                    R0 = R[n - 1]
#                    t0 = t[n - 1]
                    if hasInternal==1:
                        S0 = stretch[n - 1]
                        B0 = bend[n - 1]                        
                else:
                    w = (fidx[k] - j)/(fidx[k] - fidx[k - 1])
                    theta = w*twist[k - 1] + (1.0 - w)*twist[k]
                    R0 = w*R[k - 1] + (1.0 - w)*R[k]
#                    t0 = w*t[k - 1] + (1.0 - w)*t[k]       
                    if hasInternal==1:
                        S0 = w*stretch[k - 1] + (1.0 - w)*stretch[k]
                        B0 = w*bend[k - 1] + (1.0 - w)*bend[k]
                #####
                fout1.writelines('%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f ' % (R0[0,0], R0[0,1], R0[0,2], R0[1,0], R0[1,1], R0[1,2], R0[2,0], R0[2,1], R0[2,2]) )
#                        fout1.writelines('%.8f %.8f %.8f' % (t0[0]*0.25, t0[1]*0.25, t0[2]*0.25 ) )
                if hasInternal==1:
                    fout1.writelines('%.8f %.8f %.8f ' % (S0[0], S0[1], S0[2] ) )
                    fout1.writelines('%.8f %.8f %.8f' % (B0[0], B0[1], B0[2] ) )
                fout1.writelines('\n')
                
                fout_twist.writelines('%.8f \n' %theta)
                
##############################################
#    downSample_more = 4
#    fiber_length = vrtNum *downSample_more
#    with open(twist_n, "w") as fout_twist:
#        fout_twist.writelines('%d \n' %fiber_length)
#        for j in range(0, fiber_length):
#            k = np.searchsorted(fidx, j)
#            print(k)
#            if k == 0:
#                theta = twist[0]
#            elif k >= n:
#                theta = twist[n-1]                   
#            else:
#                w = (fidx[k] - j)/(fidx[k] - fidx[k - 1])
#                theta = w*twist[k - 1] + (1.0 - w)*twist[k]
          
#            fout_twist.writelines('%.8f \n' %theta)
##############################################
                
    if (isTrain==1):    
        with open(def_obj, "w") as fout:
            fout.writelines('%d \n' %fiberNum) ## python3: print >> fout, fiberNum
            for i in range(0, fiberNum):
                fout.writelines('%d\n' %(vrtNum*downSample) )  ## python3: print >> fout, vrtNum*downSample
                for j in range(0, vrtNum*downSample):
                    v = bundle.fiber_vertex(i, j)                    
#                    v = np.dot(R0, v) + t0  #Don't deform the yarn
#                    print >> fout, '%.6f %.6f %.6f' % (v[0], v[1], v[2])
                    fout.writelines('%.6f %.6f %.6f \n' % (v[0]*0.25, v[1]*0.25, v[2]*0.25) ) ## python3: print >> fout, '%.6f %.6f %.6f' % (v[0]*0.25, v[1]*0.25, v[2]*0.25)
    #FiberBundle(path + "frame_0015000fiber_00.obj").output_mitsuba("test0.txt")
 
# In[]:
###########################
def main (path, dataset, vrtNum, fiberNum, isTrain, firstFrame, lastFrame, yarn0, yarn1, skipFactor, downSample):
    print (dataset)
    print('frames:', firstFrame/skipFactor, lastFrame/skipFactor + 1)
    for i in range (int(firstFrame/skipFactor), int(lastFrame/skipFactor) + 1):  
        f = i * skipFactor
        print(f)
        frameNum = 'frame_'+str(f).zfill(7)
        for y in range(yarn0,yarn1):
            yarnNum = 'fiber_' + str(y).zfill(2)

            #input file: (Change for testData)
            if (isTrain):
                src_obj = path + '../fiber/' + frameNum + yarnNum+'.obj'
            else:
                src_obj = ''
            
            cntr_n_obj = path + frameNum + yarnNum + '.obj'
            dg_n = path + frameNum + yarnNum + '.fe' 
            internal_n = path + frameNum + yarnNum + '.sforce'
            def_obj = "F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/fibersim/" + dataset + '/simul_frame_' + str(f) + '_' + str(y) +'.txt'

            #output file:
            wrtPath = "F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/" + dataset
            physicalParam_trans = wrtPath + "/physicalParam/physical_"+str(f) + '_' + str(y)+"_world.txt"
            cntr_n = wrtPath + '/centerYarn_' + str(f) + '_' + str(y) + '_ds.txt'
            twist_n = wrtPath + '/twist_' + str(f) + '_' + str(y) + '_us.txt'
            cntr_n_mitsuba = "F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/output/" + dataset + '/centerlines/curve_' + str(f) + '_' + str(y) + '_ds_mitsuba.txt'
            
#            print (dg_n)
            transform(vrtNum, fiberNum, cntr_n_obj, cntr_n, cntr_n_mitsuba, twist_n, dg_n, internal_n, physicalParam_trans, def_obj, src_obj, isTrain, downSample)
# In[]: 
# set paramters for the dataset:
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
#
########################### read input
#path = "F:/sandbox/fiberSimulation/dataSets/" + dataset+"/yarn/"
#main (path, dataset, vrtNum, isTrain, firstFrame, lastFrame)



# In[]:  
            
#skipFactor = 200
#downSample = 1
#vrtNum = 397 ###before upsampling
#fiberNum = 111
#yarn0 = 150
#yarn1 = 250
#isTrain = 0 #####
##dataset = 'pattern/yarn4/spacing0.5x/10/Raymond'
##dataset = 'woven/yarn4/spacing1.0x/00011'
##dataset = 'stretch/yarn4/stretch'
##dataset = 'fall/yarn4/fall'
##dataset = 'ball_fall'
##dataset = 'shear/shear0403'
##dataset = 'yarn_stretch'
##dataset = 'single_yarn/yarn4/teeth/4_1.6'
##dataset = 'woven/push/yarn8/100x100'
##dataset = 'woven/stretch/yarn4/100x100'
##dataset = 'woven/6x6'
##dataset = 'woven/arbitrary_pattern/100x100'
#dataset = 'woven/arbitrary_pattern/150x100'
##dataset = 'woven/arbitrary_pattern/512x512'
#
#path = "F:/sandbox/fiberSimulation/dataSets/" + dataset+"/yarn/"
#
#restFrame = '0010000'
#firstFrame = 500
#lastFrame = 500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame, yarn0, yarn1, skipFactor)

# In[]: 
#skipFactor = 1000
#downSample = 1
#vrtNum = 3907 #2507 ###before upsampling
#fiberNum = 111
#yarn0 = 600
#yarn1 = 1000 #first generate for horizontal yarns with 2507 vertices) then generate for vertical for 3907 vertices
#isTrain = 0 #####
#dataset = 'woven/arbitrary_pattern/600x400'
#
#path = "F:/sandbox/fiberSimulation/dataSets/" + dataset+"/yarn/"
#
#restFrame = '0010000'
#firstFrame = 0
#lastFrame = 0
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)

# In[]:
#skipFactor = 500
#downSample = 1
#vrtNum = 1046 ###before upsampling
#fiberNum = 111
#yarn0 = 0
#yarn1 = 200
#isTrain = 0 #####
##dataset = 'pattern/yarn4/spacing0.5x/10/Raymond'
##dataset = 'woven/yarn4/spacing1.0x/00011'
##dataset = 'stretch/yarn4/stretch'
##dataset = 'fall/yarn4/fall'
##dataset = 'ball_fall'
##dataset = 'shear/shear0403'
##dataset = 'yarn_stretch'
##dataset = 'single_yarn/yarn4/teeth/4_1.6'
##dataset = 'woven/push/yarn8/100x100'
##dataset = 'woven/stretch/yarn4/100x100'
##dataset = 'woven/6x6'
#dataset = 'woven/arbitrary_pattern/100x100'
#
#path = "F:/sandbox/fiberSimulation/dataSets/" + dataset+"/yarn/"
#
#restFrame = '0010000'
#firstFrame = 0
#
#lastFrame = 15000
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)


#############################   
# In[]:            
#skipFactor = 500
#downSample = 2
#vrtNum = 51
#fiberNum = 111 #160
#yarn0 = 0
#yarn1 = 12
#isTrain = 0
#firstFrame = 6000
##dataset = 'speedtest/1.0_both'
##dataset = 'single_yarn/yarn100'
#dataset = 'woven/6x6'
#
#path = "F:/sandbox/fiberSimulation/dataSets/" + dataset+"/yarn/"
#lastFrame = 6000
#main (path, dataset, vrtNum, isTrain, firstFrame, lastFrame)
#            
#############################   
##            
#skipFactor = 1000
#downSample = 1 ############### can change from 2 only if not the training data because of fiber-sim
#vrtNum = 34607 #without upsample
#fiberNum = 111
#yarn0 = 0
#yarn1 = 1
#isTrain = 0
#firstFrame = 0
##dataset =  'single_yarn/yarn8/teeth/4_1.2_00110'
##dataset = 'pattern/yarn8/spacing1.0x/00011'
##dataset = 'single_yarn/yarn4/stretch'
##dataset = 'single_yarn/yarn11/teeth/4_1.6'
##dataset = 'single_yarn/yarn11/teeth/4_1.2'
##dataset = 'pattern/yarn4/spacing1.0x/00011'
##dataset = 'single_yarn/yarn9/stretch'
##dataset = 'single_yarn/yarn4/teeth/1.2_110'
#dataset = 'woven/knitted'
#path = "F:/sandbox/fiberSimulation/dataSets/" +dataset+ "/yarn/"
#lastFrame = 40000
#main (path, dataset, vrtNum, isTrain, firstFrame, lastFrame)
##########################
#skipFactor = 1000
#downSample = 2 ############### can change from 2 only if not the training data because of fiber-sim
#vrtNum = 150 #without upsample
#fiberNum = 120
#yarn0 = 0
#yarn1 = 1
#isTrain = 0
#firstFrame = 20000
##dataset =  'single_yarn/yarn8/teeth/4_1.2_00110'
##dataset = 'pattern/yarn8/spacing1.0x/00011'
##dataset = 'single_yarn/yarn4/stretch'
#dataset = 'single_yarn/yarn11/teeth/4_1.6'
##dataset = 'single_yarn/yarn11/teeth/4_1.2'
##dataset = 'pattern/yarn4/spacing1.0x/00011'
##dataset = 'single_yarn/yarn9/stretch'
##dataset = 'single_yarn/yarn11/teeth/1.2_110'
##dataset = 'woven/knitted'
#path = "F:/sandbox/fiberSimulation/dataSets/" +dataset+ "/yarn/"
#lastFrame = 20000
#main (path, dataset, vrtNum, isTrain, firstFrame, lastFrame)


############################
# In[]: 
skipFactor = 500
downSample = 2 ############
vrtNum = 150
yarn0 = 0
yarn1 = 1
isTrain = 1
firstFrame = 7000  
yarn_type = 'pattern/yarn100'
fiberNum = 150
###########################

#dataset =  yarn_type + '/spacing0.5x/10'
#path = "F:/sandbox/fiberSimulation/dataSets/"  +'/'+dataset+"/yarn/"
#lastFrame = 18500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing0.5x/00011'
#path = "F:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 20500
#main (path, dataset, vrtNum, isTrain, firstFrame, lastFrame)
##
#dataset =  yarn_type + '/spacing0.5x/10100'
#path = "F:/sandbox/fiberSimulation/dataSets/" + dataset+"/yarn/"
#lastFrame = 18000
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing0.5x/11110' 
#path = "F:/sandbox/fiberSimulation/dataSets/" + dataset+"/yarn/"
#lastFrame = 20500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#########
#dataset =  yarn_type + '/spacing1.0x/10'
#path = "F:/sandbox/fiberSimulation/dataSets/" +'/'+dataset+"/yarn/"
#lastFrame = 19500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame) 
#
#dataset =  yarn_type + '/spacing1.0x/00011' 
#path = "F:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 21500
#main (path, dataset, vrtNum, isTrain, firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing1.0x/10100'
#path = "F:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 19000
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing1.0x/11110'
#path = "F:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 21500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)  
#
#########
#dataset =  yarn_type + '/spacing1.5x/10'
#path = "F:/sandbox/fiberSimulation/dataSets/" +'/'+dataset+"/yarn/"
#lastFrame = 20500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing1.5x/00011'
#path = "F:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 22500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#
#dataset =  yarn_type + '/spacing1.5x/10100'
#path = "F:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 20000
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#
#dataset =  yarn_type + '/spacing1.5x/11110'
#path = "F:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 22500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#########
#dataset =  yarn_type + '/spacing2.0x/10'
#path = "F:/sandbox/fiberSimulation/dataSets/" +'/'+dataset+"/yarn/"
#lastFrame = 21500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing2.0x/00011'
#path = "F:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 23500
#main (path, dataset, vrtNum, isTrain, firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing2.0x/10100'
#path = "F:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 21000
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing2.0x/11110'
#path = "F:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 23500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
###########
#dataset =  yarn_type + '/spacing2.5x/10'
#path = "F:/sandbox/fiberSimulation/dataSets/" +'/'+dataset+"/yarn/"
#lastFrame = 22500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing2.5x/00011'
#path = "F:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 24500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#
#dataset =  yarn_type + '/spacing2.5x/10100'
#path = "F:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 22000
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#
#dataset =  yarn_type + '/spacing2.5x/11110'
#path = "F:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 24500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
###########
#
#dataset =  yarn_type + '/spacing3.0x/10'
#path = "F:/sandbox/fiberSimulation/dataSets/" +'/'+dataset+"/yarn/"
#lastFrame = 23500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing3.0x/00011'
#path = "F:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 25500
#main (path, dataset, vrtNum, isTrain, firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing3.0x/10100'
#path = "F:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 23000
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing3.0x/11110'
#path = "F:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 18000
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
