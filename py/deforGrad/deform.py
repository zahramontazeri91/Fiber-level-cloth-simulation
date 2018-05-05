import numpy as np
from fiber_bundle import FiberBundle
     
def rotationFromVectors(a, b):
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    if (s==0): #if it's already aligned with x-axis
        return np.identity(3)
    # TODO: for rotation 180 deg
    c = np.dot(a, b)
    V = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
    return np.identity(3) + V + np.dot(V, V)/(1.0 + c)

def transform(vrtNum, cntr_n_obj, cntr_n, twist_n, dg_n, internal_n, physicalParam_trans, def_obj, src_obj, isTrain):
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
#                if (abs(val[3]) > 10**-5):
                pre = val[3] 
                twist.append(pre)
                fout.writelines('%.8f %.8f %.8f \n' % (val[0]*0.25, val[1]*0.25, val[2]*0.25) )
                if len(pts1) == vrtNum:
                    break 
                
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
                
                
    
    if (isTrain==1):    
        with open(def_obj, "w") as fout:
            print >> fout, fiberNum
            for i in range(0, fiberNum):
                print >> fout, vrtNum*downSample
                for j in range(0, vrtNum*downSample):
                    v = bundle.fiber_vertex(i, j)
                    
#                    v = np.dot(R0, v) + t0  #Don't deform the yarn
#                    print >> fout, '%.6f %.6f %.6f' % (v[0], v[1], v[2])
                    print >> fout, '%.6f %.6f %.6f' % (v[0]*0.25, v[1]*0.25, v[2]*0.25)
    #FiberBundle(path + "frame_0015000fiber_00.obj").output_mitsuba("test0.txt")
 
# In[]:
###########################
def main (path, dataset, vrtNum, isTrain, firstFrame, lastFrame):
    print (dataset)
    for i in range (firstFrame/skipFactor,lastFrame/skipFactor + 1):  
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
            def_obj = "D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/data/" + dataset + '/simul_frame_' + str(f) + '_' + str(y) +'.txt'

            #output file:
            wrtPath = "D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/" + dataset
            physicalParam_trans = wrtPath + "/physicalParam/physical_"+str(f) + '_' + str(y)+"_world.txt"
            cntr_n = wrtPath + '/centerYarn_' + str(f) + '_' + str(y) + '_ds.txt'
            twist_n = wrtPath + '/twist_' + str(f) + '_' + str(y) + '_us.txt'
            transform(vrtNum, cntr_n_obj, cntr_n, twist_n, dg_n, internal_n, physicalParam_trans, def_obj, src_obj, isTrain)
# In[]:  
            
skipFactor = 1000
downSample = 1
vrtNum = 500 ###before upsampling
fiberNum = 160
yarn0 = 0
yarn1 = 200
isTrain = 0 #####
#dataset = 'pattern/yarn4/spacing0.5x/10/Raymond'
#dataset = 'woven/yarn4/spacing1.0x/00011'
#dataset = 'stretch/yarn4/stretch'
#dataset = 'fall/yarn4/fall'
#dataset = 'ball_fall'
#dataset = 'shear/shear0403'
#dataset = 'yarn_stretch'
#dataset = 'single_yarn/yarn4/teeth/4_1.6'
dataset = 'woven/release/yarn4/100x100'


path = "D:/sandbox/fiberSimulation/dataSets/" + dataset+"/yarn/"

restFrame = '0010000'
firstFrame = 0
lastFrame = 9500
main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)


#############################   
# In[]:            
#skipFactor = 1000
#downSample = 2
#vrtNum = 50
#fiberNum = 99 #160
#yarn0 = 0
#yarn1 = 16
#isTrain = 0
#firstFrame = 0
##dataset = 'speedtest/1.0_both'
##dataset = 'single_yarn/yarn11/stretch'
##dataset = 'single_yarn/yarn4/teeth/4_1.6'
##dataset = 'single_yarn/yarn9/teeth/4_1.2'
##dataset = 'single_yarn/yarn100'
#dataset = 'woven/release/yarn9/8x8'
#
#path = "D:/sandbox/fiberSimulation/dataSets/" + dataset+"/yarn/"
#lastFrame = 12000
#main (path, dataset, vrtNum, isTrain, firstFrame, lastFrame)
##            
#############################   
#            
#skipFactor = 100
#downSample = 1 ############### can change from 2 only if not the training data because of fiber-sim
#vrtNum = 150 #without upsample
#fiberNum = 160
#yarn0 = 0
#yarn1 = 1
#isTrain = 0
#firstFrame = 30000
##dataset = 'pattern/yarn4/spacing1.0x/00011'
#dataset = 'single_yarn/yarn4/stretch'
#path = "D:/sandbox/fiberSimulation/dataSets/" +dataset+ "/yarn/"
#lastFrame = 30000
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
yarn_type = 'pattern/yarn11'
fiberNum = 120
###########################

#dataset =  yarn_type + '/spacing0.5x/10'
#path = "D:/sandbox/fiberSimulation/dataSets/"  +'/'+dataset+"/yarn/"
#lastFrame = 14000
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing0.5x/00011'
#path = "D:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 16000
#main (path, dataset, vrtNum, isTrain, firstFrame, lastFrame)
##
#dataset =  yarn_type + '/spacing0.5x/10100'
#path = "D:/sandbox/fiberSimulation/dataSets/" + dataset+"/yarn/"
#lastFrame = 15000
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing0.5x/11110' 
#path = "D:/sandbox/fiberSimulation/dataSets/" + dataset+"/yarn/"
#lastFrame = 15000
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#########
#dataset =  yarn_type + '/spacing1.0x/10'
#path = "D:/sandbox/fiberSimulation/dataSets/" +'/'+dataset+"/yarn/"
#lastFrame = 14500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame) 
#
#dataset =  yarn_type + '/spacing1.0x/00011' 
#path = "D:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 17000
#main (path, dataset, vrtNum, isTrain, firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing1.0x/10100'
#path = "D:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 15500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing1.0x/11110'
#path = "D:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 16000
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)  
#
#########
#dataset =  yarn_type + '/spacing1.5x/10'
#path = "D:/sandbox/fiberSimulation/dataSets/" +'/'+dataset+"/yarn/"
#lastFrame = 15000
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing1.5x/00011'
#path = "D:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 17500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#
#dataset =  yarn_type + '/spacing1.5x/10100'
#path = "D:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 16000
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#
#dataset =  yarn_type + '/spacing1.5x/11110'
#path = "D:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 16500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#########
#dataset =  yarn_type + '/spacing2.0x/10'
#path = "D:/sandbox/fiberSimulation/dataSets/" +'/'+dataset+"/yarn/"
#lastFrame = 16000
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing2.0x/00011'
#path = "D:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 18000
#main (path, dataset, vrtNum, isTrain, firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing2.0x/10100'
#path = "D:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 17000
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing2.0x/11110'
#path = "D:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 17500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
###########
#dataset =  yarn_type + '/spacing2.5x/10'
#path = "D:/sandbox/fiberSimulation/dataSets/" +'/'+dataset+"/yarn/"
#lastFrame = 16500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing2.5x/00011'
#path = "D:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 18500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#
#dataset =  yarn_type + '/spacing2.5x/10100'
#path = "D:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 17500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#
#dataset =  yarn_type + '/spacing2.5x/11110'
#path = "D:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 18500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
###########
#
#dataset =  yarn_type + '/spacing3.0x/10'
#path = "D:/sandbox/fiberSimulation/dataSets/" +'/'+dataset+"/yarn/"
#lastFrame = 16500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing3.0x/00011'
#path = "D:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 19000
#main (path, dataset, vrtNum, isTrain, firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing3.0x/10100'
#path = "D:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 19500
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
#
#dataset =  yarn_type + '/spacing3.0x/11110'
#path = "D:/sandbox/fiberSimulation/dataSets/" +dataset+"/yarn/"
#lastFrame = 19000
#main (path, dataset, vrtNum, isTrain,  firstFrame, lastFrame)
