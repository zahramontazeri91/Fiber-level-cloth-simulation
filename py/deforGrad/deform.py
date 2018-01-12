import numpy as np
from fiber_bundle import FiberBundle
    
def rotationFromVectors(a, b):
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    assert s > 0
    c = np.dot(a, b)
    V = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
    return np.identity(3) + V + np.dot(V, V)/(1.0 + c)

def transform(vrtNum, cntr_0_obj, cntr_n_obj, cntr_n, dg_n, internal_n, physicalParam_trans, def_obj, src_obj):
    pts0 = []
    with open(cntr_0_obj, "r") as fin:
        while True:
            line = fin.readline()
            if line[0] != 'v':
                break
            val = np.array([float(x) for x in line[1:].strip().split()])
            pts0.append(val[0 : 3])
            if len(pts0) == vrtNum:
                break               
                
    pts1 = []
    with open(cntr_n_obj, "r") as fin:
        with open(cntr_n, "w") as fout:
            fout.writelines('%d \n' % (vrtNum) )
            while True:
                line = fin.readline()
                if line[0] != 'v':
                    break
                val = np.array([float(x) for x in line[1:].strip().split()])
                pts1.append(val[0 : 3])
                fout.writelines('%.8f %.8f %.8f \n' % (val[0]*0.25, val[1]*0.25, val[2]*0.25) )
                if len(pts1) == vrtNum:
                    break 
                
    dg = []
    with open( dg_n, "r") as fin:
        for line in fin.readlines():
            val = [float(x) for x in line.strip().split()]
            val = np.reshape(val, (3, 3))
            dg.append(val)
    n = len(dg)
        
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
            
    bundle = FiberBundle(src_obj)
    bundle.gen_central_line()
    idx = np.array(bundle.downsample(downSample))
    n1 = len(stretch)
    n2 = len(bend)    
    
    assert n1==n2 and n1==n+1
    assert len(pts0) == n + 1 and len(idx) == n + 1
    
    fidx = 0.5*(idx[0 : -1] + idx[1 :])
#    print(idx)
#    print(fidx)
    
    centers0 = [None]*n
    centers1 = [None]*n
    for i in range(0, n):
        centers0[i] = 0.5*(pts0[i] + pts0[i + 1])
        centers1[i] = 0.5*(pts1[i] + pts1[i + 1])
    
    R = [None]*n
    t = [None]*n
    
    # np.set_printoptions(formatter={"float" : lambda x : "%.6f" % x})
    for i in range(0, n):
        e = pts0[i + 1] - pts0[i]
        e /= np.linalg.norm(e)
        R[i] = np.dot(dg[i], rotationFromVectors(e, np.array([1.0, 0.0, 0.0])))
        t[i] = centers1[i] - np.dot(R[i], centers0[i])

                    
    with open(physicalParam_trans, "w") as fout1:
        with open(def_obj, "w") as fout:
            print >> fout, bundle.numFibers
            for i in range(0, bundle.numFibers):
                print >> fout, bundle.fiberLen
                for j in range(0, bundle.fiberLen):
                    k = np.searchsorted(fidx, j)
                    if k == 0:
                        R0 = R[0]
                        t0 = t[0]
                        S0 = stretch[0]
                        B0 = bend[0] 
                    elif k >= n:
                        R0 = R[n - 1]
                        t0 = t[n - 1]
                        S0 = stretch[n - 1]
                        B0 = bend[n - 1]                        
                    else:
                        w = (fidx[k] - j)/(fidx[k] - fidx[k - 1])
                        R0 = w*R[k - 1] + (1.0 - w)*R[k]
                        t0 = w*t[k - 1] + (1.0 - w)*t[k]                 
                        S0 = w*stretch[k - 1] + (1.0 - w)*stretch[k]
                        B0 = w*bend[k - 1] + (1.0 - w)*bend[k]
                    #####
                    if i==0:
                        fout1.writelines('%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f ' % (R0[0,0], R0[0,1], R0[0,2], R0[1,0], R0[1,1], R0[1,2], R0[2,0], R0[2,1], R0[2,2]) )
#                        fout1.writelines('%.8f %.8f %.8f' % (t0[0]*0.25, t0[1]*0.25, t0[2]*0.25 ) )
                        fout1.writelines('%.8f %.8f %.8f ' % (S0[0], S0[1], S0[2] ) )
                        fout1.writelines('%.8f %.8f %.8f' % (B0[0], B0[1], B0[2] ) )
                        fout1.writelines('\n')
                    ######
                    v = bundle.fiber_vertex(i, j)
#                    v = np.dot(R0, v) + t0  Don't deform the yarn
                    print >> fout, '%.6f %.6f %.6f' % (v[0]*0.25, v[1]*0.25, v[2]*0.25)

#FiberBundle(path + "frame_0015000fiber_00.obj").output_mitsuba("test0.txt")

# In[]:
#dataset = 'spacing3.0x_rotate_test'
dataset = 'spacing0.5x'
vrtNum = 150
downSample = 2
isTrain = 1
if (isTrain):
    path = "D:/sandbox/fiberSimulation/dataSets/spacing/train/"+dataset+"/yarn/"
else:
    path = "D:/sandbox/fiberSimulation/dataSets/spacing/test/"+dataset+"/yarn/"

skipFactor = 5
#for i in range (0,150/skipFactor + 1):
for i in range (1,2):  
    f = i * skipFactor
    f = 140
    if f < 10 :
        frameNum = '0000'+ str(f) + '00'
    elif f <100 :
        frameNum = '000'+ str(f) + '00' 
    else:
       frameNum = '00'+ str(f) + '00' 
    #input file: (Change for testData)
    if (isTrain):
        src_obj = path + '../fiber/frame_' + frameNum + 'fiber_00.obj'
    else:
        src_obj = path + '../fiber/frame_0000000fiber_00.obj'
    cntr_0_obj = path + 'frame_0000000fiber_00.obj'
    cntr_n_obj = path + 'frame_' + frameNum + 'fiber_00.obj'
    cntr_n = path + 'frame_' + frameNum + 'fiber_00.obj'
    dg_n = path + 'frame_' + frameNum + 'fiber_00.fe' 
    internal_n = path + 'frame_' + frameNum + 'fiber_00.sforce' 
    def_obj = "D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/data/" + dataset + '/simul_frame_' + str(f) + '_0.txt'
    
    #output file:
    wrtPath = "D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/" + dataset
    physicalParam_trans = wrtPath + "/physicalParam/physical_"+str(f)+"_world.txt"
    cntr_n = wrtPath + '/centerYarn_' + str(f) + '_ds.txt'
    
    transform(vrtNum, cntr_0_obj, cntr_n_obj, cntr_n, dg_n, internal_n, physicalParam_trans, def_obj, src_obj)