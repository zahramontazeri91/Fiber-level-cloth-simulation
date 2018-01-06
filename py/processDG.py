import numpy as np

def printNormalized(a):
    b = a/np.linalg.norm(a)
    print(b)


def rotationFromVectors(a, b):
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    assert s > 0
    c = np.dot(a, b)
    V = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
    return np.identity(3) + V + np.dot(V, V)/(1.0 + c)


dataset = '1231'
wrt_path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/' + dataset
f = 150 
vrtNum = 300

pts0 = []
with open(wrt_path + '/centerYarn_' + str(0) + '.txt', "r") as fin:
    while True:
        line = fin.readline()
#        if line[0] != 'v':
#            break
        val = np.array([float(x) for x in line[0:].strip().split()])
        pts0.append(val[0 : 3])
        if len(pts0) == vrtNum:
            break

pts1 = []
with open(wrt_path + '/centerYarn_' + str(f) + '.txt', "r") as fin:
    while True:
        line = fin.readline()
#        if line[0] != 'v':
#            break
        val = np.array([float(x) for x in line[0:].strip().split()])
        pts1.append(val[0 : 3])
        if len(pts1) == vrtNum:
            break

dg = []
with open(wrt_path + '/physicalParam/physical_' + str(f) + '_world.txt', "r") as fin:
    for line in fin.readlines():
        val = [float(x) for x in line.strip().split()[0:9]]
        val = np.reshape(val, (3, 3))
        dg.append(val)

print(len(dg))
n = len(dg)
assert len(pts0) == n 

#centers0 = [None]*n
#centers1 = [None]*n
#for i in range(0, n):
#    centers0[i] = 0.5*(pts0[i] + pts0[i + 1])
#    centers1[i] = 0.5*(pts1[i] + pts1[i + 1])


np.set_printoptions(formatter={'float':lambda x: '%.6f' % x})


f = 150
with open(wrt_path + '/deformGrad_' + str(f) + '.txt', 'w') as fout:
    for i in range (0,vrtNum) :
        e = pts0[i + 1] - pts0[i]
        e /= np.linalg.norm(e)
        R = rotationFromVectors(e, np.array([1.0, 0.0, 0.0]))
        R = np.dot(dg[i], R)
        t = pts1[i] - np.dot(R, pts0[i])
        fout.writelines('%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f ' % (R[0,0], R[0,1], R[0,2], R[1,0], R[1,1], R[1,2], R[2,0], R[2,1], R[2,2]) )
        fout.writelines('%.6f %.6f %.6f\n' %(t[0], t[1], t[2]))
fout.close()

#    print(pts1[i])
#    print(np.dot(R, pts0[i]) + t)
#    
#    print(pts1[i + 1])
#    print(np.dot(R, pts0[i + 1]) + t)