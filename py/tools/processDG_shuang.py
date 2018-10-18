import numpy as np
from fiber_bundle import FiberBundle

path = "C:/Users/zahra/Dropbox/fiber_data/training_data/train_teeth1231_ready/"

def rotationFromVectors(a, b):
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    assert s > 0
    c = np.dot(a, b)
    V = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
    return np.identity(3) + V + np.dot(V, V)/(1.0 + c)

pts0 = []
with open(path + "frame_0000000fiber_00.obj", "r") as fin:
    while True:
        line = fin.readline()
        if line[0] != 'v':
            break
        val = np.array([float(x) for x in line[1:].strip().split()])
        pts0.append(val[0 : 3])
        if len(pts0) == 100:
            break

pts1 = []
with open(path + "frame_0015000fiber_00.obj", "r") as fin:
    while True:
        line = fin.readline()
        if line[0] != 'v':
            break
        val = np.array([float(x) for x in line[1:].strip().split()])
        pts1.append(val[0 : 3])
        if len(pts1) == 100:
            break

dg = []
with open(path + "frame_0015000fiber_00.fe", "r") as fin:
    for line in fin.readlines():
        val = [float(x) for x in line.strip().split()]
        val = np.reshape(val, (3, 3))
        dg.append(val)

bundle = FiberBundle(path + "frame_0000000fiber_00.obj")
bundle.gen_central_line()
idx = np.array(bundle.downsample(3))

n = len(dg)
assert len(pts0) == n + 1 and len(idx) == n + 1

fidx = 0.5*(idx[0 : -1] + idx[1 :])
print(idx)
print(fidx)

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

with open("test1.txt", "w") as fout:
    print >> fout, bundle.numFibers
    for i in range(0, bundle.numFibers):
        print >> fout, bundle.fiberLen
        for j in range(0, bundle.fiberLen):
            k = np.searchsorted(fidx, j)
            if k == 0:
                R0 = R[0]
                t0 = t[0]
            elif k >= n:
                R0 = R[n - 1]
                t0 = t[n - 1]
            else:
                w = (fidx[k] - j)/(fidx[k] - fidx[k - 1])
                R0 = w*R[k - 1] + (1.0 - w)*R[k]
                t0 = w*t[k - 1] + (1.0 - w)*t[k]

            v = bundle.fiber_vertex(i, j)
            v = np.dot(R0, v) + t0
            print >> fout, '%.6f %.6f %.6f' % (v[0], v[1], v[2])

#FiberBundle(path + "frame_0015000fiber_00.obj").output_mitsuba("test0.txt")
