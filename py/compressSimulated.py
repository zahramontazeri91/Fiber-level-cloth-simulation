import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from perlin import SimplexNoise

m = 0
#with open("frame00001_hairs.txt", "r") as fin:
with open("../genYarn_frame1.txt", "r") as fin:

    line = fin.readline()
    m=0
    while True: 
        line = fin.readline()
        if not line:
            break
        m += 1
        if m == 1:
            n = int(line.strip())
            avgCurve = np.zeros((n, 3))
            curves = []
        else:
            assert n == int(line.strip())
        curCurve = np.zeros((n, 3))
        for i in range(0, n):
            pos = np.array([float(x) for x in fin.readline().strip().split()])
#            assert len(pos) == 4
            curCurve[i, :] = pos[0 : 3]
        avgCurve += curCurve
        curves.append(curCurve)
avgCurve /= m

# print(n)
# for i in range(0, n):
#     print("%.6f %.6f %.6f" % (avgCurve[i, 0], avgCurve[i, 1], avgCurve[i, 2]))

noise = SimplexNoise()
theta = np.linspace(0.0, 4.0*np.pi, n)
thetaR = np.linspace(0.0, -2.0*np.pi, n)
R1 = np.linspace(1.25, 0.75, n)
R2 = np.linspace(0.75, 1.25, n)

for i in range(0, n):
    offset = 0.5*noise.noise2(0.75*theta[i], 0.0)
    theta[i] += offset

    offset = 0.1*noise.noise2(5.0*R1[i], 10.0)
    R1[i] += offset

    offset = 0.1*noise.noise2(5.0*R2[i], 20.0)
    R2[i] += offset

    offset = noise.noise2(0.5*thetaR[i], 0.0)
    thetaR[i] += offset

with open("../compress_info.txt", "w") as fout:
    fout.write(str(n) + "\n")
    for i in range(0, n):
        fout.write("%.6f %.6f %.6f\n" % (R1[i], R2[i], theta[i]))
    fout.close()

for i in range(0, m):
    for j in range(0, n):
        offset = curves[i][j] - avgCurve[j]
        offset[2] = 0.0
        axisX = np.array([np.cos(theta[j]), np.sin(theta[j]), 0.0])
        axisY = np.array([-np.sin(theta[j]), np.cos(theta[j]), 0.0])
        x = np.dot(offset, axisX)*R1[j]
        y = np.dot(offset, axisY)*R2[j]
        offset = x*axisX + y*axisY

        rot = np.array([[np.cos(thetaR[j]), -np.sin(thetaR[j])], [np.sin(thetaR[j]), np.cos(thetaR[j])]])
        offset[0 : 2] = np.reshape(rot.dot([[offset[0]], [offset[1]]]), 2)

        curves[i][j] = avgCurve[j] + offset


with open("../genYarn_frame1_compressed_R.txt", "w") as fout:
    fout.write(str(m) + "\n")
    for i in range(0, m):
        fout.write(str(n) + "\n")
        for j in range(0, n):
            fout.write("%.6f %.6f %.6f\n" % (curves[i][j, 0], curves[i][j, 1], curves[i][j, 2]))
