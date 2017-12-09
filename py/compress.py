import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from perlin import SimplexNoise

m = 0
with open("../genYarn_frame1.txt", "r") as fin:
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
            pos = np.array([float(x) for x in fin.readline().strip().split()])*0.25
            assert len(pos) == 4
            curCurve[i, :] = pos[0 : 3]
        avgCurve += curCurve
        curves.append(curCurve)
avgCurve /= m

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

with open("../genYarn_frame1_compressed.txt", "w") as fout:
    fout.write(str(m) + "\n")
    for i in range(0, m):
        fout.write(str(n) + "\n")
        for j in range(0, n):
            fout.write("%.6f %.6f %.6f\n" % (curves[i][j, 0], curves[i][j, 1], curves[i][j, 2]))


#print(m)
#for i in range(0, m):
#    print(n)
#    for j in range(0, n):
#        print("%.6f %.6f %.6f" % (curves[i][j, 0], curves[i][j, 1], curves[i][j, 2]))

# with open("compress_info.txt", "w") as fout:
#     fout.write(str(n) + "\n")
#     for i in range(0, n):
#         fout.write("%.6f %.6f %.6f\n" % (R1[i], R2[i], theta[i]))
#     fout.close()

plt.figure(figsize=(9, 4))

plt.subplot(121)
plt.plot(R1, label='Sx')
plt.plot(R2, label='Sy')
plt.legend()

plt.subplot(122)
plt.plot(theta, label=r'$\theta_S$')
plt.plot(thetaR, label=r'$\theta_R$')
plt.legend()

plt.tight_layout()
plt.savefig('ref.png', dpi=200)
plt.show()
