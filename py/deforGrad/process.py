import numpy as np
import io

m = 0
with io.open("frame00001_hairs.txt", "r") as fin:
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

print(n)
for i in range(0, n):
    print("%.6f %.6f %.6f" % (avgCurve[i, 0], avgCurve[i, 1], avgCurve[i, 2]))
