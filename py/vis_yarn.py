import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect('equal')

lower = np.array([1e10, 1e10, 1e10])
upper = np.array([-1e10, -1e10, -1e10])

with open('../gen_yarn.txt', 'r') as fin:
    n = int(fin.readline())
    for i in range(0, n):
        m = int(fin.readline())
        fiber = np.zeros((m, 3))
        for j in range(0, m):
            pos = [float(val) for val in fin.readline().strip().split(' ')]
            fiber[j, :] = pos

            lower = np.minimum(lower, pos)
            upper = np.maximum(upper, pos)

        ax.plot(fiber[:, 0], fiber[:, 1], fiber[:, 2], alpha=0.5)

span = np.max(upper - lower)
lower -= 0.05*span
upper += 0.05*span
mid = 0.5*(lower + upper)

# Create cubic bounding box to simulate equal aspect ratio
max_range = np.max(upper - lower)
Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + mid[0]
Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + mid[1]
Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + mid[2]
# Comment or uncomment following both lines to test the fake bounding box:
for xb, yb, zb in zip(Xb, Yb, Zb):
   ax.plot([xb], [yb], [zb], 'w')

plt.tight_layout()
plt.show()