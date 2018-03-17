import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})
import matplotlib.patches as mpatches

fig = plt.figure(figsize=(10,20))
ax = fig.gca(projection='3d')
#ax.set_aspect('equal')

# data = np.loadtxt('../junk_single.txt')
data = np.loadtxt('../../data/TNB.txt')
#cntr.x cntr.y cntr.z tg.x tg.y tg.z norm.x norm.y norm.z

n = data.shape[0]
L = 0.25

ax.plot(data[:, 0], data[:, 1], data[:, 2], color='black', linewidth=1.0)
for i in range(0, n):
    tang = data[i, 3 : 6]
    tang /= np.linalg.norm(tang)
    norm = data[i, 6 : 9]
    norm /= np.linalg.norm(norm)
    binorm = np.cross(tang, norm)
    binorm /= np.linalg.norm(binorm)
#    assert (np.dot(binorm, norm) < 0.0001 and np.dot(binorm, tang) < 0.0001 and np.dot(tang, norm)<0.001)
#    print(np.dot(binorm, norm), np.dot(binorm, tang), np.dot(tang, norm))

    ax.plot([data[i, 0], data[i, 0] + L*tang[0]], \
            [data[i, 1], data[i, 1] + L*tang[1]], \
            [data[i, 2], data[i, 2] + L*tang[2]], \
            color='red', alpha=0.5)
    ax.plot([data[i, 0], data[i, 0] + L*norm[0]], \
            [data[i, 1], data[i, 1] + L*norm[1]], \
            [data[i, 2], data[i, 2] + L*norm[2]], \
            color='green', alpha=0.5)
    ax.plot([data[i, 0], data[i, 0] + L*binorm[0]], \
            [data[i, 1], data[i, 1] + L*binorm[1]], \
            [data[i, 2], data[i, 2] + L*binorm[2]], \
            color='blue', alpha=0.5)

lower = np.amin(data[:, 0:3], axis=0) - 2.0*L
upper = np.amax(data[:, 0:3], axis=0) + 2.0*L
mid = 0.5*(lower + upper)

# Create cubic bounding box to simulate equal aspect ratio
max_range = np.max(upper - lower)
Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + mid[0]
Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + mid[1]
Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + mid[2]
# Comment or uncomment following both lines to test the fake bounding box:
for xb, yb, zb in zip(Xb, Yb, Zb):
   ax.plot([xb], [yb], [zb], 'w')

#plt.tight_layout()
red = mpatches.Patch(color='red', label='tangents')
blue = mpatches.Patch(color='blue', label='binormals')
green = mpatches.Patch(color='green', label='normals')

plt.legend(handles=[red, green, blue],prop={'size': 15}) 
ax.view_init(azim=10., elev=50)
plt.xlim(-.5,.5)
plt.ylim(-.5,.5)
plt.savefig("../../data/vis_crossSections/frame29_hermite_2.png" )
plt.show()
