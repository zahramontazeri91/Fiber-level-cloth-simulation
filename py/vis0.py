import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d')

data = np.loadtxt('../junk.txt')

n = data.shape[0]
L = 0.15

ax.plot(data[:, 0], data[:, 1], data[:, 2], color='black', linewidth=5.0)
for i in range(0, n):
    tang = data[i, 3 : 6]
    # tang /= np.linalg.norm(tang)
    norm = data[i, 6 : 9]
    # norm /= np.linalg.norm(norm)
    binorm = np.cross(norm, tang)
    # binorm /= np.linalg.norm(binorm)

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

plt.axis('equal')

plt.tight_layout()
plt.show()
