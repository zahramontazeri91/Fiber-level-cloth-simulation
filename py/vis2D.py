
import numpy as np
import matplotlib

from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import math


def visualize(S00, S01, S10, S11, title = "", showWindow = True):
    
    plt.figure(figsize=(8,10))

#    plt.subplot(211)
    plt.plot(S00, label='S(0,0)')
    plt.plot(S01, label='S(0,1)')
    plt.plot(S10, label='S(1,0)')
    plt.plot(S11, label='S(1,1)')
    plt.legend()
    if title != "":
        plt.title(title, fontsize=14)
#    plt.ylim([0,1.1])

#    plt.subplot(212)
#    plt.plot(S10, label='S(0,0)')
#    plt.plot(S11, label='S(1,1)')
#    plt.legend()

    
    plt.tight_layout()
    if showWindow:
        plt.show()

       

################################

S00 = []
S01 = []
S10 = []
S11 = []
fname = '../input/spacing1.0x_00011/NN/trainY_17000_0.txt'

with open(fname, 'r') as fin:
#    N = int(fin.readline())
    N = 211
    for i in range (0,N):
        compress = fin.readline().split()

        S00.append(float(compress[0]))
        S01.append(float(compress[1]))
        S10.append(float(compress[2]))
        S11.append(float(compress[3]))

    fin.close()

title = 'NN output'
visualize(S00, S01, S10, S11, title)