
import numpy as np
import matplotlib

from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import math


def visualize(S00, S01, S10, S11, S00_NN, S01_NN, S10_NN, S11_NN, title = "", showWindow = True):
    
    plt.figure(figsize=(8,10))

    plt.subplot(411)
    plt.plot(S00, label='S(0,0) train')
    plt.plot(S00_NN, label='S(0,0) NN', color='r')
#    plt.legend()
    plt.ylim([0,6.28])
    
    plt.subplot(412)
    plt.plot(S01, label='S(0,1) train', color='b')
    plt.plot(S01_NN, label='S(0,1) NN', color='r')
    plt.legend()
    
    plt.subplot(413)
    plt.plot(S10, label='S(1,0) train', color='y')
    plt.plot(S10_NN, label='S(1,0) NN', color='r')
    plt.legend()
    
    plt.subplot(414)
    plt.plot(S11, label='S(1,1) train', color='black')
    plt.plot(S11_NN, label='S(1,1) NN', color='r')
    plt.legend()
    if title != "":
        plt.title(title, fontsize=14)
#    plt.ylim([0,1.1])
    
    plt.tight_layout()
    if showWindow:
        plt.show()

       

################################

S00 = []
S01 = []
S10 = []
S11 = []
#fname = '../input/spacing1.0x_00011/NN/testY_17000_0.txt'
fname = '../input/spacing1.0x_00011_woven/NN/testY_200_0.txt'
data = np.loadtxt(fname)
N = data.shape[0]
with open(fname, 'r') as fin:
    for i in range (0,N):
        compress = fin.readline().split()

        S00.append(float(compress[0]))
#        S01.append(float(compress[1]))
#        S10.append(float(compress[2]))
#        S11.append(float(compress[3]))
    fin.close()

S00_NN = []
S01_NN = []
S10_NN = []
S11_NN = []    
#fname = '../input/spacing1.0x_00011/NN/testY_NN_full_17000_0.txt'
fname = '../input/spacing1.0x_00011_woven/NN/testY_NN_full_200_0.txt'
data = np.loadtxt(fname)
N = data.shape[0]
with open(fname, 'r') as fin:
    for i in range (0,N):
        compress = fin.readline().split()

#        S00_NN.append(float(compress[0]))
#        S01_NN.append(float(compress[1]))
#        S10_NN.append(float(compress[2]))
#        S11_NN.append(float(compress[3]))
    fin.close()

title = 'NN output'
visualize(S00, S01, S10, S11, S00_NN, S01_NN, S10_NN, S11_NN, title)