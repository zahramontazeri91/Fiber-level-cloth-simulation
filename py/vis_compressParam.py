"""
visualize compress parameter

@author: zahra
"""
import numpy as np
import matplotlib

from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'font.size': 16})
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import math


def nextTheta(theta0, theta1):
    k = np.floor((theta0 - theta1)*0.5/np.pi)
    ans = theta1 + (k - 1)*2.0*np.pi
    best = np.abs(theta0 - ans)
    for i in range(0, 2):
        cur = theta1 + (k + i)*2.0*np.pi
        val = np.abs(theta0 - cur)
        if best > val:
            best = val;
            ans = cur
    return ans

def visualize(Sx, Sy, thetaS, thetaR, title = "", fname = "", showWindow = True):
#    assert Sx.shape == Sy.shape and Sx.shape == thetaS.shape and Sx.shape == thetaR.shape

#    m = len(thetaS)
#    for i in range(1, m):
#        thetaS[i] = nextTheta(thetaS[i - 1], thetaS[i])
##        thetaR[i] = nextTheta(thetaR[i - 1], thetaR[i])
#    thetaS -= np.floor(0.25*(thetaS[0] + thetaS[-1])/np.pi + 0.5)*2.0*np.pi
##    thetaR -= np.floor(0.25*(thetaR[0] + thetaR[-1])/np.pi + 0.5)*2.0*np.pi
    
    
    plt.figure(figsize=(8,10))

    plt.subplot(211)
    plt.plot(Sx, label='S(0,0)')
    plt.plot(Sy, label='S(1,1)')
    plt.legend()
    if title != "":
        plt.title(title, fontsize=14)
#    plt.ylim([0,1.1])

    plt.subplot(212)
    plt.plot(thetaS, label=r'S(0,1)', color = 'g')
    plt.legend()
    
#    plt.subplot(413)
#    plt.plot(thetaR, label=r'$\theta_R$')
#    plt.legend()
#
#    plt.subplot(414)
#    plt.plot(L2, label='L2 error', color = 'r')
#    plt.legend()
#    plt.ylim([0,0.1])
    
    plt.tight_layout()
    if fname != "":
        plt.savefig('../../data/vis_crossSections/'+ fname + '.png', dpi=200)
    if showWindow:
        plt.show()

       


path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/NN/'

for w in range (0,1):
#for w in range (0,36):
    lng = []
    shrt = []
    theta = []
    rot = []
    Tx = []
    Ty = []
    L2 = []

    f = w*5

#    if f < 10 :
#        fname = path + 'param_1220_0000' + str(f) +'00.txt'
#    elif f < 100 :
#        fname = path + 'param_1220_000' + str(f) +'00.txt'
#    else :
#        fname = path + 'param_1220_00' + str(f) +'00.txt'
#    fname = path + 'trainY_' + str(f) + '.txt'    
    fname = path + 'testY_NN.txt'
    #fname = '../compressParams.txt'
    
    writefile = 'frame_' + str(f) 
    with open(fname, 'r') as fin:
    #    N = int(fin.readline())
#        N = 299
        N = 122
        for i in range (0,N):
            compress = fin.readline().split()
    #        if i>65 and i<75:
            lng.append(float(compress[0]))
            shrt.append(float(compress[1]))
            theta.append(float(compress[2]))
    #        rot.append(float(compress[3]))
    
        fin.close()
     
    with open('D:/sandbox/fiberSimulation/yarn_generation_project/data/L2.txt', 'r') as finL2:
        N2 = int(finL2.readline())
        for i in range (0,N2):   
            e = finL2.readline()
            L2.append(float(e))
        finL2.close()
    
    title = 'dataset 1120\n compressed: simulated frame ' + str(f) +' \n reference: simulated frame 0'
    visualize(lng, shrt, theta, rot, title, writefile)

##################################################
#ind = np.arange(N) 
#plt.figure(figsize=(15,10))
#
#plt.subplot(311)
#plt.plot(ind, lng, color='r', label='Sx')
#plt.plot(ind, shrt, color='b', label='Sx')
#plt.legend()
#plt.title('Simulated data frame 29')
#          
#plt.subplot(312)         
#plt.plot(ind, theta, color='g', label=r'$\theta_S$')
#plt.legend()
#
#plt.subplot(313)         
#plt.plot(ind, rot, color='y', label=r'$\theta_R$')
#plt.legend()
#
#plt.savefig("../../data/vis_crossSections/shapeMatch.png")
#plt.show()