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

    m = len(thetaS)
    for i in range(1, m):
        thetaS[i] = nextTheta(thetaS[i - 1], thetaS[i])
        thetaR[i] = nextTheta(thetaR[i - 1], thetaR[i])
    thetaS -= np.floor(0.25*(thetaS[0] + thetaS[-1])/np.pi + 0.5)*2.0*np.pi
    thetaR -= np.floor(0.25*(thetaR[0] + thetaR[-1])/np.pi + 0.5)*2.0*np.pi

    plt.figure(figsize=(8,10))

    plt.subplot(411)
    plt.plot(Sx, label='Sx')
    plt.plot(Sy, label='Sy')
    plt.legend()
    if title != "":
        plt.title(title, fontsize=14)
#    plt.ylim([0,1.1])

    plt.subplot(412)
    plt.plot(thetaS, label=r'$\theta_S$', color = 'g')
    plt.legend()

    plt.subplot(413)
    plt.plot(thetaR, label=r'$\theta_R$')
    plt.legend()

    plt.subplot(414)
    plt.plot(L2, label='L2 error', color = 'r')
    plt.legend()
#    plt.ylim([0,300])
    
    plt.tight_layout()
    if fname != "":
        plt.savefig("../../data/vis_crossSections/shapeMatch.png", dpi=200)
    if showWindow:
        plt.show()

#def visualize_T(Sx, Sy, thetaS, thetaR, Tx, Ty, title = "", fname = "", showWindow = True):
##    assert Sx.shape == Sy.shape and Sx.shape == thetaS.shape and Sx.shape == thetaR.shape
#
#    m = len(thetaS)
#    for i in range(1, m):
#        thetaS[i] = nextTheta(thetaS[i - 1], thetaS[i])
#        thetaR[i] = nextTheta(thetaR[i - 1], thetaR[i])
#    thetaS -= np.floor(0.25*(thetaS[0] + thetaS[-1])/np.pi + 0.5)*2.0*np.pi
#    thetaR -= np.floor(0.25*(thetaR[0] + thetaR[-1])/np.pi + 0.5)*2.0*np.pi
#
#    plt.figure(figsize=(8, 8))
#
#    plt.subplot(411)
#    plt.plot(Sx, label='Sx')
#    plt.plot(Sy, label='Sy')
#    plt.legend()
#    if title != "":
#        plt.title(title, fontsize=14)
##    plt.ylim([0,1.1])
#
#    plt.subplot(412)
#    plt.plot(thetaS, label=r'$\theta_S$')
#    plt.legend()
#
#    plt.subplot(413)
#    plt.plot(thetaR, label=r'$\theta_R$')
#    plt.legend()
#
#    plt.subplot(414)
#    plt.plot(Tx, label='Tx')
#    plt.plot(Ty, label='Ty')
#    plt.legend()
#    
#    plt.tight_layout()
#    if fname != "":
#        plt.savefig("../../data/vis_crossSections/shapeMatch.png", dpi=200)
#    if showWindow:
#        plt.show()
        
lng = []
shrt = []
theta = []
rot = []
Tx = []
Ty = []
L2 = []

#fname = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/x64/Release/compressParams.txt'
fname = '../compressParams_seg.txt'

with open(fname, 'r') as fin:
    N = int(fin.readline())
    for i in range (0,N):
        compress = fin.readline().split()
#        if i>65 and i<75:
        lng.append(float(compress[0]))
        shrt.append(float(compress[1]))
        theta.append(float(compress[2]))
        rot.append(float(compress[3]))
#         Tx.append(float(compress[4]))
#         Ty.append(float(compress[5]))
    fin.close()
 
with open('D:/sandbox/fiberSimulation/yarn_generation_project/data/L2.txt', 'r') as finL2:
    N = int(finL2.readline())
    for i in range (0,N):   
        e = finL2.readline()
        L2.append(float(e))
    finL2.close()

title = 'dataset3\n compressed: simulated frame 11400 \n reference: simulated frame 0'
#title = 'compressed: synthetic frame 1 compressed \n reference: synthetic frame 1'
visualize(lng, shrt, theta, rot, title, fname)

#ind = np.arange(N) 
#plt.figure(figsize=(15,10))
#
#plt.subplot(311)
#plt.plot([200,200], [0, 0.06], color='black') 
#plt.plot([1200,1200], [0, 0.06], color='black') 
##plt.ylim(-0.2, 2.5)
#plt.ylim(0.0, 2.5)
#plt.plot(ind, lng, color='r', label='Sx')
#plt.plot(ind, shrt, color='b', label='Sx')
#plt.legend()
#plt.title('Simulated data frame 29')
#          
#plt.subplot(312)         
#plt.plot([200,200], [0, 3], color='black') 
#plt.plot([1200,1200], [0, 3], color='black') 
#plt.plot(ind, theta, color='g', label=r'$\theta_S$')
#plt.legend()
#
#plt.subplot(313)         
#plt.plot([200,200], [0, 3], color='black') 
#plt.plot([1200,1200], [0, 3], color='black') 
#plt.plot(ind, rot, color='y', label=r'$\theta_R$')
#plt.legend()
#
#plt.savefig("../../data/vis_crossSections/shapeMatch.png")
#plt.show()