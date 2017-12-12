import numpy as np
import matplotlib
import matplotlib.pyplot as plt

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
    assert Sx.shape == Sy.shape and Sx.shape == thetaS.shape and Sx.shape == thetaR.shape

    m = len(thetaS)
    for i in range(1, m):
        thetaS[i] = nextTheta(thetaS[i - 1], thetaS[i])
        thetaR[i] = nextTheta(thetaR[i - 1], thetaR[i])
    thetaS -= np.floor(0.25*(thetaS[0] + thetaS[-1])/np.pi + 0.5)*2.0*np.pi
    thetaR -= np.floor(0.25*(thetaR[0] + thetaR[-1])/np.pi + 0.5)*2.0*np.pi

    plt.figure(figsize=(8, 8))

    plt.subplot(311)
    plt.plot(Sx, label='Sx')
    plt.plot(Sy, label='Sy')
    plt.legend()
    if title != "":
        plt.title(title, fontsize=14)

    plt.subplot(312)
    plt.plot(thetaS, label=r'$\theta_S$')
    plt.legend()

    plt.subplot(313)
    plt.plot(thetaR, label=r'$\theta_R$')
    plt.legend()

    plt.tight_layout()
    if fname != "":
        plt.savefig(fname, dpi=200)
    if showWindow:
        plt.show()
