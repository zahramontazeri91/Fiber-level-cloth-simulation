import numpy as np
import matplotlib
import matplotlib.pyplot as plt

w1 = 1.0
w2 = 1.0
w3 = 0.5/np.pi

def readData(fname):
    with open(fname, "r") as fin:
        n = int(fin.readline().strip())
        m = int(fin.readline().strip())
        data = np.zeros((n, m, 3))
        for i in range(0, n):
            if i > 0:
                assert int(fin.readline().strip()) == m
            for j in range(0, m):
                data[i, j, :] = [float(x) for x in fin.readline().strip().split()]
        return data


# Solving for a linear transformation to Q so that the result matches P.
# (i.e., P = A*Q + t)
def shape_matching(P, Q, allowTranslation = True):
    assert P.shape == Q.shape
    n = P.shape[1]

    if allowTranslation:
        p0 = np.mean(P, axis=1)
        q0 = np.mean(Q, axis=1)
    else:
        p0 = np.zeros((2, 1))
        q0 = np.zeros((2, 1))

    Apq = np.matrix(np.zeros((2, 2)))
    Aqq = np.matrix(np.zeros((2, 2)))
    for i in range(0, n):
        Apq += (P[:, i] - p0)*(Q[:, i] - q0).transpose()/n
        Aqq += (Q[:, i] - q0)*(Q[:, i] - q0).transpose()/n
    Aqq = np.linalg.inv(Aqq)
    A = Apq*Aqq

    U, sig, V = np.linalg.svd(A) # A == U*np.diag(sig)*V
    return U, sig, V, p0 - A*q0


def thetaDiff(theta0, theta1):
    diff = np.abs(np.mod(theta0 - theta1, 2.0*np.pi))
    diff = min(diff, 2.0*np.pi - diff)
    assert diff >= 0.0 and diff < 2.0*np.pi
    return diff


def cost(v0, v1):
    global w1, w2, w3
    assert len(v0) == 3 and len(v1) == 3
    return w1*np.power(v0[0] - v1[0], 2) + \
           w2*np.power(v0[1] - v1[1], 2) + \
           w3*np.power(thetaDiff(v0[2], v1[2]), 2)


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


if __name__ == "__main__":
    useTranslation = 1
    # useTranslation = 0

    baseName = "../frame00001"
#    baseName = "../genYarn_frame1"

    data0 = readData(baseName + "_scaled.txt")
    data1 = readData(baseName + "_compressed.txt")
    assert data0.shape == data1.shape
    m = data0.shape[1]

    S = np.zeros((2, m))
    theta = np.zeros((2, m))
    tansMag = np.zeros(m)
    idx = 100
    for i in range(0, m):
        A = np.matrix(data0[:, i, 0:2]).transpose()
        B = np.matrix(data1[:, i, 0:2]).transpose()

        U, sig, V, t0 = shape_matching(B, A, useTranslation)
        rot = U*V
        scl = V.transpose()*np.diag(sig)*V

        # Making sure V does not involve reflection
        if np.linalg.det(V) < 0.0:
            V = np.diag([1.0, -1.0])*V
        # Transpose V to match the convention of S = V*sig*V.transpose()
        V = V.transpose()

        assert np.linalg.det(rot) > 0.0
        assert np.allclose(V*np.diag(sig)*V.transpose(), scl)

        if i == 0:
            S[0, 0] = sig[0]
            S[1, 0] = sig[1]
            theta[0, 0] = np.arctan2(V[1, 0], V[0, 0])
        else:
            prev = np.array([S[0, i - 1], S[1, i - 1], theta[0, i - 1]])
            cur = np.zeros(3)
            best = np.inf
            for j in range(0, 4):
                cur[0 : 2] = sig if j % 2 == 0 else [sig[1], sig[0]]
                cur[2] = np.arctan2(V[1, 0], V[0, 0]) + j*0.5*np.pi
                if cur[2] >= np.pi:
                    cur[2] -= 2.0*np.pi
                val = cost(prev, cur)
                if best > val:
                    best = val
                    ans = np.copy(cur)
            S[:, i] = ans[0 : 2]
            theta[0, i] = ans[2]
        theta[1, i] = np.arctan2(rot[1, 0], rot[0, 0])
        tansMag[i] = np.linalg.norm(t0)

        if i == idx:
            P = np.asarray(B)
            Q0 = np.asarray(A)
            Q1 = np.asarray(rot*scl*A + t0)

    for i in range(1, m):
        for j in range(0, 2):
            theta[j, i] = nextTheta(theta[j, i - 1], theta[j, i])
    for i in range(0, 2):
        k = np.floor(0.25*(theta[i, 0] + theta[i, -1])/np.pi + 0.5)
        theta[i, :] -= k*2.0*np.pi

    plt.figure(figsize=(9, 8))

    plt.subplot(221)
    plt.plot(S[0, :], label='Sx')
    plt.plot(S[1, :], label='Sy')
    plt.legend()
    if useTranslation:
        plt.suptitle("%s, with translation" % baseName, fontsize=14)
    else:
        plt.suptitle("%s, without translation" % baseName, fontsize=14)

    plt.subplot(222)
    plt.plot(theta[0, :], label=r'$\theta_S$')
    plt.plot(theta[1, :], label=r'$\theta_R$')
    plt.legend()

    plt.subplot(223)
    plt.scatter(Q0[0, :], Q0[1, :], marker='.', alpha=0.25, label='Src')
    plt.scatter(Q1[0, :], Q1[1, :], marker='.', alpha=0.5, label='Src-def')
    plt.scatter(P[0, :], P[1, :], marker='.', alpha=0.5, label='Trg')
    plt.axis('equal')
    plt.axis([-0.05, 0.05, -0.05, 0.05])
    plt.legend()
    plt.title('Cross section #%d' % idx)

    plt.subplot(224)
    plt.plot(tansMag)
    plt.title('Amount of translation')

    # plt.tight_layout()
    plt.savefig('test%d.png' % useTranslation, dpi=200)
    # plt.show()
