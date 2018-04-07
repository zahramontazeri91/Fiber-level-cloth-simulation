# -*- coding: utf-8 -*-
"""
Created on Sat Mar 17 10:06:03 2018

@author: zahra
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
import numpy as np


path = '../input/'
skipFactor = 500
firstFrame = 8000
lastFrame = 17000
isFirst = 1
for i in range (firstFrame/skipFactor,lastFrame/skipFactor + 1):  
    f = i * skipFactor
    fn_train = path + 'pattern/yarn4/spacing1.0x/00011//NN/trainX_' + str(f) + '_0.txt'
    if (isFirst):
        X_train = np.loadtxt(fn_train)
        isFirst = 0
    else:
        X = np.loadtxt(fn_train)
        X_train = np.concatenate((X_train, X), axis=0)
    
    
skipFactor = 10000
firstFrame = 0
lastFrame = 98000
isFirst = 1
for i in range (firstFrame/skipFactor,lastFrame/skipFactor + 1):  
    f = i * skipFactor
#    fn_test = path + 'pattern/yarn4/spacing1.5x/10/NN/trainX_' + str(f) + '_0.txt'
    fn_test = path + 'twist/yarn4/damp2_500/NN/trainX_' + str(f) + '_0.txt'
#    fn_test = path + 'woven/yarn4/spacing1.0x/00011/NN/trainX_' + str(f) + '_0.txt'
#    fn_test = path + 'stretch/yarn4/stretch/NN/trainX_' + str(f) + '_0.txt'
    if (isFirst):
        X_test = np.loadtxt(fn_test)
        isFirst = 0
    else:
        X = np.loadtxt(fn_test)
        X_test = np.concatenate((X_test, X), axis=0)
        


# To getter a better understanding of interaction of the dimensions
# plot the first three PCA dimensions
fig = plt.figure(1, figsize=(8, 6))
ax = Axes3D(fig, elev=0, azim=0)

pca = PCA(n_components=3)
pca.fit(X_train)
X_reduced_train = pca.fit_transform(X_train)
#train = ax.scatter(X_reduced_train[:, 0], X_reduced_train[:, 1], X_reduced_train[:, 2],
#           cmap=plt.cm.Set1, edgecolor='k', s=60, label='train data', alpha=0.3)
X_reduced_test = pca.fit_transform(X_test)
test = ax.scatter(X_reduced_test[:, 0], X_reduced_test[:, 1], X_reduced_test[:, 2], c='r',
           cmap=plt.cm.Set1, edgecolor='k', s=70, label='test data - woven', alpha=0.2)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)

#ax.set_xlim3d(-1,1)
#ax.set_ylim3d(-1,1)
#ax.set_zlim3d(-1,1) 
ax.set_title("local-DG dimention reduction spacing1.0x, 00011")
ax.set_xlabel("1st eigenvector")
ax.w_xaxis.set_ticklabels([])
ax.set_ylabel("2nd eigenvector")
ax.w_yaxis.set_ticklabels([])
ax.set_zlabel("3rd eigenvector")
ax.w_zaxis.set_ticklabels([])

plt.savefig('../../data/test_damp2_500_space_side.jpg') 
plt.show()