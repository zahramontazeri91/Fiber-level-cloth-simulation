# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 11:30:51 2018

@author: zahra
"""

import matplotlib.pyplot as plt
import numpy as np

fn = 'F:/YarnGeneration/input/yarn11/single_teeth_1.6_10/globalRot_0_0.txt'
data = np.loadtxt(fn, delimiter=None)
#data = data.reshape(data.shape[0], 1)

plt.plot(data[:,0])
plt.plot(data[:,1])
plt.plot(data[:,2])
plt.plot(data[:,3])
plt.show()
    