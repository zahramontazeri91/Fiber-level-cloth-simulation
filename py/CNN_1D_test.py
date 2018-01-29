# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 08:54:29 2018

@author: zahra
"""

### LOAD PACKAGES 
from numpy.random import seed
from pandas import read_csv, DataFrame
from sklearn.preprocessing import minmax_scale
from keras.layers.convolutional import Conv1D, MaxPooling1D
from keras.optimizers import SGD
from keras.models import Sequential
from keras.layers import Dense, Flatten
import numpy as np
### PREPARE THE DATA 
# nb_samples = 1000
# window_size = 3
# nb_channels (input) = 9
# nb_channels (output) = 4
Y = np.ones(1000*3*4)
y_train = Y.reshape(1000,3, 4)

X = np.zeros(1000*3*9)
x_train = X.reshape(1000,3, 9)


### FIT A 1D CONVOLUTIONAL NEURAL NETWORK
seed(2017)
conv = Sequential()
conv.add(Conv1D(32, 2, input_shape = x_train.shape[1:3], activation = 'relu', padding = 'same'))
conv.add(Conv1D(64, 2, input_shape = x_train.shape[1:3], activation = 'relu', padding = 'same'))
conv.add(Conv1D(4, 2, input_shape = x_train.shape[1:3], activation = 'relu', padding = 'same'))
#conv.add(MaxPooling1D(2))
#conv.add(Flatten())
#conv.add(Dense(1, activation = 'sigmoid'))
sgd = SGD(lr = 0.1, momentum = 0.9, decay = 0, nesterov = False)
conv.compile(loss = 'binary_crossentropy', optimizer = sgd, metrics = ['accuracy'])
conv.fit(x_train, y_train, batch_size = 1, epochs = 3, verbose = 0)


X = np.zeros(10*3*9)
x_test = X.reshape(10,3, 9)
y_test = conv.predict(x_test)
print(y_test)