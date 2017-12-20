# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 16:59:25 2017

@author: zahra
"""
from __future__ import print_function
import keras
from keras.datasets import mnist
from keras.layers import Dense, Flatten
from keras.layers import Conv2D, MaxPooling2D
from keras.models import Sequential
import matplotlib.pylab as plt

import numpy as np

batch_size = 128
epochs = 10

# input image dimensions
img_c, img_r = 9, 30
target_c, target_r = 2, 30
nb_data = 10
wnd_size = 3

# read data
data = np.zeros(img_r*img_c).reshape(img_c, img_r)
target = np.zeros(target_r*target_c)
path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/NN/'
nb_train = int(nb_data*0.7)
nb_test = int(nb_data - nb_data*0.7 )
x_train = np.zeros((nb_train, img_c, img_r ))
y_train = np.zeros((nb_train, target_c*target_r ))
x_test = np.zeros((nb_test, img_c, img_r ))
y_test = np.zeros((nb_test, target_c*target_r ))
for i in range (0,nb_data):
    fn = path + 'trainX_' + str(i) + '.txt'
    data = np.loadtxt(fn)

    fn = path + 'trainY_' + str(i) + '.txt'
    tmp = np.loadtxt(fn)
    target = tmp[:, 0:target_c]
    target = target.reshape(target_r*target_c)
    if i < nb_train : 
        x_train[i] = np.transpose(data)
        y_train[i] = target
    else:
        x_test[i-nb_train] = np.transpose(data)
        y_test[i-nb_train] = target       

#x_train = np.random.rand(nb_data, img_c, img_r)
#y_train = np.random.rand(nb_data, target_c * target_r)
#x_test = np.random.rand(nb_data, img_c, img_r)
#y_test = np.random.rand(nb_data, target_c * target_r)



# reshape the data into a 4D tensor - (sample_number, x_img_size, y_img_size, num_channels)
# because the MNIST is greyscale, we only have a single channel - RGB colour images would have 3
x_train = x_train.reshape(x_train.shape[0], img_c, img_r, 1)
x_test = x_test.reshape(x_test.shape[0], img_c, img_r, 1)
input_shape = (img_c, img_r, 1)

print('x_train shape:', x_train.shape)
print(x_train.shape[0], 'train samples')
print(x_test.shape[0], 'test samples')



model = Sequential()
layer1 = model.add(Conv2D(1, kernel_size=(wnd_size, img_c), strides=(1, 1),
                 activation='relu',
                 input_shape=input_shape))
layer2 = model.add(MaxPooling2D(pool_size=(1,3), strides=(1,3)))
layer3 = model.add(Flatten())
layer4 = model.add(Dense(target_c * target_r, activation='softmax'))

model.compile(optimizer='sgd', loss='mse', metrics=['accuracy'])



class AccuracyHistory(keras.callbacks.Callback):
    def on_train_begin(self, logs={}):
        self.acc = []

    def on_epoch_end(self, batch, logs={}):
        self.acc.append(logs.get('acc'))

history = AccuracyHistory()

model.fit(x_train, y_train,
          batch_size=batch_size,
          epochs=epochs,
          verbose=1,
          validation_data=(x_test, y_test),
          callbacks=[history]
          )
score = model.evaluate(x_test, y_test, verbose=0)
print('Test loss:', score[0])
print('Test accuracy:', score[1])
plt.plot(range(1, 11), history.acc)
plt.xlabel('Epochs')
plt.ylabel('Accuracy')
plt.show()