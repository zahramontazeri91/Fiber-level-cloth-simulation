# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 08:45:56 2018

@author: zahra
"""

import numpy as np
from keras.layers import Convolution1D, Dense, MaxPooling1D, Flatten
from keras.models import Sequential
import math

def buildModel(window_size, filter_length, nb_features, nb_filter):
    model = Sequential()

    
    model.add(Convolution1D(nb_filter=nb_filter, filter_length=filter_length, activation='relu', input_shape=(window_size, nb_features)))
    print('model output layer CNN1: ', model.output_shape)
    model.add(MaxPooling1D()) #default pool_size is 2
    print('model output layer MP1: ', model.output_shape)
    model.add(Convolution1D(nb_filter=nb_filter, filter_length=filter_length, activation='relu'))
    print('model output layer CNN2: ', model.output_shape)
    model.add(MaxPooling1D())
    print('model output layer MP2: ', model.output_shape)
#    model.add(Convolution1D(nb_filter=nb_filter, filter_length=filter_length, activation='relu'))
#    print('model output layer CNN3: ', model.output_shape)    
#    model.add(MaxPooling1D())   
#    print('model output layer MP3: ', model.output_shape)   
    model.add(Flatten())
    print('model output layer flatten: ', model.output_shape)
    model.add(Dense(4, activation='linear'))
    print('model output layer Dense: ', model.output_shape)
    
    model.compile(loss='mse', optimizer='adam', metrics=['mae'])
    # To perform (binary) classification instead:
    # model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['binary_accuracy'])

    return model

def loadData(fn_trainX, fn_trainY, nb_features, window_size):
    
    X_train_all = np.loadtxt(fn_trainX, delimiter=None)
    Y_train_all = np.loadtxt(fn_trainY, delimiter=None)
            
    print("TrainingX shape: ", X_train_all.shape)
    print("TrainingY shape: ", Y_train_all.shape)

    assert (X_train_all.shape[0] == Y_train_all.shape[0] )
    nb_instance = Y_train_all.shape[0]
    nb_output = Y_train_all.shape[1]
    assert (X_train_all.shape[1] == window_size*nb_features)
    
    split = 0.8
    nb_halfdata = round(nb_instance*split)
    all_train = np.concatenate((X_train_all,Y_train_all), axis=1) 
    np.random.shuffle(all_train)
    print("all_train shape: ", all_train.shape)
    
    X_train = all_train[0:nb_halfdata,0: window_size*nb_features]
    Y_train = all_train[0:nb_halfdata,-nb_output:]   
    X_valid = all_train[nb_halfdata:,0: window_size*nb_features]
    Y_valid = all_train[nb_halfdata:,-nb_output:] 
    
    X_train_ = X_train.reshape(nb_halfdata, window_size, nb_features)
    Y_train_ = Y_train.reshape(nb_halfdata, nb_output)
    X_valid_ = X_valid.reshape(nb_instance-nb_halfdata, window_size, nb_features)
    Y_valid_ = Y_valid.reshape(nb_instance-nb_halfdata, nb_output)

    print("trainX shape: ", X_train_.shape)
    print("trainY shape: ", Y_train_.shape)
    print("validX shape: ", X_valid_.shape)
    print("validY shape: ", Y_valid_.shape)
    
    return X_train_, Y_train_, X_valid_, Y_valid_, nb_instance, nb_output

def extrapolate(predicted, totalNum, nb_output, stride):
    total = np.zeros(totalNum*nb_output).reshape(totalNum,nb_output)
    assert(predicted.shape[1] == total.shape[1])
    # extrapolate predicted by stride_num
    vrtNum = predicted.shape[0] * stride
    predictExtr = np.zeros(vrtNum*nb_output).reshape(vrtNum,nb_output)
    fidx = np.linspace (1,vrtNum-stride-1, num=vrtNum/stride, dtype=int)

    for i in range(0,vrtNum):
        k = np.searchsorted(fidx,i)
        if k==0:
#            print(i)
            predictExtr[i] = predicted[0]
        elif k == int(vrtNum/stride):
            predictExtr[i] = predicted[int(vrtNum/stride) - 1 ]
#            print(i)
        else: 
            w = (fidx[k] - i)/(fidx[k] - fidx[k - 1])
            predictExtr[i] = w*predicted[k - 1] + (1.0 - w)*predicted[k]
#            print(i,k, w)
    
    r = math.floor (( totalNum - predictExtr.shape[0] )/2)
    total[r:r+predictExtr.shape[0], :] = predictExtr
    return total
    
def evaluate(nb_features, window_size):
    
    # load data
    fn_trainX = "D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/all/NN/trainX_all.txt"
    fn_trainY = "D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/all/NN/trainY_all.txt"
    X_train, Y_train, X_valid, Y_valid, nb_instance, nb_output = loadData(fn_trainX, fn_trainY, nb_features, window_size)
    
    # train network
    filter_length = 3
    nb_filter = 8
    model = buildModel(window_size, filter_length, nb_features, nb_filter)
    model.fit(X_train, Y_train, epochs=25, batch_size=2, validation_data=(X_valid, Y_valid))
    
    return nb_output, model            

nb_features = 9
window_size = 21                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
nb_output, model = evaluate(nb_features, window_size)

# predict test data
yarnNum = 1
skipFactor = 500        
firstFrame = 8500
totalNum = 300
dataset = 'spacing0.5x_00011'
path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
frame0 = int(firstFrame/skipFactor)
frame1 = int(16000/skipFactor + 1)
for i in range (frame0, frame1):
    f = i*skipFactor
    for y in range (0,yarnNum):
        X_test = np.loadtxt(path + "testX_" + str(f-firstFrame) + '_' + str(y) + ".txt",delimiter=None)
        X_test_ = X_test.reshape(X_test.shape[0], window_size, nb_features)
        Y_test_NN = model.predict(X_test_) 
        Y_test_NN_ = Y_test_NN.reshape(Y_test_NN.shape[0], Y_test_NN.shape[2])
        Y_test_NN_total = extrapolate(Y_test_NN_, totalNum, nb_output, 1)
        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
        np.savetxt(path + filename, Y_test_NN_total, fmt='%.6f', delimiter=' ')