# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 08:45:56 2018

@author: zahra
"""

import matplotlib.pyplot as plt
import numpy as np
from keras.layers import Convolution1D, Dense, MaxPooling1D, Flatten
from keras.models import Sequential
import math
from shutil import copyfile

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
    model.add(Convolution1D(nb_filter=nb_filter, filter_length=filter_length, activation='relu'))
    print('model output layer CNN3: ', model.output_shape)    
    model.add(MaxPooling1D())   
    print('model output layer MP3: ', model.output_shape)  
    
    
    model.add(Flatten())
    print('model output layer flatten: ', model.output_shape)
    model.add(Dense(4, activation='linear'))
    print('model output layer Dense: ', model.output_shape)
    
    model.compile(loss='mse', optimizer='adam', metrics=['mse'])
    # To perform (binary) classification instead:
    # model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['binary_accuracy'])

    return model

def loadData(fn_trainX, fn_trainY, nb_features):
    
    X_train_all = np.loadtxt(fn_trainX, delimiter=None)
    Y_train_all = np.loadtxt(fn_trainY, delimiter=None)
            
    print("TrainingX shape: ", X_train_all.shape)
    print("TrainingY shape: ", Y_train_all.shape)

    assert (X_train_all.shape[0] == Y_train_all.shape[0] )
    nb_instance = Y_train_all.shape[0]
    nb_output = Y_train_all.shape[1]
    window_size = int(X_train_all.shape[1] / nb_features )
#    if (X_train_all.shape[1] != window_size*nb_features):
#        print (X_train_all.shape[1], window_size*nb_features)
#    assert (X_train_all.shape[1] == window_size*nb_features)
    
    split = 0.99
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
    
    return X_train_, Y_train_, X_valid_, Y_valid_, nb_instance, nb_output, window_size

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
    
def evaluate(nb_features, fn_trainX, fn_trainY):
    
    # load data
    X_train, Y_train, X_valid, Y_valid, nb_instance, nb_output, window_size = loadData(fn_trainX, fn_trainY, nb_features)
    
    # train network
    filter_length = 3
    nb_filter = 256
    epochs=50
    model = buildModel(window_size, filter_length, nb_features, nb_filter)
    history = model.fit(X_train, Y_train, epochs=epochs, batch_size=2, validation_data=(X_valid, Y_valid))
    
    # Plot loss trajectory throughout training.
    plt.figure()
    plt.subplot(1,2,1)
    plt.plot(history.history['mean_squared_error'], label='train')
    plt.plot(history.history['val_mean_squared_error'], label='valid')
    plt.xlabel('Epoch')
    plt.ylabel('mse')
    plt.legend()
    plt.show()
    
    return nb_output, model, window_size            

## rotate the NN output back to RM frames
# In[]:
def rotate(predicted, angles):
    predicted_rot = predicted
    n = len(angles)
    for i in range (0,n):
#        c, s = np.cos(-1*angles[i]), np.sin(-1*angles[i])
        # should be angle not -angle because R_T * S * R 
        c, s = np.cos(angles[i]), np.sin(angles[i])
        R = np.matrix('{} {}; {} {}'.format(c, -s, s, c))
        S = predicted[i].reshape([2,2])
#        rot = R*S*R.transpose()
        rot = R.transpose()*S*R
        predicted_rot[i] = np.array([rot[0,0] , rot[0,1] , rot[1,0] , rot[1,1] ])
    return predicted_rot

# In[]:
def append2sets(dataset1, dataset2, w_path):
    path1 = w_path + dataset1 + '/NN/'
    path2 = w_path + dataset2 + '/NN/'
    X_train_all_1 = np.loadtxt(path1 + "trainX_all.txt",delimiter=None)
    Y_train_all_1 = np.loadtxt(path1 + "trainY_all.txt",delimiter=None)
    X_train_all_2 = np.loadtxt(path2 + "trainX_all.txt",delimiter=None)
    Y_train_all_2 = np.loadtxt(path2 + "trainY_all.txt",delimiter=None)
    
    X_train_all = np.concatenate((X_train_all_1,X_train_all_2), axis=0)
    Y_train_all = np.concatenate((Y_train_all_1,Y_train_all_2), axis=0)
    
    np.savetxt(w_path + 'all/NN/trainX_all.txt', X_train_all, fmt='%.6f', delimiter=' ')
    np.savetxt(w_path + 'all/NN/trainY_all.txt', Y_train_all, fmt='%.6f', delimiter=' ')
    
## append training data from different datasets
# In[] 
def appendTrainingData(datasets, w_path, fn_trainX, fn_trainY):
    for i in range (0, len(datasets)):
        if (i==0):
            p0x = w_path + datasets[i] + '/NN/trainX_all.txt'
            p0y = w_path + datasets[i] + '/NN/trainY_all.txt'
            copyfile (p0x, fn_trainX)
            copyfile (p0y, fn_trainY)
        else:
            append2sets('all', datasets[i],w_path)
        
# In[] 
datasets = []
w_path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'
#datasets.append('spacing0.5x_00011')
#datasets.append('spacing0.5x_10100')
#datasets.append('spacing0.5x_11110')
datasets.append('spacing1.0x_00011')
#datasets.append('spacing1.0x_10100')
#datasets.append('spacing1.0x_11110')
#datasets.append('spacing1.5x_00011')
#datasets.append('spacing1.5x_10100')
#datasets.append('spacing1.5x_11110')

fn_trainX = w_path + "all/NN/trainX_all.txt"
fn_trainY = w_path + "all/NN/trainY_all.txt"
appendTrainingData(datasets, w_path, fn_trainX, fn_trainY)

nb_features = 9                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
nb_output, model, window_size = evaluate(nb_features, fn_trainX, fn_trainY)

# In[]:
# predict test data
yarnNum = 1
skipFactor = 500        
firstFrame = 17000
lastFrame = 17000
totalNum = 300 ################# NOTE: downsampled
dataset = 'yarn4/spacing1.0x/00011'
path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
frame0 = int(firstFrame/skipFactor)
frame1 = int(lastFrame/skipFactor + 1)
for i in range (frame0, frame1):
    f = i*skipFactor
    for y in range (0,yarnNum):
        X_test = np.loadtxt(path + "trainX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
        X_test_ = X_test.reshape(X_test.shape[0], window_size, nb_features)
        Y_test_NN = model.predict(X_test_) 
#        Y_test_NN = np.loadtxt(path + "trainY_" + str(f) + '_' + str(y) + ".txt", delimiter=None)
        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
        angles = np.loadtxt(anglesFile, delimiter=None)
        Y_test_NN_rot = rotate(Y_test_NN, angles)
#        Y_test_NN_rot = Y_test_NN
        Y_test_NN_total = extrapolate(Y_test_NN_rot, totalNum, nb_output, 1)
        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
        np.savetxt(path + filename, Y_test_NN_total, fmt='%.6f', delimiter=' ')
#####################################
        
# In[]:
# predict test data
#yarnNum = 1
#skipFactor = 1        
#firstFrame = 240
#lastFrame = 240
#totalNum = 300 ################# NOTE: downsampled
#dataset = 'twist_only'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "trainX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        X_test_ = X_test.reshape(X_test.shape[0], window_size, nb_features)
#        Y_test_NN = model.predict(X_test_) 
##        Y_test_NN = np.loadtxt(path + "trainY_" + str(f) + '_' + str(y) + ".txt", delimiter=None)
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        angles = np.loadtxt(anglesFile, delimiter=None)
#        Y_test_NN_rot = rotate(Y_test_NN, angles)
##        Y_test_NN_rot = Y_test_NN
#        Y_test_NN_total = extrapolate(Y_test_NN_rot, totalNum, nb_output, 1)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        np.savetxt(path + filename, Y_test_NN_total, fmt='%.6f', delimiter=' ')
#################################
        
### predict test data
#yarnNum = 1
#skipFactor = 100        
#firstFrame = 200
#lastFrame = 200
#totalNum = 300 ################# NOTE: downsampled
#dataset = 'spacing1.0x_00011_woven'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "trainX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        X_test_ = X_test.reshape(X_test.shape[0], window_size, nb_features)
##        Y_test_NN = model.predict(X_test_) 
#        Y_test_NN = np.loadtxt(path + "trainY_" + str(f) + '_' + str(y) + ".txt", delimiter=None)
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        angles = np.loadtxt(anglesFile, delimiter=None)
#        Y_test_NN_rot = rotate(Y_test_NN, angles)
##        Y_test_NN_rot = Y_test_NN
#        Y_test_NN_total = extrapolate(Y_test_NN_rot, totalNum, nb_output, 1)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        np.savetxt(path + filename, Y_test_NN_total, fmt='%.6f', delimiter=' ')
        
# In[]:
#dataset = 'spacing0.5x_00011'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(16000/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f-firstFrame) + '_' + str(y) + ".txt",delimiter=None)
#        X_test_ = X_test.reshape(X_test.shape[0], window_size, nb_features)
#        Y_test_NN = model.predict(X_test_) 
#        Y_test_NN_total = extrapolate(Y_test_NN, totalNum, nb_output, 1)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        np.savetxt(path + filename, Y_test_NN_total, fmt='%.6f', delimiter=' ')
#
#dataset = 'spacing0.5x_10100'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(15000/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f-firstFrame) + '_' + str(y) + ".txt",delimiter=None)
#        X_test_ = X_test.reshape(X_test.shape[0], window_size, nb_features)
#        Y_test_NN = model.predict(X_test_) 
#        Y_test_NN_total = extrapolate(Y_test_NN, totalNum, nb_output, 1)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        np.savetxt(path + filename, Y_test_NN_total, fmt='%.6f', delimiter=' ')
#
#dataset = 'spacing0.5x_11110'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(15000/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f-firstFrame) + '_' + str(y) + ".txt",delimiter=None)
#        X_test_ = X_test.reshape(X_test.shape[0], window_size, nb_features)
#        Y_test_NN = model.predict(X_test_) 
#        Y_test_NN_total = extrapolate(Y_test_NN, totalNum, nb_output, 1)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        np.savetxt(path + filename, Y_test_NN_total, fmt='%.6f', delimiter=' ')
#        
#################
#dataset = 'spacing1.0x_00011'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(17000/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f-firstFrame) + '_' + str(y) + ".txt",delimiter=None)
#        X_test_ = X_test.reshape(X_test.shape[0], window_size, nb_features)
#        Y_test_NN = model.predict(X_test_) 
#        Y_test_NN_total = extrapolate(Y_test_NN, totalNum, nb_output, 1)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        np.savetxt(path + filename, Y_test_NN_total, fmt='%.6f', delimiter=' ')
#
#
#dataset = 'spacing1.0x_10100'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(15500/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f-firstFrame) + '_' + str(y) + ".txt",delimiter=None)
#        X_test_ = X_test.reshape(X_test.shape[0], window_size, nb_features)
#        Y_test_NN = model.predict(X_test_) 
#        Y_test_NN_total = extrapolate(Y_test_NN, totalNum, nb_output, 1)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        np.savetxt(path + filename, Y_test_NN_total, fmt='%.6f', delimiter=' ')
#
#dataset = 'spacing1.0x_11110'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(16000/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f-firstFrame) + '_' + str(y) + ".txt",delimiter=None)
#        X_test_ = X_test.reshape(X_test.shape[0], window_size, nb_features)
#        Y_test_NN = model.predict(X_test_) 
#        Y_test_NN_total = extrapolate(Y_test_NN, totalNum, nb_output, 1)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        np.savetxt(path + filename, Y_test_NN_total, fmt='%.6f', delimiter=' ')
#
##################
#dataset = 'spacing1.5x_00011'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(17500/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f-firstFrame) + '_' + str(y) + ".txt",delimiter=None)
#        X_test_ = X_test.reshape(X_test.shape[0], window_size, nb_features)
#        Y_test_NN = model.predict(X_test_) 
#        Y_test_NN_total = extrapolate(Y_test_NN, totalNum, nb_output, 1)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        np.savetxt(path + filename, Y_test_NN_total, fmt='%.6f', delimiter=' ')
#
#
#dataset = 'spacing1.5x_10100'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(16000/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f-firstFrame) + '_' + str(y) + ".txt",delimiter=None)
#        X_test_ = X_test.reshape(X_test.shape[0], window_size, nb_features)
#        Y_test_NN = model.predict(X_test_) 
#        Y_test_NN_total = extrapolate(Y_test_NN, totalNum, nb_output, 1)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        np.savetxt(path + filename, Y_test_NN_total, fmt='%.6f', delimiter=' ')
#
#dataset = 'spacing1.5x_11110'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(16500/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f-firstFrame) + '_' + str(y) + ".txt",delimiter=None)
#        X_test_ = X_test.reshape(X_test.shape[0], window_size, nb_features)
#        Y_test_NN = model.predict(X_test_) 
#        Y_test_NN_total = extrapolate(Y_test_NN, totalNum, nb_output, 1)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        np.savetxt(path + filename, Y_test_NN_total, fmt='%.6f', delimiter=' ')
#######################
#dataset = 'spacing0.5x'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(14000/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f-firstFrame) + '_' + str(y) + ".txt",delimiter=None)
#        X_test_ = X_test.reshape(X_test.shape[0], window_size, nb_features)
#        Y_test_NN = model.predict(X_test_) 
#        Y_test_NN_total = extrapolate(Y_test_NN, totalNum, nb_output, 1)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        np.savetxt(path + filename, Y_test_NN_total, fmt='%.6f', delimiter=' ')
#
#
#dataset = 'spacing1.0x'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(14500/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f-firstFrame) + '_' + str(y) + ".txt",delimiter=None)
#        X_test_ = X_test.reshape(X_test.shape[0], window_size, nb_features)
#        Y_test_NN = model.predict(X_test_) 
#        Y_test_NN_total = extrapolate(Y_test_NN, totalNum, nb_output, 1)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        np.savetxt(path + filename, Y_test_NN_total, fmt='%.6f', delimiter=' ')
#
#dataset = 'spacing1.5x'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(15000/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f-firstFrame) + '_' + str(y) + ".txt",delimiter=None)
#        X_test_ = X_test.reshape(X_test.shape[0], window_size, nb_features)
#        Y_test_NN = model.predict(X_test_) 
#        Y_test_NN_total = extrapolate(Y_test_NN, totalNum, nb_output, 1)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        np.savetxt(path + filename, Y_test_NN_total, fmt='%.6f', delimiter=' ')