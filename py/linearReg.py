
import numpy as np
import matplotlib.pyplot as plt
#import mltools as ml
from keras.models import Sequential
from keras.layers import Dense, Activation
from sklearn.preprocessing import MinMaxScaler
from keras.layers.normalization import BatchNormalization
import math
from shutil import copyfile
from keras.models import load_model
from sklearn import linear_model
from sklearn.preprocessing import PolynomialFeatures

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

## Load data
# In[]
def loadScaler(fn_trainX, fn_trainY):
    X_train_all = np.loadtxt(fn_trainX,delimiter=None)
    Y_train_all = np.loadtxt(fn_trainY,delimiter=None)

    print("Original training data shape (X): ", X_train_all.shape)
    print("Original training data shape (Y): ", Y_train_all.shape)
    
    #rescale the output data
    scaler = MinMaxScaler(feature_range=(0, 1))
    scaler.fit(Y_train_all)
    Y_train_all = scaler.transform(Y_train_all)
    
    nb_features = X_train_all.shape[1]
    nb_traindata = X_train_all.shape[0]
    split = 1.0
    nb_halfdata = round(nb_traindata*split)
    nb_outputs = Y_train_all.shape[1]
    
    return (nb_features,nb_outputs, scaler)

    
def loadData(fn_trainX, fn_trainY, fn_validX, fn_validY):
    X_train_all = np.loadtxt(fn_trainX,delimiter=None)
    Y_train_all = np.loadtxt(fn_trainY,delimiter=None)

    print("Original training data shape (X): ", X_train_all.shape)
    print("Original training data shape (Y): ", Y_train_all.shape)
    
    #rescale the output data
    scaler = MinMaxScaler(feature_range=(0, 1))
    scaler.fit(Y_train_all)
    Y_train_all = scaler.transform(Y_train_all)
    
    nb_features = X_train_all.shape[1]
    nb_traindata = X_train_all.shape[0]
    split = 0.99
    nb_halfdata = round(nb_traindata*split)
    nb_outputs = Y_train_all.shape[1]
    
    # using subset data as training and validation
    all_train = np.concatenate((X_train_all,Y_train_all), axis=1)  
    np.random.shuffle(all_train)
    X_train = all_train[0:nb_halfdata,0:nb_features]
    Y_train = all_train[0:nb_halfdata,nb_features:]  
    
    ### prepare the validation data
#    X_valid = all_train[nb_halfdata:,0:nb_features]
#    Y_valid = all_train[nb_halfdata:,nb_features:] 
    
    X_valid_all = np.loadtxt(fn_validX,delimiter=None)
    Y_valid_all = np.loadtxt(fn_validY,delimiter=None)
    all_valid = np.concatenate((X_valid_all,Y_valid_all), axis=1)
    np.random.shuffle(all_valid)
    split_val = 0.2
    nb_validdata = round(X_valid_all.shape[0] * split_val)
    X_valid = all_valid[0:nb_validdata,0:nb_features]
    Y_valid = all_valid[0:nb_validdata,nb_features:]
    
    # polynomio
#    X_train_,params = ml.transforms.rescale(X_train);
#    X_valid_,_ = ml.transforms.rescale( X_valid, params);
#    X_test_,_ = ml.transforms.rescale( X_test, params);
#    nb_features = X_train_.shape[1]
    
    # Represent the targets as one-hot vectors: e.g. 0 -> [1,0];  1 -> [0, 1].
    print("Training Y matrix shape: ", Y_train.shape)
    print("Validation Y matrix shape: ", Y_valid.shape)

     
#    return (X_train_, Y_train, X_valid_, Y_valid, nb_features,nb_outputs, X_test_)
    return (X_train, Y_train, X_valid, Y_valid, nb_features,nb_outputs, scaler)

## Build neural network model
# In[]:
def buildModel(input_dim, output_dim, neurons):

    # Simple fully-connected neural network with 2 hidden layers.
    # Including dropout layer helps avoid overfitting.
    model = Sequential()      
    
    model.add(Dense(neurons, input_dim=input_dim))
    model.add(Activation('relu'))   
    model.add(BatchNormalization())

    model.add(Dense(neurons))
    model.add(Activation('relu')) 
    model.add(BatchNormalization())
    
    model.add(Dense(neurons))
    model.add(Activation('relu')) 
    model.add(BatchNormalization())
    
#    model.add(Dense(neurons))
#    model.add(Activation('relu')) 
#    model.add(BatchNormalization())
    
    model.add(Dense(output_dim))
    model.add(Activation('linear'))
#    model.add(Activation('sigmoid'))
    
    model.compile(optimizer='adam', loss='mse', metrics=['mse'])
    
    return model

## Train the model
# In[]
def trainModel(model, X_train, Y_train, X_valid, Y_valid):
    
    # Weights are updated one mini-batch at a time. A running average of the training loss is computed in real time, which is useful for identifying problems (e.g. the loss might explode or get stuck right). The validation loss is evaluated at the end of each epoch (without dropout).
#    history = model.fit(X_train, Y_train, batch_size = 16, epochs = 30, verbose = 2,
#                        validation_data=(X_valid, Y_valid))
    ################ linear regression
    print ("fitting the model...")
#    model = linear_model.LinearRegression()
# Train the model using the training sets
#    history = model.fit(X_train_, Y_train_)   

    ############### polynomial regression degree 2
    # create a Linear Regressor 
    # Polynomial regression is a special case of linear regression
#    model = linear_model.LinearRegression()
#    
#    # pass the order of your polynomial here  
#    poly = PolynomialFeatures(degree=10)
#    
#    # convert to be used further to linear regression
#    X_train_ = poly.fit_transform(X_train)
#
#    # Train the model using the training sets
#    history = model.fit(X_train_, Y_train)
    ###############
    
    # Instantiate a Gaussian Process model
    kernel = C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e2))
    #RBF( length_scale, length_scale_bounds )
    # RBF is covarrience function
    model = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=9)
    # Train the model using the training sets
    history = model.fit(X_train, Y_train)

    # Plot loss trajectory throughout training.
#    plt.figure()
#    plt.subplot(1,2,1)
#    plt.plot(history.history['mean_squared_error'], label='train')
#    plt.plot(history.history['val_mean_squared_error'], label='valid')
#    plt.xlabel('Epoch')
#    plt.ylabel('mse')
#    plt.legend()
#    plt.show()
    
    # ## Evaluate performance    
    # Note: when calling evaluate, dropout is automatically turned off.
#    score = model.evaluate(X_valid, Y_valid, verbose=0)
#    print('Validation loss: %0.5f' % score[0])

    return model

## extrapolate
# In[]:
def extrapolate(predicted, totalNum, nb_outputs, stride):
    total = np.zeros(totalNum*nb_outputs).reshape(totalNum,nb_outputs)
    assert(predicted.shape[1] == total.shape[1])
    # extrapolate predicted by stride_num
    vrtNum = predicted.shape[0] * stride
    predictExtr = np.zeros(vrtNum*nb_outputs).reshape(vrtNum,nb_outputs)
    fidx = np.linspace (1,vrtNum-stride-1, num=int(vrtNum/stride), dtype=int)

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
    #fill up the two ends with the first and last value
    for i in range (0,r):
        total[i, : ] = predictExtr[0]
    for j in range (r+predictExtr.shape[0], totalNum):
        total[j, : ] = predictExtr[-1]
    return total
    

## rotate the NN output back to RM frames
# In[]:
def rotate(predicted, angles):
    predicted_rot = predicted
    n = len(angles)
    for i in range (0,n):
        c, s = np.cos(-1*angles[i]), np.sin(-1*angles[i])
        R = np.matrix('{} {}; {} {}'.format(c, -s, s, c))
        S = predicted[i].reshape([2,2])
        rot = R*S*R.transpose()
        predicted_rot[i] = np.array([rot[0,0] , rot[0,1] , rot[1,0] , rot[1,1] ])
    return predicted_rot

## upsample the output
# In[]:
def upsample(predicted, upsample):
    predicted_us = []
    sz = predicted.shape[0]
    for i in range (0, sz-1):
        for u in range (0, upsample):
            w = u * 1.0/upsample            
            tmp = (1-w) * predicted[i] + w*predicted[i+1]
            predicted_us.append(tmp)
    predicted_us.append(predicted[-1])
    predicted_us_np = np.array(predicted_us)
    return predicted_us_np

## sweep NN-output for potential noise using median filter
# In[]:
def sweep(predicted, window_sweep):
    sz = predicted.shape[0]
    predicted_reg = []
    for i in range (0, sz - window_sweep + 1):
        window_arr_tmp = []
        for w in range (0, window_sweep):
            window_arr_tmp.append(predicted[i+w][0])
        window_arr = np.array(window_arr_tmp)
#        window_arr = np.array([predicted[i][0], predicted[i+1][0], predicted[i+2][0] ]) 
        med = np.median(window_arr)
        idx = (np.abs(window_arr - med)).argmin()
        predicted_reg.append(predicted[idx])

    predicted_reg_np = np.array(predicted_reg)
    return predicted_reg_np        
        
    
## regularize NN-output
# In[]:
def regularize(predicted, window_reg):
    sz = predicted.shape[0]
    predicted_reg = []
    for i in range (0, sz - window_reg + 1):
        sum_window = np.zeros(predicted.shape[1])
        sum_w = 0       
        
#        s = 0
#        if (window_reg == 7):
#            gauss = [0.01, 0.05, 0.32, 1, 0.32, 0.05, 0.01]
#        if (window_reg == 5):
#            gauss = [0.05, 0.32, 1, 0.32, 0.05]
#        if (window_reg ==3):
#            gauss = [ 0.32, 1, 0.32 ]
#        for w in gauss:
#            sum_window = w*predicted[i+s] + sum_window
#            sum_w = sum_w + w
#            s = s+1
        
        
        for s in range (0, window_reg):         
            cntr = int (window_reg/2)  #index of middle point
            w = abs( abs(s - cntr) - (cntr +1) )        
#            w=1  #for mean filter
            sum_window = w*predicted[i+s] + sum_window
            sum_w = sum_w + w
            
            
        avg_window = sum_window/sum_w
        predicted_reg.append(avg_window)
    
    predicted_reg_np = np.array(predicted_reg)
      
    return predicted_reg_np
#    return predicted_reg_np_again


## Prediction
# In[]:
def predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, isRot, upsample_rate, path):
    
    predicted = model.predict(X_test)
#    all_test = np.concatenate((X_test,predicted), axis=1)
    predicted = scaler.inverse_transform(predicted)
#    np.savetxt(path + 'testY_NN.txt', all_test[:, -3:], fmt='%.6f', delimiter=' ')
    
    # rotate the shape-match back to original frame (window-Rotation-minimizing to yarn-rotation-minimizing)
#    predicted = np.loadtxt(path + 'trainY_15000_0.txt')
    if isRot:
        angles = np.loadtxt(anglesFile, delimiter=None)
        predicted = rotate(predicted, angles)

    ########## Regularize synthesized example ##########
#    predicted_sweep = sweep(predicted, window_sweep)
#    predicted_sweep2 = sweep(predicted_sweep, window_sweep)
#    predicted_reg = regularize(predicted_sweep2, window_reg) 
#    predicted_reg2 = regularize(predicted_reg, window_reg+2) 
    ####################
#    predicted_sweep = predicted
#    predicted_sweep = sweep(predicted, window_sweep)
#    predicted_reg = regularize(predicted_sweep, window_reg) 
#    predicted_reg2 = regularize(predicted_reg, window_reg+2) 
#    predicted_reg2 = predicted_sweep
    
    predicted_us = upsample(predicted, upsample_rate)
    np.savetxt(path + 'testY_NN.txt', predicted_us, fmt='%.6f', delimiter=' ')
    
    predicted_total = extrapolate(predicted_us, vrtxNum, nb_outputs, stride)
#    np.savetxt(path + filename, predicted_total, fmt='%.6f', delimiter=' ')
    return predicted_total
## Main
# In[]    
def test(neurons,fn_trainX, fn_trainY, fn_validX, fn_validY, reTrain, w_path): 
    storeModel = w_path + 'train_all/model_ws5.h5'
    if (reTrain):
        (X_train, Y_train, X_valid, Y_valid, nb_features, nb_outputs, scaler) = loadData(fn_trainX, fn_trainY, fn_validX, fn_validY)
        model = buildModel(nb_features, nb_outputs, neurons)
        model = trainModel(model, X_train, Y_train, X_valid, Y_valid)
#        model.save(storeModel)  # creates a HDF5 file 'my_model.h5'
    else:
        (nb_features, nb_outputs, scaler) = loadScaler(fn_trainX, fn_trainY)
#        del model  # deletes the existing model
#        model = load_model(storeModel) 
        
    return model, scaler, nb_outputs

# In[]:
def append2sets(dataset2, w_path):
    path1 = w_path + 'train_all/'
    path2 = w_path + dataset2 + '/NN/'
    X_train_all_1 = np.loadtxt(path1 + "trainX_all.txt",delimiter=None)
    Y_train_all_1 = np.loadtxt(path1 + "trainY_all.txt",delimiter=None)
    X_train_all_2 = np.loadtxt(path2 + "trainX_all.txt",delimiter=None)
    Y_train_all_2 = np.loadtxt(path2 + "trainY_all.txt",delimiter=None)
    
    X_train_all = np.concatenate((X_train_all_1,X_train_all_2), axis=0)
    Y_train_all = np.concatenate((Y_train_all_1,Y_train_all_2), axis=0)
    
    np.savetxt(w_path + 'train_all/trainX_all.txt', X_train_all, fmt='%.6f', delimiter=' ')
    np.savetxt(w_path + 'train_all/trainY_all.txt', Y_train_all, fmt='%.6f', delimiter=' ')  
    print('all training are written to ', w_path + 'train_all' )
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
            append2sets(datasets[i],w_path)
        
# In[]
stride = 1
window_reg = 3
window_sweep = 3
def main_NN(yarn_type,upsample_rate, dataset, firstFrame, lastFrame, yarn0, yarn1, skipFactor, vrtxNum):
    datasets = []
    
    config = 'pattern/'+ yarn_type +'/'
    loc = 'F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'
    w_path =  loc + config
#    datasets.append('spacing0.5x/10')
#    datasets.append('spacing0.5x/00011')
#    datasets.append('spacing0.5x/10100')
#    datasets.append('spacing0.5x/11110')
    datasets.append('spacing1.0x/10')
    datasets.append('spacing1.0x/00011')
    datasets.append('spacing1.0x/10100')
    datasets.append('spacing1.0x/11110')
    datasets.append('spacing1.5x/10')
    datasets.append('spacing1.5x/00011')
    datasets.append('spacing1.5x/10100')
    datasets.append('spacing1.5x/11110')
    datasets.append('spacing2.0x/10')
    datasets.append('spacing2.0x/00011')
    datasets.append('spacing2.0x/10100')
    datasets.append('spacing2.0x/11110')
    datasets.append('spacing2.5x/10')
    datasets.append('spacing2.5x/00011')
    datasets.append('spacing2.5x/10100')
    datasets.append('spacing2.5x/11110')
    datasets.append('spacing3.0x/10')
    datasets.append('spacing3.0x/00011')
    datasets.append('spacing3.0x/10100')
    datasets.append('spacing3.0x/11110')
    
    datasets.append('../../yarn_stretch')
    
    fn_trainX = w_path + "train_all/trainX_all.txt"
    fn_trainY = w_path + "train_all/trainY_all.txt"
    
    fn_validX = loc + 'single_yarn/' + yarn_type + '/stretch/trainX_all.txt'
    fn_validY = loc + 'single_yarn/' + yarn_type + '/stretch/trainY_all.txt'
    
    reTrain = 1
    if (reTrain):
        appendTrainingData(datasets, w_path, fn_trainX, fn_trainY)
        
    model, scaler, nb_outputs = test(256, fn_trainX, fn_trainY, fn_validX, fn_validY, reTrain, w_path)
    #upsample_rate = 1
    
    
    
    #yarn0 = 150
    #yarn1 = 250
    #skipFactor = 200       
    #vrtxNum = 397  #300*2 ## after upsampling
    #firstFrame = 400
    #lastFrame = 6500
    
    
    #dataset = 'stretch/yarn4/stretch'
    #dataset = 'fall/yarn4/fall'
    #dataset = 'ball_fall'
    #dataset = 'twist/yarn4/damp2_500'
    #dataset = 'woven/yarn4/spacing1.0x/00011/shear'
    #dataset = 'pattern/yarn4/spacing0.5x/10/Raymond'
    #dataset = 'single_yarn/yarn11/stretch'
    #dataset = 'single_yarn/yarn4/teeth/4_1.6'
    #dataset = 'woven/6x6' 
    #dataset = 'woven/arbitrary_pattern/150x100'
    #dataset = 'woven/stretch/yarn4/100x100'
    #dataset = 'woven/push/yarn8/100x100'
    
    path = 'F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/'
    frame0 = int(firstFrame/skipFactor)
    frame1 = int(lastFrame/skipFactor) + 1
    for i in range (frame0, frame1):
        f = i*skipFactor
        print(i)
        for y in range (yarn0, yarn1):
            X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
            
            #####
#            # pass the order of your polynomial here  
#            poly = PolynomialFeatures(degree=10)
#            # convert to be used further to linear regression
#            X_test_ = poly.fit_transform(X_test)
#            X_test = X_test_
            #####
    
            filename = "testY_NN_" + str(f) + '_' + str(y) +  ".txt"
            anglesFile = path + "angle_" + str(f) + '_' + str(y) + ".txt"
            isRot = 1
            predicted_total = predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, isRot, upsample_rate, path)
            np.savetxt(path + filename, predicted_total, fmt='%.6f', delimiter=' ') 
  
    
# In[] 
#skipFactor = 1000
#downSample = 2
#vrtNum = 150 ###before upsampling
#fiberNum = 160
#yarn0 = 0
#yarn1 = 1
#isTrain = 1
##dataset = 'woven/arbitrary_pattern/150x100'
#dataset = 'single_yarn/yarn8/stretch'
#firstFrame = 25000
#lastFrame = 35000
#
########################### NN
#vrtx_us = vrtNum*downSample
#main_NN('yarn8',downSample, dataset, firstFrame, lastFrame, yarn0, yarn1, skipFactor, vrtx_us )
############################
# In[] 

#yarn0 = 0
#yarn1 = 1
#
#stride = 1
#skipFactor = 2000       
#vrtxNum = 34607 #3291 #300*2 ## after upsampling
#firstFrame = 2000
#lastFrame = 36000
#
##dataset = 'stretch/yarn4/stretch'
##dataset = 'fall/yarn4/fall'
##dataset = 'ball_fall'
##dataset = 'twist/yarn4/damp2_500'
##dataset = 'woven/yarn4/spacing1.0x/00011/shear'
##dataset = 'pattern/yarn4/spacing0.5x/10/Raymond'
##dataset = 'single_yarn/yarn11/stretch'
##dataset = 'single_yarn/yarn4/teeth/4_1.6'
##dataset = 'woven/6x6' 
#dataset = 'woven/knitted'
##dataset = 'woven/arbitrary_pattern/512x512'
##dataset = 'woven/arbitrary_pattern/100x100'
##dataset = 'woven/stretch/yarn4/100x100'
##dataset = 'woven/push/yarn8/100x100'
#
#path = 'F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    print(i)
#    for y in range (yarn0, yarn1):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predicted_total = predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        np.savetxt(path + filename, predicted_total, fmt='%.6f', delimiter=' ') 
##  
# In[] 
#yarn0 = 0
#yarn1 = 12
#
#stride = 1
#skipFactor = 500       
#vrtxNum = 102 #300*2 ## after upsampling
#firstFrame = 15000
#lastFrame = 15000
#
#
##dataset = 'stretch/yarn4/stretch'
##dataset = 'fall/yarn4/fall'
##dataset = 'ball_fall'
##dataset = 'twist/yarn4/damp2_500'
##dataset = 'woven/yarn4/spacing1.0x/00011/shear'
##dataset = 'pattern/yarn4/spacing0.5x/10/Raymond'
##dataset = 'single_yarn/yarn11/stretch'
##dataset = 'single_yarn/yarn4/teeth/4_1.6'
#dataset = 'woven/6x6' 
##dataset = 'woven/arbitrary_pattern/100x100'
##dataset = 'woven/stretch/yarn4/100x100'
##dataset = 'woven/push/yarn4/100x100'
#
#path = 'F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    print(i)
#    for y in range (yarn0, yarn1):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predicted_total = predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        np.savetxt(path + filename, predicted_total, fmt='%.6f', delimiter=' ') 
#  
# In[] 
#yarn0 = 0
#yarn1 = 1
#
#stride = 1
#skipFactor = 1000        
#vrtxNum = 150 ## after upsampling
#firstFrame = 20000
#lastFrame = 20000
#
#
##dataset = 'fall/yarn4/fall'
##dataset = 'ball_fall'
##dataset = 'twist/yarn4/damp2_500'
##dataset = 'woven/yarn4/spacing1.0x/00011/shear'
##dataset = 'pattern/yarn100/spacing2.5x/10'
#dataset = 'single_yarn/yarn9/stretch'
##dataset = 'single_yarn/yarn11/teeth/4_1.6'
##dataset = 'single_yarn/yarn4/teeth/4_1.2'
##dataset = 'single_yarn/yarn8/teeth/4_1.2_00110'
##dataset = 'single_yarn/yarn100'
##dataset = 'woven/release/yarn9/8x8' 
##dataset = 'woven/release/yarn4/100x100'
##dataset = 'single_yarn/yarn4/teeth/1.2_110'
#
#path = 'F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    print(f)
#    for y in range (yarn0, yarn1):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predicted_total = predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        np.savetxt(path + filename, predicted_total, fmt='%.6f', delimiter=' ') 
#       

# In[]
## Regularize the global rotation for 6x6 result:
#fn_global_rot = '../input/woven/6x6/'
#for i in range (0,12):
#    global_rot = np.loadtxt(fn_global_rot + "global_rot_15000_" +str(i)+".txt",delimiter=None)
##    global_rot_sweep = global_rot
#    global_rot_sweep = sweep(global_rot, 3)
#    global_rot_reg = regularize(global_rot_sweep, 3)
#    
#    diff = global_rot.shape[0] - global_rot_reg.shape[0]
#    for j in range (0,diff):
#        add = np.array([ [1,0, 0,1] ])
#        global_rot_reg = np.concatenate((global_rot_reg,add), axis=0)
#    np.savetxt(fn_global_rot + "global_rot_15000_" +str(i)+".txt", global_rot_reg, fmt='%.6f', delimiter=' ')
#            