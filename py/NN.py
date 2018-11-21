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

from keras import backend as K

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
    nb_outputs = Y_train_all.shape[1]
    
    return (nb_features,nb_outputs, scaler)

    
def loadData(fn_trainX, fn_trainY, fn_validX, fn_validY):
    X_train_all = np.loadtxt(fn_trainX,delimiter=None)
    Y_train_all = np.loadtxt(fn_trainY,delimiter=None)
    
    #to take only part of data
#    X_train_all = X_train_all2[0:20, :]
#    Y_train_all = Y_train_all2[0:20, :]

    print("Original training data shape (X): ", X_train_all.shape)
    print("Original training data shape (Y): ", Y_train_all.shape)
    
    #rescale the output data
    scaler = MinMaxScaler(feature_range=(0, 1))
    scaler.fit(Y_train_all)
    Y_train_all = scaler.transform(Y_train_all)
    
    nb_features = X_train_all.shape[1]
    nb_traindata = X_train_all.shape[0]
    split = 0.99
    nb_halfdata = int ( round(nb_traindata*split) / 16 ) * 16 #make the number of data to be divisible by 16 so all batches has exastly same number of data
    print (round(nb_traindata*split), nb_halfdata)
    nb_outputs = Y_train_all.shape[1]
    
    # using subset data as training and validation
    all_train = np.concatenate((X_train_all,Y_train_all), axis=1)  
    np.random.shuffle(all_train)
    X_train = all_train[0:nb_halfdata,0:nb_features]
    Y_train = all_train[0:nb_halfdata,nb_features:]  
    
    ### prepare the validation data
    X_valid_all = all_train[nb_halfdata:,0:nb_features]
    Y_valid_all = all_train[nb_halfdata:,nb_features:] 
    
#    X_valid_all = np.loadtxt(fn_validX,delimiter=None)
#    Y_valid_all = np.loadtxt(fn_validY,delimiter=None)
#    
    all_valid = np.concatenate((X_valid_all,Y_valid_all), axis=1)
    np.random.shuffle(all_valid)
    split_val = 0.2
    nb_validdata = int (  round(X_valid_all.shape[0] * split_val) / 16 ) * 16 #make the number of data to be divisible by 16
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


## Build a customized loss function
# In[]:
def customLoss(nbatch, pnts, nPnt, lenPnt):
    def loss(y_true, y_pred):
        w = 0
        err = 0 
        for i in range (0, nbatch):
            y_true_mat = K.reshape(y_true[i], [-1,2])
            y_pred_mat = K.reshape(y_pred[i], [-1,2]) #-1 means: reshape so that columns will be 2 with whatever number of row
            for j in range (0, nPnt):
                pnt_v = K.reshape(pnts[j], [-1,1])
                pnt_pred = K.dot(y_pred_mat, pnt_v)
                pnt_true = K.dot(y_true_mat, pnt_v)
                err = err + K.mean( K.square(pnt_pred - pnt_true ) )
 
        err /= ( nbatch * nPnt )
        c = 0.001
        w = c/lenPnt
        print('w for regularization term is ', w )
        return  K.mean(K.square(y_true-y_pred) ) + w*err
    return loss

## Build neural network model
# In[]:
def buildModel(input_dim, output_dim, neurons, w_path, nbatch, pnts, nPnts, lenPnts):

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
    
    model.add(Dense(output_dim))
    model.add(Activation('linear'))
#    model.add(Activation('sigmoid'))

    
#    model.compile(loss='mse', optimizer='adam', metrics=['mse'])
    model.compile(loss=customLoss( nbatch, pnts, nPnts, lenPnts), optimizer='adam', metrics=['mse'])
    
    return model

## Train the model
# In[]
def trainModel(model, X_train, Y_train, X_valid, Y_valid, nbatch):
    
    # Weights are updated one mini-batch at a time. A running average of the training loss is computed in real time, which is useful for identifying problems (e.g. the loss might explode or get stuck right). The validation loss is evaluated at the end of each epoch (without dropout).
    history = model.fit(X_train, Y_train, batch_size = nbatch, epochs = 2, verbose = 2,
                        validation_data=(X_valid, Y_valid))

    # Plot loss trajectory throughout training.
    plt.figure()
    plt.subplot(1,2,1)
    plt.plot(history.history['mean_squared_error'], label='train')
    plt.plot(history.history['val_mean_squared_error'], label='valid')
    plt.xlabel('Epoch')
    plt.ylabel('mse')
    plt.legend()
    plt.show()
    
    # ## Evaluate performance    
    # Note: when calling evalate, dropout is automatically turned off.
    score = model.evaluate(X_valid, Y_valid, verbose=0)
    print('Validation loss: %0.5f' % score[0])

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
            predictExtr[i] = predicted[0]
        elif k == int(vrtNum/stride):
            predictExtr[i] = predicted[int(vrtNum/stride) - 1 ]
        else: 
            w = (fidx[k] - i)/(fidx[k] - fidx[k - 1])
            predictExtr[i] = w*predicted[k - 1] + (1.0 - w)*predicted[k]
    
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
    

    predicted = model.predict(X_test, verbose=0)
    predicted = scaler.inverse_transform(predicted)

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
#    uncomment sweep only for stretching
    predicted_sweep = predicted
#    predicted_sweep = sweep(predicted, window_sweep)
    predicted_reg = regularize(predicted_sweep, window_reg)
    predicted_reg1 = regularize(predicted_reg, window_reg+2)
    predicted_reg2 = regularize(predicted_reg1, window_reg) 
#    predicted_reg2 = predicted_sweep
    
    predicted_us = upsample(predicted_reg2, upsample_rate)
    predicted_total = extrapolate(predicted_us, vrtxNum, nb_outputs, stride)
#    np.savetxt(path + filename, predicted_total, fmt='%.6f', delimiter=' ')
    
    return predicted_total
## Main
# In[]    
def test(neurons,fn_trainX, fn_trainY, fn_validX, fn_validY, reTrain, w_path): 
    storeModel = w_path + 'train_all/model_ws5_new.h5'
    nbatch = 16
    # Read point cloud needed for loss-function
    pnts = np.loadtxt(w_path+'train_all/ref2D.txt')
    lenPnts = (np.linalg.norm(pnts))
    nPnts = pnts.shape[0]#number of points in one cross-section
    pnt_K = K.variable(value=pnts)
    print('point cloud shape: ', pnts.shape)
    
    if (reTrain):
        print ('<<<<< train the model >>>>>>')        
        
        (X_train, Y_train, X_valid, Y_valid, nb_features, nb_outputs, scaler) = loadData(fn_trainX, fn_trainY, fn_validX, fn_validY)
        model = buildModel(nb_features, nb_outputs, neurons, w_path, nbatch, pnt_K, nPnts, lenPnts)
        model = trainModel(model, X_train, Y_train, X_valid, Y_valid, nbatch)
        model.save(storeModel)
    else:
        print ('<<<<< load the stored model >>>>>')
        (nb_features, nb_outputs, scaler) = loadScaler(fn_trainX, fn_trainY)
#        del model  # deletes the existing model
        model = load_model(storeModel, custom_objects={'loss': customLoss( nbatch, pnt_K, nPnts, lenPnts) }) 
#        model = load_model(storeModel)
        
    return model, scaler, nb_outputs

# In[]:
def append2sets(dataset2, w_path):
    path1 = w_path + 'train_all/'
    path2 = w_path + dataset2 
    X_train_all_1 = np.loadtxt(path1 + "trainX_all.txt",delimiter=None)
    Y_train_all_1 = np.loadtxt(path1 + "trainY_all.txt",delimiter=None)
    X_train_all_2 = np.loadtxt(path2 + "trainX_all.txt",delimiter=None)
    Y_train_all_2 = np.loadtxt(path2 + "trainY_all.txt",delimiter=None)
    
    X_train_all = np.concatenate((X_train_all_1,X_train_all_2), axis=0)
    Y_train_all = np.concatenate((Y_train_all_1,Y_train_all_2), axis=0)
    
    np.savetxt(w_path + 'train_all/trainX_all.txt', X_train_all, fmt='%.6f', delimiter=' ')
    np.savetxt(w_path + 'train_all/trainY_all.txt', Y_train_all, fmt='%.6f', delimiter=' ')    
    
## append training data from different datasets
# In[] 
def appendTrainingData(datasets, w_path, fn_trainX, fn_trainY):
    for i in range (0, len(datasets)):
        if (i==0):
            p0x = w_path + datasets[i] + '/trainX_all.txt'
            p0y = w_path + datasets[i] + '/trainY_all.txt'
            copyfile (p0x, fn_trainX)
            copyfile (p0y, fn_trainY)
        else:
            append2sets(datasets[i],w_path)
            print('appended ' + datasets[i])
        
# In[]
stride = 1
#window_reg = 5 #for stretching
window_reg = 5 #for non-stretching
window_sweep = 5
def main_NN(yarn_type,upsample_rate, dataset, firstFrame, lastFrame, yarn0, yarn1, skipFactor, vrtxNum):
    datasets = []
    
    config = yarn_type + '/train/'
    loc = 'F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'
    w_path =  loc + config
    datasets.append('spacing0.5x/10/')
    datasets.append('spacing0.5x/00011/')
    datasets.append('spacing0.5x/10100/')
    datasets.append('spacing0.5x/11110/')
    datasets.append('spacing1.0x/10/')
    datasets.append('spacing1.0x/00011/')
    datasets.append('spacing1.0x/10100/')
    datasets.append('spacing1.0x/11110/')
    datasets.append('spacing1.5x/10/')
    datasets.append('spacing1.5x/00011/')
    datasets.append('spacing1.5x/10100/')
    datasets.append('spacing1.5x/11110/')
    datasets.append('spacing2.0x/10/')
    datasets.append('spacing2.0x/00011/')
    datasets.append('spacing2.0x/10100/')
    datasets.append('spacing2.0x/11110/')
    datasets.append('spacing2.5x/10/')
    datasets.append('spacing2.5x/00011/')
    datasets.append('spacing2.5x/10100/')
    datasets.append('spacing2.5x/11110/')
    datasets.append('spacing3.0x/10/')
    datasets.append('spacing3.0x/00011/')
    datasets.append('spacing3.0x/10100/')
    datasets.append('spacing3.0x/11110/')
    
    fn_trainX = w_path + "train_all/trainX_all.txt"
    fn_trainY = w_path + "train_all/trainY_all.txt"
    
    fn_validX = loc + 'single_yarn/' + yarn_type + '/stretch/trainX_all.txt'
    fn_validY = loc + 'single_yarn/' + yarn_type + '/stretch/trainY_all.txt'
    
    reTrain = 1
#    if (reTrain):
#        appendTrainingData(datasets, w_path, fn_trainX, fn_trainY)
       
    model, scaler, nb_outputs = test(256, fn_trainX, fn_trainY, fn_validX, fn_validY, reTrain, w_path)
    
    print ('***** predict for the test data ******')
    path = 'F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset + '/'
    frame0 = int(firstFrame/skipFactor)
    frame1 = int(lastFrame/skipFactor) + 1
    for i in range (frame0, frame1):
        f = i*skipFactor
        print("predicted frame ", f)
        for y in range (yarn0, yarn1):
            X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
            #X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + "_temporal.txt",delimiter=None)
            filename = "testY_NN_" + str(f) + '_' + str(y) +  ".txt"
            anglesFile = path + "angle_" + str(f) + '_' + str(y) + ".txt"
            isRot = 1 
            predicted_total = predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, isRot, upsample_rate, path)
            np.savetxt(path + filename, predicted_total, fmt='%.6f', delimiter=' ') 
  
