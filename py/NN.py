
import numpy as np
import matplotlib.pyplot as plt
#import mltools as ml
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout
from sklearn.preprocessing import MinMaxScaler
from keras.layers.normalization import BatchNormalization
import math

## Load data
# In[]
 
def loadData():
    
#    X_train_all = np.loadtxt(path + "trainX_all.txt",delimiter=None)
#    Y_train_all = np.loadtxt(path + "trainY_all.txt",delimiter=None)
    X_train_all = np.loadtxt("../input/all/NN/trainX_all.txt",delimiter=None)
    Y_train_all = np.loadtxt("../input/all/NN/trainY_all.txt",delimiter=None)

    print("Original training data shape (X): ", X_train_all.shape)
    print("Original training data shape (Y): ", Y_train_all.shape)
    
    #rescale the output data
    scaler = MinMaxScaler(feature_range=(0, 1))
    scaler.fit(Y_train_all)
    Y_train_all = scaler.transform(Y_train_all)
    
    nb_features = X_train_all.shape[1]
    nb_traindata = X_train_all.shape[0]
    split = 0.75
    nb_halfdata = round(nb_traindata*split)
    nb_outputs = Y_train_all.shape[1]
    
    # using subset data as training and validation
    all_train = np.concatenate((X_train_all,Y_train_all), axis=1) 
    
    
    np.random.shuffle(all_train)
    X_train = all_train[0:nb_halfdata,0:nb_features]
    Y_train = all_train[0:nb_halfdata,nb_features:]   
    X_valid = all_train[nb_halfdata:,0:nb_features]
    Y_valid = all_train[nb_halfdata:,nb_features:] 
     
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
    
#    model.add(Dense(126))
#    model.add(Activation('relu')) 
#    model.add(BatchNormalization())
    
    model.add(Dense(output_dim))
    model.add(Activation('linear'))
    
    model.compile(optimizer='adam', loss='mse', metrics=['mse'])
    
    return model

## Train the model
# In[]
def trainModel(model, X_train, Y_train, X_valid, Y_valid):
    
    # Weights are updated one mini-batch at a time. A running average of the training loss is computed in real time, which is useful for identifying problems (e.g. the loss might explode or get stuck right). The validation loss is evaluated at the end of each epoch (without dropout).
    history = model.fit(X_train, Y_train, batch_size = 16, epochs = 80, verbose = 2,
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
    # Note: when calling evaluate, dropout is automatically turned off.
    score = model.evaluate(X_valid, Y_valid, verbose=0)
    print('Validation loss: %0.5f' % score[0])

    return model

## extrapolate
# In[]:
def extrapolate(predicted, totalNum, filename, nb_outputs, stride):
    total = np.zeros(totalNum*nb_outputs).reshape(totalNum,nb_outputs)
    assert(predicted.shape[1] == total.shape[1])
    # extrapolate predicted by stride_num
    vrtNum = predicted.shape[0] * stride
    predictExtr = np.zeros(vrtNum*nb_outputs).reshape(vrtNum,nb_outputs)
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
    np.savetxt(path + filename, total, fmt='%.6f', delimiter=' ')
    
## Prediction
# In[]:
def predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride):
    
    predicted = model.predict(X_test, verbose=0)
#    all_test = np.concatenate((X_test,predicted), axis=1)
    predicted = scaler.inverse_transform(predicted)
#    np.savetxt(path + 'testY_NN.txt', all_test[:, -3:], fmt='%.6f', delimiter=' ')
    np.savetxt(path + 'testY_NN.txt', predicted, fmt='%.6f', delimiter=' ')
    
    extrapolate(predicted, vrtxNum, filename, nb_outputs, stride)
## Main
# In[]    
def test(neurons): 
    (X_train, Y_train, X_valid, Y_valid, nb_features, nb_outputs, scaler) = loadData()
    model = buildModel(nb_features, nb_outputs, neurons)
    model = trainModel(model, X_train, Y_train, X_valid, Y_valid)
    return model, scaler, nb_outputs
    

model, scaler, nb_outputs = test(64)
    

yarnNum = 1
stride = 1
skipFactor = 5        
vrtxNum = 300
firstFrame = 85
dataset = 'spacing0.5x_00011_stretch'
path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
frame0 = int(firstFrame/skipFactor)
frame1 = int(290/skipFactor + 1)
for i in range (frame0, frame1):
    f = i*skipFactor
    for y in range (0,yarnNum):
        X_test = np.loadtxt(path + "testX_" + str(f-firstFrame) + '_' + str(y) + ".txt",delimiter=None)
        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride)