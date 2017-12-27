
import numpy as np
import matplotlib.pyplot as plt
import mltools as ml
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout
from sklearn.preprocessing import MinMaxScaler
from keras.layers.normalization import BatchNormalization
import math

## Load data
# In[]
 
def loadData():

    X_train_all = np.loadtxt(path + "trainX_all.txt",delimiter=None)
    Y_train_all = np.loadtxt(path + "trainY_all.txt",delimiter=None)
    Y_train_all = Y_train_all[:, 0:3]
    
    #duplicate data
#    X_train_all = np.concatenate((X_train_all,X_train_all), axis=0)
#    X_train_all = np.concatenate((X_train_all,X_train_all), axis=0)
#    Y_train_all = np.concatenate((Y_train_all,Y_train_all), axis=0)
#    Y_train_all = np.concatenate((Y_train_all,Y_train_all), axis=0)
 
    print("Original training data shape (X): ", X_train_all.shape)
    print("Original training data shape (Y): ", Y_train_all.shape)
    
    #rescale the data
#    scaler = MinMaxScaler(feature_range=(0, 1))
#    scaler.fit(X_train_all)
#    X_train_all = scaler.transform(X_train_all)
#    print(scaler.data_max_)
#    scaler = MinMaxScaler(feature_range=(0, 1))
#    scaler.fit(X_test_all)
#    X_test_all = scaler.transform(X_test_all)

    scaler = MinMaxScaler(feature_range=(0, 1))
    scaler.fit(Y_train_all)
    Y_train_all = scaler.transform(Y_train_all)
    
    nb_features = X_train_all.shape[1]
    nb_traindata = X_train_all.shape[0]
    nb_halfdata = round(nb_traindata*0.95)
    nb_outputs = Y_train_all.shape[1]
    
    # using subset data as training and validation
    all_train = np.concatenate((X_train_all,Y_train_all), axis=1) 
    
    
    np.random.shuffle(all_train)
    X_train = all_train[0:nb_halfdata,0:nb_features]
    Y_train = all_train[0:nb_halfdata,nb_features:]   
    X_valid = all_train[nb_halfdata:,0:nb_features]
    Y_valid = all_train[nb_halfdata:,nb_features:] 
     
    # polynomio
    X_train_,params = ml.transforms.rescale(X_train);
    X_valid_,_ = ml.transforms.rescale( X_valid, params);
    X_test_,_ = ml.transforms.rescale( X_test, params);
    nb_features = X_train_.shape[1]
    
    # Represent the targets as one-hot vectors: e.g. 0 -> [1,0];  1 -> [0, 1].
    print("Training Y matrix shape: ", Y_train.shape)
    print("Validation Y matrix shape: ", Y_valid.shape)

     
#    return (X_train_, Y_train, X_valid_, Y_valid, nb_features,nb_outputs, X_test_)
    return (X_train, Y_train, X_valid, Y_valid, nb_features,nb_outputs, X_test, scaler)

## Build neural network model
# In[]:
def buildModel(input_dim, output_dim):

    # Simple fully-connected neural network with 2 hidden layers.
    # Including dropout layer helps avoid overfitting.
    model = Sequential()
    
#    model.add(Dense(64, input_dim=input_dim))
#    model.add(Activation('relu'))            
    
    model.add(Dense(64, input_dim=input_dim))
    model.add(Activation('relu'))   
    model.add(BatchNormalization())
    
    model.add(Dense(64))
    model.add(Activation('relu')) 
    model.add(BatchNormalization())
    
    model.add(Dense(output_dim))
    model.add(Activation('linear'))
    
    model.compile(optimizer='adam', loss='mse', metrics=['mse'])
    
    return model

## Train the model
# In[]
def trainModel(model, X_train, Y_train, X_valid, Y_valid):
    
    # Weights are updated one mini-batch at a time. A running average of the training loss is computed in real time, which is useful for identifying problems (e.g. the loss might explode or get stuck right). The validation loss is evaluated at the end of each epoch (without dropout).

    history = model.fit(X_train, Y_train, batch_size = 16, epochs = 100, verbose = 2,
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
def extrapolate(predicted, totalNum, filename, nb_outputs):
    total = np.zeros(totalNum*nb_outputs).reshape(totalNum,nb_outputs)
    assert(predicted.shape[1] == total.shape[1])
    r = math.floor (( totalNum - predicted.shape[0] )/2)
    total[r:r+predicted.shape[0], :] = predicted
    np.savetxt(path + filename, total, fmt='%.6f', delimiter=' ')
    
## Prediction
# In[]:
def predict(model, X_test, scaler, nb_outputs, filename):
    
    predicted = model.predict(X_test, verbose=0)
#    all_test = np.concatenate((X_test,predicted), axis=1)
    predicted = scaler.inverse_transform(predicted)
#    np.savetxt(path + 'testY_NN.txt', all_test[:, -3:], fmt='%.6f', delimiter=' ')
    np.savetxt(path + 'testY_NN.txt', predicted, fmt='%.6f', delimiter=' ')
    
    extrapolate(predicted, 300, filename, nb_outputs)
## Main
# In[]
path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/1220/NN/'
    
(X_train, Y_train, X_valid, Y_valid, nb_features, nb_outputs, X_test, scaler) = loadData()

model = buildModel(nb_features, nb_outputs)

model = trainModel(model, X_train, Y_train, X_valid, Y_valid)


## Test all frames
# In[]
frame0 = 3
frame1 = 36
for i in range (0 , frame1-frame0):
    f = i*5
    X_test = np.loadtxt(path + "trainX_" + str(f) + ".txt",delimiter=None)
    filename = "testY_NN_full" + str(f + frame0*5) + ".txt"
    predict(model, X_test, scaler, nb_outputs, filename)