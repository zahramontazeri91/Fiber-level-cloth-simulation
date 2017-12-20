
import numpy as np
import matplotlib.pyplot as plt
import mltools as ml
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout

## Load data
# In[] 
def loadData():
    #    X_train_all = np.random.rand(10,2)
    #    Y_train_all = np.random.rand(10,1)
    #    X_test_all = np.random.rand(10,2)
    path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/NN/'
    X_train_all = np.loadtxt(path + "trainX_all.txt",delimiter=None)
    Y_train_all = np.loadtxt(path + "trainY_all.txt",delimiter=None)
    Y_train_all = Y_train_all[:, 0:3]
    X_test_all = np.loadtxt(path + "testX_0.txt",delimiter=None)
    print("Original training data shape (X): ", X_train_all.shape)
    print("Original training data shape (Y): ", Y_train_all.shape)
    
    
    nb_features = X_train_all.shape[1]
    nb_traindata = X_train_all.shape[0]
    nb_halfdata = round(nb_traindata*0.6)
    nb_outputs = Y_train_all.shape[1]
    
    # using subset data as training and validation
    X_train = X_train_all[0:nb_halfdata,:]
    Y_train = Y_train_all[0:nb_halfdata]
#    X_train = X_train_all
#    Y_train = Y_train_all
    X_valid = X_train_all[nb_halfdata:,:]
    Y_valid = Y_train_all[nb_halfdata:]
    X_test = X_test_all
     
    # polynomio
    X_train_,params = ml.transforms.rescale(X_train);
    X_valid_,_ = ml.transforms.rescale( X_valid, params);
    X_test_,_ = ml.transforms.rescale( X_test, params);
    nb_features = X_train_.shape[1]
    
    # Represent the targets as one-hot vectors: e.g. 0 -> [1,0];  1 -> [0, 1].
    print("Training Y matrix shape: ", Y_train.shape)
    print(Y_train[0:10])
    print("Validation Y matrix shape: ", Y_valid.shape)
    print(Y_valid[0:10])
     
    return (X_train_, Y_train, X_valid_, Y_valid, nb_features,nb_outputs, X_test_)
#    return (X_train, Y_train, X_valid, Y_valid, nb_features,nb_outputs, X_test)

## Build neural network model
# In[]:
def buildModel(input_dim, output_dim):

    print(input_dim, output_dim)
    # Simple fully-connected neural network with 2 hidden layers.
    # Including dropout layer helps avoid overfitting.
    model = Sequential()
    
    model.add(Dense(4, input_dim=input_dim)) # Use input_shape=(28,28) for unflattened data.
    model.add(Activation('relu'))    
#    model.add(Dropout(0.1));
    
#    model.add(Dense(32)) # Use input_shape=(28,28) for unflattened data.
#    model.add(Activation('relu'))
#    model.add(Dropout(0.1));
#        
#    model.add(Dense(256)) # Use input_shape=(28,28) for unflattened data.
#    model.add(Activation('relu'))
#    model.add(Dropout(0.1));

    
    model.add(Dense(output_dim))
    model.add(Activation('linear'))
    # Use softmax layer for multi-class problems.
    
    model.compile(optimizer='sgd', loss='mse', metrics=['accuracy'])
    
    return model

## Train the model
# In[]
def trainModel(model, X_train, Y_train, X_valid, Y_valid):
    
    # Weights are updated one mini-batch at a time. A running average of the training loss is computed in real time, which is useful for identifying problems (e.g. the loss might explode or get stuck right). The validation loss is evaluated at the end of each epoch (without dropout).

    history = model.fit(X_train, Y_train, batch_size = 16, epochs = 500, verbose = 2,
                        validation_data=(X_valid, Y_valid))
        
    # Plot loss trajectory throughout training.
    plt.figure()
    plt.subplot(1,2,1)
    plt.plot(history.history['loss'], label='train')
    plt.plot(history.history['val_loss'], label='valid')
    plt.xlabel('Epoch')
    plt.ylabel('Cross-Entropy Loss')
    plt.legend()
    plt.show()
    
    plt.subplot(1,2,2)
    plt.plot(history.history['acc'], label='train')
    plt.plot(history.history['val_acc'], label='valid')
    plt.xlabel('Epoch')
    plt.ylabel('Accuracy')
    plt.legend()
    plt.show()
    
    
    # ## Evaluate performance    
    # Note: when calling evaluate, dropout is automatically turned off.
    score = model.evaluate(X_valid, Y_valid, verbose=0)
    print('Validation cross-entropy loss: %0.5f' % score[0])
    print('Validation accuracy: %0.2f' % score[1])

    return model

## Prediction
# In[]:
def predict(model, X_test):
    
    predicted = model.predict(X_test, verbose=0)
#    print(predicted)
#    ID = np.arange(0,X_test.shape[0])
#    np.savetxt('data/Y_test.txt', np.c_[ID.conj().T, predicted[:,1]], fmt='%i, %1.2f', delimiter=',')
    np.savetxt('testResult.txt', np.vstack( (np.arange(len(predicted)) , predicted[:,0]) ).T, '%d, %.2f',header='ID,Prob1',comments='',delimiter=',');

## Main
# In[]
(X_train, Y_train, X_valid, Y_valid, nb_features, nb_outputs, X_test) = loadData()

model = buildModel(nb_features, nb_outputs)

model = trainModel(model, X_train, Y_train, X_valid, Y_valid)

predict(model, X_test)
