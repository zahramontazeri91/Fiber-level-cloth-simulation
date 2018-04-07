
import numpy as np
import matplotlib.pyplot as plt
#import mltools as ml
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout
from sklearn.preprocessing import MinMaxScaler
from keras.layers.normalization import BatchNormalization
import math
from shutil import copyfile
from keras.models import load_model

## Load data
# In[]
 
def loadData(fn_trainX, fn_trainY):
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
    split = 0.95
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
    
    model.add(Dense(neurons))
    model.add(Activation('relu')) 
    model.add(BatchNormalization())
    
#    model.add(Dense(neurons))
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
    history = model.fit(X_train, Y_train, batch_size = 16, epochs = 3, verbose = 2,
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
## Prediction
# In[]:
def predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, isRot):
    
    predicted = model.predict(X_test, verbose=0)
#    all_test = np.concatenate((X_test,predicted), axis=1)
    predicted = scaler.inverse_transform(predicted)
#    np.savetxt(path + 'testY_NN.txt', all_test[:, -3:], fmt='%.6f', delimiter=' ')
    
    # rotate the shape-match back to original frame (window-Rotation-minimizing to yarn-rotation-minimizing)
#    predicted = np.loadtxt(path + 'trainY_15000_0.txt')
    if isRot:
        angles = np.loadtxt(anglesFile, delimiter=None)
        predicted = rotate(predicted, angles)
    
    np.savetxt(path + 'testY_NN.txt', predicted, fmt='%.6f', delimiter=' ')
    
    extrapolate(predicted, vrtxNum, filename, nb_outputs, stride)
## Main
# In[]    
def test(neurons,fn_trainX, fn_trainY): 
    (X_train, Y_train, X_valid, Y_valid, nb_features, nb_outputs, scaler) = loadData(fn_trainX, fn_trainY)
    model = buildModel(nb_features, nb_outputs, neurons)
    model = trainModel(model, X_train, Y_train, X_valid, Y_valid)
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
datasets = []
config = 'pattern/yarn11/'
w_path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/' + config
datasets.append('spacing0.5x/10')
datasets.append('spacing0.5x/00011')
datasets.append('spacing0.5x/10100')
datasets.append('spacing0.5x/11110')
datasets.append('spacing1.0x/10')
datasets.append('spacing1.0x/00011')
datasets.append('spacing1.0x/10100')
datasets.append('spacing1.0x/11110')
datasets.append('spacing1.5x/10')
datasets.append('spacing1.5x/00011')
datasets.append('spacing1.5x/10100')
datasets.append('spacing1.5x/11110')
#datasets.append('spacing2.0x/10')
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

fn_trainX = w_path + "train_all/trainX_all.txt"
fn_trainY = w_path + "train_all/trainY_all.txt"

appendTrainingData(datasets, w_path, fn_trainX, fn_trainY)

model, scaler, nb_outputs = test(256, fn_trainX, fn_trainY)
model.save('savedNN/model_yarn11_ws9.h5')  # creates a HDF5 file 'my_model.h5'
del model  # deletes the existing model
model = load_model('savedNN/model_yarn11_ws9.h5')        

# In[] 
#yarnNum = 1
#stride = 1
#skipFactor = 2000        
#vrtxNum = 300
#firstFrame = 22000
#lastFrame = 99500
#dataset = 'twist/yarn4/damp2_500'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)

# In[] 
yarn0 = 10
yarn1 = 11
stride = 1
skipFactor = 500        
vrtxNum = 250
firstFrame = 10000
lastFrame = 74500
dataset = 'stretch/yarn4/stretch'
path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
frame0 = int(firstFrame/skipFactor)
frame1 = int(lastFrame/skipFactor + 1)
for i in range (frame0, frame1):
    f = i*skipFactor
    for y in range (yarn0, yarn1):
        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)


# In[] 
#yarnNum = 46
#stride = 1
#skipFactor = 500        
#vrtxNum = 250
#firstFrame = 0
#lastFrame = 8000
#dataset = 'woven/yarn4/spacing1.0x/00011'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)

# In[] 
#yarnNum = 1
#stride = 1
#skipFactor = 500        
#vrtxNum = 300
#
#firstFrame = 12000
#lastFrame = 14000
#dataset = 'pattern/yarn11/spacing0.5x/10'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        
#firstFrame = 14000
#lastFrame = 16000
#dataset = 'pattern/yarn11/spacing0.5x/00011'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        
#        
#firstFrame = 13000
#lastFrame = 15000
#dataset = 'pattern/yarn11/spacing0.5x/10100'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        
#firstFrame = 13000
#lastFrame = 15000
#dataset = 'pattern/yarn11/spacing0.5x/11110'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
## In[]
#firstFrame = 12500
#lastFrame = 14500
#dataset = 'pattern/yarn11/spacing1.0x/10'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        
#firstFrame = 15000
#lastFrame = 17000
#dataset = 'pattern/yarn11/spacing1.0x/00011'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        
#        
#firstFrame = 13500
#lastFrame = 15500
#dataset = 'pattern/yarn11/spacing1.0x/10100'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        
#firstFrame = 14000
#lastFrame = 16000
#dataset = 'pattern/yarn11/spacing1.0x/11110'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1) 
#
## In[]
#firstFrame = 13000
#lastFrame = 15000
#dataset = 'pattern/yarn11/spacing1.5x/10'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        
#firstFrame = 15500
#lastFrame = 17500
#dataset = 'pattern/yarn11/spacing1.5x/00011'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        
#        
#firstFrame = 14000
#lastFrame = 16000
#dataset = 'pattern/yarn11/spacing1.5x/10100'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        
#firstFrame = 14500
#lastFrame = 16500
#dataset = 'pattern/yarn11/spacing1.5x/11110'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1) 
        
# In[]
#firstFrame = 14000
#lastFrame = 16000
#dataset = 'pattern/yarn11/spacing2.0x/10'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        
#firstFrame = 16000
#lastFrame = 18000
#dataset = 'pattern/yarn11/spacing2.0x/00011'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        
#        
#firstFrame = 15000
#lastFrame = 17000
#dataset = 'pattern/yarn11/spacing2.0x/10100'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        
#firstFrame = 15500
#lastFrame = 17500
#dataset = 'pattern/yarn11/spacing2.0x/11110'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1) 
#
## In[]
#firstFrame = 14500
#lastFrame = 16500
#dataset = 'pattern/yarn11/spacing2.5x/10'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        
#firstFrame = 16500
#lastFrame = 18500
#dataset = 'pattern/yarn11/spacing2.5x/00011'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        
#        
#firstFrame = 15500
#lastFrame = 17500
#dataset = 'pattern/yarn11/spacing2.5x/10100'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        
#firstFrame = 16500
#lastFrame = 18500
#dataset = 'pattern/yarn11/spacing2.5x/11110'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1) 
#
## In[]
#firstFrame = 14500
#lastFrame = 16500
#dataset = 'pattern/yarn11/spacing3.0x/10'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        
#firstFrame = 17000
#lastFrame = 19000
#dataset = 'pattern/yarn11/spacing3.0x/00011'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        
#        
#firstFrame = 18000
#lastFrame = 20000
#dataset = 'pattern/yarn11/spacing3.0x/10100'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1)
#        
#firstFrame = 17000
#lastFrame = 19000
#dataset = 'pattern/yarn11/spacing3.0x/11110'
#path = 'D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/input/'+dataset+'/NN/'
#frame0 = int(firstFrame/skipFactor)
#frame1 = int(lastFrame/skipFactor + 1)
#for i in range (frame0, frame1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        X_test = np.loadtxt(path + "testX_" + str(f) + '_' + str(y) + ".txt",delimiter=None)
#        filename = "testY_NN_full_" + str(f) + '_' + str(y) +  ".txt"
#        anglesFile = path + "angles_" + str(f) + '_' + str(y) + ".txt"
#        predict(model, X_test, scaler, nb_outputs, filename, vrtxNum, stride, anglesFile, 1) 