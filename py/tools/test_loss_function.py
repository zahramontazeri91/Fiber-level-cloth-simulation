from keras.layers import Input, Dense
from keras.models import Model
from keras import backend as K
import numpy as np

def customLoss(pnt_ref='', pnt_simul='', angles='', nbatch=16):
    def loss(y_true, y_pred):   
        err = 0
        for i in range (0, nbatch):
            #Here we want to do something with each y_pred
            pnt = np.random.rand(2,1)
            pnt_K = K.variable(value=pnt) #convert the point to tensor
            y_pred_mat = K.reshape(y_pred[i], [-1,2])
            new_pnt = K.dot(y_pred_mat, pnt_K)
            err = err + K.mean( K.square(new_pnt-pnt_K) )
            print(i, y_pred[i])
        return err
    return loss


inputs = Input(shape=(1,))
preds = Dense(4,activation='linear')(inputs) #output 4 values

model = Model(inputs=inputs,outputs=preds)
model.compile(loss=customLoss(nbatch=16), optimizer='adam' ,metrics=['mse'])
#model.fit(x,y, batch_size=1, epochs=30, shuffle=False)