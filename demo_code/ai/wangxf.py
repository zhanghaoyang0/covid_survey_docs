import xml.dom.minidom
import numpy as np
import pandas as pd
import sys
import os
from scipy.stats.mstats import winsorize
import tensorflow as tf
from tensorflow import keras
from keras_transformer import get_encoders
from keras_transformer import gelu
from tensorflow.keras import backend as K
from sklearn.metrics import f1_score
from python_utils import converters
from tensorflow.keras.models import Sequential,load_model
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score,fbeta_score, roc_auc_score, roc_curve, roc_curve, auc
from sklearn import metrics
import neurokit2 as nk
import random

os.environ['TF_KERAS']='1' #其中key和value均为string类型
os.putenv('TF_KERAS', '1')
os.environ.setdefault('TF_KERAS', '1')

####################数据增强
#基准校正
def ecg_signal_baseline_process(all_ecg_data): # [x, 12,5000]
    for sample_id in range(0, all_ecg_data.shape[0]):
        ecg_data = all_ecg_data[sample_id]
        #基准校正
        for lead_id in range(0, 12):
            ecg_data[lead_id] = nk.ecg_clean(ecg_data[lead_id], sampling_rate=500, method='neurokit')
    return all_ecg_data

#加载数据集，同时进行打乱操作
def data_augment(X_all, Y_all, add_number):
    permutation = list(np.random.permutation(add_number))
    X_augment = X_all[permutation, :].copy()
    Y_augment = Y_all[permutation].copy()
    lens = X_augment.shape[2]
    for i in range(0,add_number):
        location = np.random.randint(-40, 40)
        #print("loact", location)
        if(location>0): #左移
            X_augment[i,:,0: lens-location] = X_augment[i,:,location:]
            for row in range(0, X_augment.shape[1]):
                X_augment[i,row,lens-location:] = X_augment[i,row, lens-location]
        elif (location<=0): #右移
            location = location * -1
            X_augment[i,:, location:] = X_augment[i,:,0: lens-location]
            for row in range(0,X_augment.shape[1]):
                X_augment[i,row, 0:location] = X_augment[i][row][location]
    X_all = np.concatenate((X_all, X_augment), axis=0)
    Y_all = np.concatenate((Y_all, Y_augment), axis=0)
    return X_all, Y_all

def add_gauss(X):
    for i in range(X.shape[0]):
        for row in range(X.shape[1]):
            for column in range(X.shape[2]):
                X[i,row,column] = X[i,row,column] + random.gauss(0,1)
    return X

############################
def FCN_bone(input_layber, features, ksize):
    # conv1 = keras.layers.Conv1D(filters=features, kernel_size=ksize, padding='same')(input_layber)
    # conv1 = keras.layers.BatchNormalization()(conv1)
    # conv1 = keras.layers.Activation(activation='relu')(conv1)
    conv1 = keras.layers.Conv1D(filters=features, kernel_size=ksize, padding='same')(input_layber)
    conv1 = keras.layers.BatchNormalization()(conv1)
    conv1 = keras.layers.Activation(activation='relu')(conv1)
    conv1 = keras.layers.MaxPool1D()(conv1)
    # conv1 = keras.layers.Conv1D(filters=features, kernel_size=ksize, padding='same')(conv1)
    # conv1 = keras.layers.BatchNormalization()(conv1)
    # conv1 = keras.layers.Activation('relu')(conv1)
    return  conv1

def Dense(shape):
    inputlayer = keras.layers.Input(shape)
    dense_layer = keras.layers.Dense(units=512, activation='relu')(inputlayer)
    #dense_layer = keras.layers.Dropout(0.3)(dense_layer)
    #dense_layer = keras.layers.Dense(units=128, activation='relu')(dense_layer)
    #dense_layer = keras.layers.Dropout(0.5)(dense_layer)
    dense_layer = keras.layers.Dense(units=256, activation='relu')(dense_layer)
    dense_layer = keras.layers.Dropout(0.5)(dense_layer)
    dense_layer = keras.layers.Dense(units=32, activation='relu')(dense_layer)
    outputlayer = keras.layers.Dense(1, activation='sigmoid')(dense_layer)
    model = keras.Model(inputs=inputlayer, outputs=outputlayer)
    model.compile(loss=keras.losses.BinaryCrossentropy(), optimizer=keras.optimizers.Adam(lr=0.001),
                  metrics=[keras.metrics.BinaryAccuracy(name='accuracy'),
                           keras.metrics.AUC(name='auc'),
                           keras.metrics.Recall(name='Recall')
                         ])
    model.summary()
    return model

def FCN(shape):
    inputlayer = keras.layers.Input(shape)
    conv1 = FCN_bone(inputlayer, 48, 8)
    conv1 = FCN_bone(conv1, 64, 5)
    conv1 = FCN_bone(conv1, 128, 5)
    conv1 = FCN_bone(conv1, 256, 3)
    gap_layer = keras.layers.GlobalAveragePooling1D()(conv1)
    dense_layer = keras.layers.Dense(units=32, activation='relu')(gap_layer)
    outputlayer = keras.layers.Dense(1, activation='sigmoid')(dense_layer)
    model = keras.Model(inputs=inputlayer, outputs=outputlayer)
    model.compile(loss=keras.losses.BinaryCrossentropy(), optimizer=keras.optimizers.Adam(lr=0.001),
                  metrics=[keras.metrics.BinaryAccuracy(name='accuracy'),
                           keras.metrics.AUC(name='auc'),
                           keras.metrics.Recall(name='Recall')
                         ])
    model.summary()
    return model
    
def Inception_bone(input_layer, n_feature_maps):
    conv_x = keras.layers.Conv1D(filters=n_feature_maps, kernel_size=1, padding='same')(input_layer)
    conv_x = keras.layers.BatchNormalization()(conv_x)
    conv_x = keras.layers.Activation('relu')(conv_x)
    conv_y = keras.layers.Conv1D(filters=n_feature_maps, kernel_size=3, padding='same')(input_layer)
    conv_y = keras.layers.BatchNormalization()(conv_y)
    conv_y = keras.layers.Activation('relu')(conv_y)
    conv_z = keras.layers.Conv1D(filters=n_feature_maps, kernel_size=5, padding='same')(input_layer)
    conv_z = keras.layers.BatchNormalization()(conv_z)
    conv_z = keras.layers.Activation('relu')(conv_z)
    cnn_concat = keras.layers.concatenate([conv_x,conv_y, conv_z])
    return  cnn_concat

def Inception_cnn_1D(shape,  feature_num):
    n_feature_maps = feature_num
    input_layer = keras.layers.Input(name='the_input', shape=shape, dtype='float32')  # (None, 128, 64, 1)
    conv1 = FCN_bone(input_layer, 48, 8)
    conv1 = FCN_bone(conv1, 64, 5)
    incep_bone = Inception_bone(conv1, n_feature_maps)
    cnn_concat = keras.layers.Conv1D(filters=n_feature_maps*3, kernel_size=(3), padding='same')(incep_bone)
    cnn_concat = keras.layers.BatchNormalization()(cnn_concat)
    cnn_concat = keras.layers.Activation('relu')(cnn_concat)
    cnn_concat = keras.layers.MaxPool1D()(cnn_concat)
    gap_layer = keras.layers.GlobalAveragePooling1D()(cnn_concat)
    gap_layer = keras.layers.Dense(units=32, activation='relu')(gap_layer)
    # gap_layer = keras.layers.Dropout(dropout)(gap_layer)
    output_layer = keras.layers.Dense(1, activation='sigmoid')(gap_layer)
    model = keras.models.Model(inputs=input_layer, outputs=output_layer)
    model.compile(loss=keras.losses.BinaryCrossentropy(), optimizer=keras.optimizers.Adam(lr=0.001),
                  metrics=[keras.metrics.BinaryAccuracy(name='accuracy'),
                           keras.metrics.AUC(name='auc'),
                           keras.metrics.Recall(name='Recall')
                         ])
    model.summary()
    return model

def Attention_encoder_model(shape):
    input_layer = keras.layers.Input(shape)
    conv1 = FCN_bone(input_layer, 64, 8)
    conv1 = FCN_bone(conv1, 128, 5)
    # conv block -1
    conv3 = keras.layers.Conv1D(filters=256, kernel_size=5, padding='same')(conv1)
    conv3 = keras.layers.BatchNormalization()(conv3)
    conv3 = keras.layers.PReLU(shared_axes=[1])(conv3)
    conv3 = keras.layers.MaxPooling1D()(conv3)
    # split for attention
    attention_data = keras.layers.Lambda(lambda x: x[:, :, :128])(conv3)
    attention_softmax = keras.layers.Lambda(lambda x: x[:, :, 128:])(conv3)
    # attention mechanism
    attention_softmax = keras.layers.Softmax()(attention_softmax)
    multiply_layer = keras.layers.Multiply()([attention_softmax, attention_data])
    # last layer
    gap_layer = keras.layers.GlobalAveragePooling1D()(multiply_layer)
    dense_layer = keras.layers.Dense(units=32, activation='relu')(gap_layer)
    # output layer
    output_layer = keras.layers.Dense(units=1, activation='sigmoid')(dense_layer)
    model = keras.models.Model(inputs=input_layer, outputs=output_layer)
    model.compile(loss=keras.losses.BinaryCrossentropy(), optimizer=keras.optimizers.Adam(lr=0.001),
                  metrics=[keras.metrics.BinaryAccuracy(name='accuracy'),
                           keras.metrics.AUC(name='auc'),
                           keras.metrics.Recall(name='Recall')
                         ])
    model.summary()
    return model
   
def LSTM_Model(shape, feature_num):
    # Make Networkw
    inputs = keras.layers.Input(shape)
    conv1 = FCN_bone(inputs, 64, 5)
    # CNN to RNN
    lstm_1 = keras.layers.LSTM(128, return_sequences=True, kernel_initializer='he_normal', name='lstm1')(
        conv1)  # (None, 32, 512)
    lstm_2 = keras.layers.LSTM(128, return_sequences=True, kernel_initializer='he_normal', name='lstm2')(lstm_1)#lstm1_merged
    # transforms RNN output to character activations:
    conv1 = FCN_bone(lstm_2, 256, 3)
    gap_layer = keras.layers.GlobalAveragePooling1D()(conv1)
    gap_layer = keras.layers.Dense(32,activation= "relu")(gap_layer)
    output_layer = keras.layers.Dense(1, activation='sigmoid')(gap_layer)  # (None, 32, 63)
    model = tf.keras.models.Model(inputs=inputs, outputs=output_layer)
    model.compile(loss=keras.losses.BinaryCrossentropy(), optimizer=keras.optimizers.Adam(lr=0.001),
                  metrics=[keras.metrics.BinaryAccuracy(name='accuracy'),
                           keras.metrics.AUC(name='auc'),
                           keras.metrics.Recall(name='Recall')
                         ])
    model.summary()
    return model

def Attention(input_layer):
    conv1 = FCN_bone(input_layer, 64, 8)
    conv1 = FCN_bone(conv1, 128, 5)
    # conv block -1
    conv3 = keras.layers.Conv1D(filters=256, kernel_size=5, padding='same')(conv1)
    conv3 = keras.layers.BatchNormalization()(conv3)
    conv3 = keras.layers.PReLU(shared_axes=[1])(conv3)
    conv3 = keras.layers.MaxPooling1D()(conv3)
    # split for attention
    attention_data = keras.layers.Lambda(lambda x: x[:, :, :128])(conv3)
    attention_softmax = keras.layers.Lambda(lambda x: x[:, :, 128:])(conv3)
    # attention mechanism
    attention_softmax = keras.layers.Softmax()(attention_softmax)
    multiply_layer = keras.layers.Multiply()([attention_softmax, attention_data])
    # last layer
    gap_layer = keras.layers.GlobalAveragePooling1D()(multiply_layer)
    dense_layer = keras.layers.Dense(units=32, activation='relu')(gap_layer)
    return dense_layer

def LSTM_Attention_Model(shape):
    # Make Networkw
    inputs = keras.layers.Input(shape)
    input_atten = keras.layers.Input(shape)
    atten_dense = Attention(input_atten)
    conv1 = FCN_bone(inputs, 64, 5)
    # CNN to RNN
    lstm_1 = keras.layers.LSTM(128, return_sequences=True, kernel_initializer='he_normal', name='lstm1')(
        conv1)  # (None, 32, 512)
    lstm_2 = keras.layers.LSTM(128, return_sequences=True, kernel_initializer='he_normal', name='lstm2')(lstm_1)#lstm1_merged
    # transforms RNN output to character activations:
    conv1 = FCN_bone(lstm_2, 256, 3)
    gap_layer = keras.layers.GlobalAveragePooling1D()(conv1)
    gap_layer = keras.layers.Dense(32,activation= "relu")(gap_layer)
    contac_layer = keras.layers.add([gap_layer,atten_dense])
    #contac_layer = keras.layers.concatenate([gap_layer,atten_dense], axis=-1)
    dense = keras.layers.Dense(8, activation= 'relu')(contac_layer)
    output_layer = keras.layers.Dense(1, activation='sigmoid')(dense)  # (None, 32, 63)
    model = tf.keras.models.Model(inputs=[inputs,input_atten], outputs=output_layer)
    model.compile(loss=keras.losses.BinaryCrossentropy(), optimizer=keras.optimizers.Adam(lr=0.001),
                  metrics=[#keras.metrics.BinaryAccuracy(name='accuracy'),
                           keras.metrics.AUC(name='auc'),
                           keras.metrics.Recall(name='Recall')
                         ])
    model.summary()
    return model

def residual_bone(input_layer,ksize, feature_num, stride = 2):
    n_feature_maps = feature_num
    # BLOCK 1
    conv_x = keras.layers.Conv1D(filters=n_feature_maps, padding='valid', kernel_size=ksize, strides= stride)(input_layer)#
    conv_x = keras.layers.BatchNormalization()(conv_x)
    conv_x = keras.layers.Activation('relu')(conv_x)
    conv_y = keras.layers.Conv1D(filters=n_feature_maps*2, padding='valid', kernel_size=ksize,strides= stride )(conv_x)
    conv_y = keras.layers.BatchNormalization()(conv_y)
    conv_y = keras.layers.Activation('relu')(conv_y)
    #print("conv_y shape", conv_y.shape)
    # conv_y = keras.layers.Conv1D(filters=n_feature_maps*2, padding='valid', kernel_size=ksize,strides= stride )(conv_y)
    # conv_y = keras.layers.BatchNormalization()(conv_y)
    # conv_y = keras.layers.Activation('relu')(conv_y)
    output_block = keras.layers.MaxPool1D()(conv_y)
    return output_block

def residual_bone_nostride(input_layer,ksize, feature_num,is_pool):
    n_feature_maps = feature_num
    # BLOCK 1
    conv_x = keras.layers.Conv1D(filters=n_feature_maps, padding='same', kernel_size=ksize)(input_layer)#
    conv_x = keras.layers.BatchNormalization()(conv_x)
    conv_x = keras.layers.Activation('relu')(conv_x)
    print("conv_x shape", conv_x.shape)
    conv_y = keras.layers.Conv1D(filters=n_feature_maps, padding='same', kernel_size=ksize)(conv_x)
    conv_y = keras.layers.BatchNormalization()(conv_y)
    conv_y = keras.layers.Activation('relu')(conv_y)
    print("conv_y shape", conv_y.shape)
    conv_z = keras.layers.Conv1D(filters=n_feature_maps, padding='same', kernel_size=ksize)(conv_y) #padding = "same"
    conv_z = keras.layers.BatchNormalization()(conv_z)
    # expand channels for the sum
    shortcut_y = keras.layers.Conv1D(filters=n_feature_maps, kernel_size=1, padding='same')(input_layer)
    shortcut_y = keras.layers.BatchNormalization()(shortcut_y)
    output_block = keras.layers.add([shortcut_y, conv_z])
    output_block = keras.layers.Activation('relu')(conv_z)
    if is_pool:
        output_block = keras.layers.MaxPool1D()(output_block)
    # BLOCK 2
    return output_block
 
def Resnet_model_group(shape):
    n_feature_maps = 128
    input_shape = (shape)
    input_list = []
    flatten_group = []
    for i in range(0, 12):
        input_layer = keras.layers.Input(input_shape)
        input_list.append(input_layer)
        net_layer = residual_bone(input_layer, 5, n_feature_maps, 1)
        group_layer = residual_bone(net_layer, 3, n_feature_maps, 1)
        #group_layer = residual_bone(group_layer, 3, n_feature_maps, 1)
        #flatten_layer = keras.layers.Flatten()(group_layer)
        flatten_layer = keras.layers.GlobalAveragePooling1D()(group_layer)
        flatten_group.append(flatten_layer)
    # input_layer1 = keras.layers.Input(input_shape)
    # input_list.append(input_layer1)
    # # input_layer2 = keras.layers.Input(input_shape)
    # group1 = residual_bone(input_list[0],3, n_feature_maps,1)
    # flatten_layer1 = keras.layers.GlobalAveragePooling1D()(group1)
    concatenated = keras.layers.concatenate(flatten_group)
    #flatten_layer = keras.layers.Flatten()(concatenated)
    #dense_layer1 = keras.layers.Dense(units=256, activation='relu')(concatenated)
    #dense_layer1 = keras.layers.Dense(units=256, activation='relu')(concatenated)
    dense_layer2 = keras.layers.Dense(units=64, activation='relu')(concatenated)
    output_layer = keras.layers.Dense(1, activation='sigmoid')(dense_layer2)
    model = keras.models.Model(inputs=input_list, outputs=output_layer)
    model.compile(loss=keras.losses.BinaryCrossentropy(), optimizer=keras.optimizers.Adam(lr=0.001),
                  metrics=[keras.metrics.BinaryAccuracy(name='accuracy'),
                           keras.metrics.AUC(name='auc'),
                           keras.metrics.Recall(name='Recall')
                         ])
    model.summary()
    return model

def Resnet_model(shape, ksize, feature_num):
    input_layer = keras.layers.Input(shape)
    #inner = residual_bone(input_layer, 3, 128, 0)
    inner = residual_bone(input_layer, 3, 128, 1)
    inner = residual_bone(inner, 3, 128, 1)
    dense_layer = keras.layers.GlobalAveragePooling1D()(inner)
    dense_layer = keras.layers.Dense(units=32, activation='relu')(dense_layer)
    #dense_layer = keras.layers.Dense(units=8, activation='relu')(dense_layer)
    #dense_layer = keras.layers.Dropout(0.2)(dense_layer)
    output_layer = keras.layers.Dense(1, activation='sigmoid')(dense_layer) #softmax, sigmoid
    # output layer
    model = keras.models.Model(inputs=input_layer, outputs=output_layer)
    model.compile(loss=keras.losses.BinaryCrossentropy(), optimizer=keras.optimizers.Adam(lr=0.001),
                  metrics=[keras.metrics.BinaryAccuracy(name='accuracy'),
                           keras.metrics.AUC(name='auc'),
                           keras.metrics.Recall(name='Recall')
                         ])
    model.summary()
    return model

def transformer_encoder(shape):
    input_layer = keras.layers.Input(shape)
    conv1 = FCN_bone(input_layer, 48, 8)
    conv1 = FCN_bone(conv1, 64, 5)
    conv1 = FCN_bone(conv1, 120, 5)
    encoded_layer = get_encoders(
        encoder_num=8,
        input_layer=conv1,
        head_num=12,
        hidden_dim=24,
        attention_activation='relu',
        feed_forward_activation=gelu,
        dropout_rate= 0,
    )
    #if (is_cnn_afterencoder==True):
    encoded_layer = keras.layers.Reshape((encoded_layer.shape[2],encoded_layer.shape[1]))(encoded_layer)
    gap_layer = keras.layers.GlobalAveragePooling1D()(encoded_layer)
    # pool = keras.layers.MaxPool1D()(encoded_layer)
    # flatten_layer = keras.layers.Flatten()(pool)
    gap_layer = keras.layers.Dense(units=32, activation='relu')(gap_layer)
    #gap_layer = keras.layers.Dropout(0)(gap_layer)
    output_layer = keras.layers.Dense(1, activation='sigmoid')(gap_layer) #softmax, sigmoid
    # output layer
    model = keras.models.Model(inputs=input_layer, outputs=output_layer)
    model.compile(loss=keras.losses.BinaryCrossentropy(), optimizer=keras.optimizers.Adam(lr=0.0001),
                  metrics=[keras.metrics.BinaryAccuracy(name='accuracy'),
                           keras.metrics.AUC(name='auc'),
                           keras.metrics.Recall(name='Recall')
                         ])
    model.summary()
    return model

def Transfer_learn_resnet():
    model_ptb_path = "/data2/users/wangxinfeng/UKB_ECG/PTB/select_sample_rating/one_wave_data/Model/"
    model = load_model(model_ptb_path + "ptb_RES_0_cv.h5")
    print(model.summary())
    # return model
    for (i, layer) in enumerate(model.layers):
        print(i, layer)
        #layer = model.layers[i]
        #layer.trainable = False
    # print(model.summary())
    # for (i) in range(0, 27): #
    #     layer = model.layers[i]
    #     layer.trainable = False
    #ptb_model = tf.keras.Model(inputs=model.input, outputs=model.get_layer('dense').output) # 'dense' #'global_average_pooling1d'
    ptb_model = tf.keras.Model(inputs=model.input, outputs=model.layers[-3].output) # 'dense' #'global_average_pooling1d'
    #model2 = tf.keras.Model(inputs=model.input, outputs=model.layers[-2].output)
    print("------------------test\n")
    #model2 = Sequential()
    #for layer in model.layers[:-1]:  # 跳过最后一层
    #    model2.add(layer)
    new_model = Sequential()
    # model2.summary()
    new_model.add(ptb_model)
    # new_model.add(keras.layers.Conv1D(filters=128, padding='same', kernel_size=5))
    # #new_model.add(keras.layers.Conv1D(keras.layers.BatchNormalization()))
    # new_model.add(keras.layers.Conv1D(keras.layers.Activation('relu')))
    # new_model.add(keras.layers.GlobalAveragePooling1D())
    new_model.add(keras.layers.Dense(32, activation='relu'))
    new_model.add(keras.layers.Dense(1, activation='sigmoid'))
    new_model.compile(loss=keras.losses.BinaryCrossentropy(), optimizer=keras.optimizers.Adam(lr=0.0001),
                  metrics=[keras.metrics.BinaryAccuracy(name='accuracy'),
                           keras.metrics.AUC(name='auc')
                         ])
    new_model.summary()
    return  new_model

GPU_num = 4#int(sys.argv[1])
Type_num = 0#int(sys.argv[2])
gpus = tf.config.experimental.list_physical_devices(device_type='GPU')
print(gpus)
if gpus:
    tf.config.experimental.set_visible_devices(devices=(gpus[GPU_num]), device_type='GPU')
    tf.config.experimental.set_memory_growth(gpus[GPU_num], True)
Add_sex_age = 0

def divide_train_val_by_kfold(k_fold, fold_sum, X_all,Y_all,Sex_age_Label,Sample_info):
    #5拆交叉取数据
    num_val_samples = int(X_all.shape[0] // fold_sum)
    X_val = X_all[num_val_samples * k_fold:num_val_samples * (k_fold + 1)]
    X_train = np.concatenate((X_all[:num_val_samples * k_fold], X_all[num_val_samples * (k_fold + 1):]),axis =0)
    y_val = Y_all[num_val_samples * k_fold:num_val_samples * (k_fold + 1)]
    y_train = np.concatenate((Y_all[:num_val_samples * k_fold], Y_all[num_val_samples * (k_fold + 1):]),axis =0)
    sex_val = Sex_age_Label[num_val_samples * k_fold:num_val_samples * (k_fold + 1)]
    sex_train = np.concatenate((Sex_age_Label[:num_val_samples * k_fold], Sex_age_Label[num_val_samples * (k_fold + 1):]),axis =0)
    sample_info_val = Sample_info.iloc[num_val_samples * k_fold:num_val_samples * (k_fold + 1)]
    #sample_id_train = np.concatenate((Sex_age_Label[:num_val_samples * k_fold], Sex_age_Label[num_val_samples * (k_fold + 1):]),axis =0)
    #X_train, X_val = preprocess_signals(X_train, X_val, outputfolder)
    print(X_train.shape, y_train.shape, X_val.shape, y_val.shape, sex_train.shape, sex_val.shape)
    return (X_train, y_train, X_val, y_val, sex_train, sex_val, sample_info_val)

def load_dataset(name,lead_id = 0):
    df_label = pd.read_csv(outputfolder + name+"_all_id.csv", delimiter=",", index_col=0, header=0)
    all_label = df_label["label"].values
    sex_age = df_label.iloc[:, [0, 1]]
    Sample_info = df_label
    sex_age_arr = sex_age.values
    diease_data = np.load(outputfolder+name+"_ecg_wave_460.npy", allow_pickle=True)
    #diease_data = np.reshape(diease_data,[diease_data.shape[0], diease_data.shape[1]* diease_data.shape[2]])
    print("diease_data shape:", diease_data.shape)
    #diease_data = ecg_signal_baseline_process(diease_data)
    #diease_data.dump(outputfolder+name+"_ecg_baseline.npy")
    noise_data = np.load(outputfolder+name+"_noise_wave_460.npy", allow_pickle=True)
    #noise_data = np.reshape(noise_data, [noise_data.shape[0], noise_data.shape[1] * noise_data.shape[2]])
    #noise_data = ecg_signal_baseline_process(noise_data)
    #noise_data.dump(outputfolder+name+"_noise_ecg_baseline.npy")
    #合并
    X_all = np.concatenate((diease_data, noise_data), axis=0)
    #以整个数据集来归一化不太好，以样本归一化试下
    Y_all = all_label#np.concatenate((part_diease_lable, part_noise_lable), axis=0)
    #发现一个问题，如果要交换两个轴，应该用a.swapaxes(1,2)
    #X_all = X_all.swapaxes(1,2)
    #print("x_all swapaxes", X_all.shape)
    #wave_len = X_all.shape[2] * X_all.shape[3]
    #把两头连接处置为0，方便网络确认
    #X_all = X_all[:,:, 0:6] = 0
    # X_all = X_all.swapaxes(1, 2)
    # for i in range(0, X_all.shape[0]):
    #     mean = X_all[i].mean()
    #     X_all[i] -= mean
    #     std = X_all[i].std()
    #     X_all[i] /= std
    # X_all = X_all.swapaxes(1, 2)
    #X_all = np.reshape(X_all, (X_all.shape[0], X_all.shape[1], 1))
    X_all = tf.reshape(X_all, (X_all.shape[0], X_all.shape[2],X_all.shape[1])) #如果是二维array,交换的不是转置操作。
    #X_all = X_all.swapaxes(2, 3)
    X_all = np.array(X_all)
    ##再将合并后的数据打乱
    np.random.seed(1000)
    permutation = list(np.random.permutation(X_all.shape[0]))
    X_shuf = X_all[permutation, :]
    Y_shuf = Y_all[permutation] # all_lable
    Sex_age_Label = sex_age_arr[permutation, :]
    Sample_info = Sample_info.iloc[permutation, :]
    return X_shuf,Y_shuf,Sex_age_Label,Sample_info

def Save_CV_predict_result(sample_info_val, y_pred, Net_id):
    #分段保存预测数值
    sample_info_val["PredictLabel"] = y_pred
    if (os.path.exists(outputfolder+"CV_result/"+Net_name[Net_id]+ "_cv.csv")):
        X_exist =  pd.read_csv(outputfolder +"CV_result/"+Net_name[Net_id]+"_cv.csv", delimiter=",", index_col = 0,header = 0)
        X_exist = pd.concat((X_exist, sample_info_val))
        X_exist.to_csv(outputfolder +"CV_result/"+Net_name[Net_id]+"_cv.csv")
    else:
        sample_info_val.to_csv(outputfolder+"CV_result/"+Net_name[Net_id]+"_cv.csv")

def run_model_in_k_fold(k_fold_sum, fold_index, Net_Id):
    (X_train, y_train, X_val, y_val, sex_train, sex_val, sample_info_val) = \
                      divide_train_val_by_kfold(fold_index, k_fold_sum, X_all, Y_all, Sex_age_Label,Sample_info)
    results = []
    reduce_lr = tf.keras.callbacks.ReduceLROnPlateau(
        monitor='val_auc', factor=0.1, patience=4, verbose=1, mode='max',
        min_delta=0.0001, cooldown=0, min_lr=0)
    early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_auc', mode='max', verbose=1, patience=10)
    # if Net_Id == 8:
    #     model = LSTM_Attention_Model(X_train[0].shape)
    #     model.fit([X_train,X_train], y_train, batch_size=batchsize, epochs=epochs_num, validation_data=([X_val,X_val], y_val),
    #               callbacks=[reduce_lr, early_stop],verbose=1)  # callbacks=[reduce_lr, early_stop],validation_split = 0.2,
    #     result = model.evaluate([X_val,X_val], y_val)
    #     results.append(result[2])
    #     y_val_pred = model.predict([X_val,X_val])
    #
    #     result = f1_score(y_val, y_val_pred.round(), average='micro')
    #     results.append(result)
    #     return results
    if Net_Id == 0:
        model = FCN(X_train[0].shape)
    elif Net_Id == 1:
        model = Inception_cnn_1D(X_train[0].shape, feature_num[0])
    elif Net_Id == 2:
        model = Resnet_model(X_train[0].shape, ksize[0], feature_num[0])
    elif Net_Id == 3 :
        model = LSTM_Model(X_train[0].shape, feature_num[0])
    elif Net_Id == 4:
        model = Attention_encoder_model(X_train[0].shape)
    elif Net_Id == 5:
        model = transformer_encoder(X_train[0].shape)
    elif Net_Id == 6:
        model = Dense(X_train[0].shape)
    elif Net_Id == 7:
        model = Transfer_learn_resnet()
    model.fit(X_train, y_train, batch_size= batchsize, epochs = epochs_num,validation_data=(X_val,y_val),
                        callbacks=[reduce_lr,early_stop], verbose = 1)  # callbacks=[reduce_lr, early_stop],validation_split = 0.2,
    # #if(Net_Id == 2):
    # #model.save(outputfolder +"/Model/"+ Net_name[Net_Id]+"_"+str(fold_index)+"_cv.h5")
    result = model.evaluate(X_val,y_val)
    results.append(result[2])
    y_val_pred = model.predict(X_val)
    result = f1_score(y_val,y_val_pred.round(),average='micro')
    results.append(result)
    #if (Net_Id == 2):
    Save_CV_predict_result(sample_info_val, y_val_pred, Net_Id)
    return results

def apply_standardizer(X_all):
    ss = StandardScaler()
    ss.fit(np.vstack(X_all).flatten()[:, np.newaxis].astype(float))
    X_tmp = []
    for x in X_all:
        x_shape = x.shape
        X_tmp.append(ss.transform(x.flatten()[:, np.newaxis]).reshape(x_shape))
    X_tmp = np.array(X_tmp)
    return X_tmp


Net_name = ["CNN", "INCep", "RES", "LSTM", "Atten", "Transf","Dense", "Transfer_Ptb", "subgroup"]
diease_name = "I251"#"I20-I251"
sample_type = "/select_sample_rating/"#random_select_sample select_sample_rating
outputfolder = "/data2/users/wangxinfeng/UKB_ECG/"+diease_name+sample_type + "one_wave_data/"
model_path = "/data2/users/wangxinfeng/UKB_ECG/"+diease_name+sample_type+"Model/"
k_fold_sum = 5
ksize = [[8, 5, 3]]  # ,[3,3,3][8,5,3],,[17, 11, 9, 5, 3]
feature_num = [64]  # 64,
batchsize = 64
epochs_num = 5

for netid in range(2,3):
    for lead_id in range(0, 12):
        (X_all,Y_all,Sex_age_Label, Sample_info) = load_dataset(diease_name, 0)
        print("X_all :", X_all.shape)    #X_all shape (2672,5000,12)
        Net_Id = netid #Type_num
        k_fold_sum = 5
        for fold_index in range(0,k_fold_sum):#k_fold_sum
            keras.backend.clear_session()
            print(fold_index, k_fold_sum)
            results = run_model_in_k_fold(k_fold_sum, fold_index,Net_Id)
            file_proce = open(outputfolder+"output/"+Net_name[Net_Id]+".txt", "a")
            pro_loss = str("k_flod_num: " +str(fold_index) + ",AUC_F1:" + str(results) + "\n")
            #pro_loss = str("lead: "+str(lead_id)+", k_flod_num: " +str(fold_index) + ",AUC_F1:" + str(results) + "\n")
            file_proce.write(pro_loss)
            file_proce.close()