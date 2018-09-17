#!/usr/bin/env python3
#
# Copyright (c) 2018 Institute for Basic Science
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

from keras.layers import Dense, Dropout, LSTM, CuDNNLSTM
from keras.models import Sequential
from keras.callbacks import Callback, EarlyStopping, CSVLogger, ModelCheckpoint
from keras.utils.training_utils import multi_gpu_model
from numpy.lib.format import open_memmap
from scipy.stats import pearsonr
import shutil
import tensorflow as tf
import numpy as np
import pandas as pd
from time import time
import random
import pickle
import os

def build_layers(input_shape, output_variables, cudnn=False):
    if cudnn:
        lstm_layer = CuDNNLSTM
        lstm_options = {}
    else:
        lstm_layer = LSTM
        lstm_options = {'recurrent_activation': 'sigmoid'}

    model = Sequential()

    model.add(lstm_layer(48, return_sequences=True, input_shape=input_shape,
                         **lstm_options))
    model.add(Dropout(0.1))

    model.add(lstm_layer(48, return_sequences=False, **lstm_options))
    model.add(Dropout(0.2))

    model.add(Dense(output_variables, activation='linear'))
    return model

def create_training_model(params, input_shape, num_outputs):
    print('Creating model...')
    with tf.device('/cpu:0'):
        model = build_layers(input_shape, num_outputs, cudnn=True)

    print('Compiling...')

    if params['ngpu'] > 1:
        pmodel = multi_gpu_model(model, gpus=params['ngpu'])
        pmodel.compile(loss='mean_squared_error',
                       optimizer=params['optimizer'],
                       metrics=['mae'])
    else:
        pmodel = model

    model.compile(loss='mean_squared_error',
                  optimizer=params['optimizer'],
                  metrics=['mae'])

    print('  - done.')

    return model, pmodel

def convert_model_to_noncudnn(input_file, output_file, input_shape, num_outputs,
                              params):
    with tf.device('/cpu:0'):
        cpu_model = build_layers(input_shape, num_outputs, cudnn=False)

        cpu_model.compile(loss='mean_squared_error',
                          optimizer=params['optimizer'], metrics=['mae'])

    cpu_model.load_weights(input_file)
    cpu_model.save(output_file)
    return cpu_model

def train_model(model, pmodel, global_params, training_input, training_output, output_dir):
    callbacks = [
        # Logging
        CSVLogger(output_dir + '/training-log.csv'),
        # Save best models
        ModelCheckpoint(output_dir + '/checkpoint-best.hdf5', save_best_only=True,
                        monitor='val_loss', save_weights_only=True),
        # Early stop at stalled learning
        EarlyStopping(monitor='val_loss',
            min_delta=global_params['earlystopping_min_delta'],
            patience=global_params['earlystopping_patience'],
            verbose=1)
    ]

    hist = pmodel.fit(training_input, training_output,
                      batch_size=global_params['batchsize_train'],
                      epochs=global_params['epochs'],
                      validation_split=global_params['validation_split'],
                      verbose=1, callbacks=callbacks)

    model.save(output_dir + '/final-cudamodel.hdf5')

def evaluate_model(model, global_params, output_dir):
    testing_input = open_memmap(global_params['testing-x'])
    testing_output = open_memmap(global_params['testing-y'])

    xfrm_params = eval(open(global_params['transform-y']).read())

    predmtx = model.predict(testing_input, global_params['batchsize_test'], verbose=1)
    expected_mtx = np.array([
        (testing_output[:, 0] * xfrm_params['scale_std']) + xfrm_params['scale_mean'],
        (testing_output[:, 1] * xfrm_params['shift_std']) + xfrm_params['shift_mean']
    ]).T
    predmtx = np.array([
        (predmtx[:, 0] * xfrm_params['scale_std']) + xfrm_params['scale_mean'],
        (predmtx[:, 1] * xfrm_params['shift_std']) + xfrm_params['shift_mean']
    ]).T

    dt = np.hstack([expected_mtx, predmtx])

    print("\t\tscale\tshift")
    print("Pearson r\t{:.5f}\t{:.5f}".format(pearsonr(dt[:, 0], dt[:, 2])[0],
                                             pearsonr(dt[:, 1], dt[:, 3])[0]))
    print("RMSD\t\t{:.5f}\t{:.5f}".format(((dt[:, 0] - dt[:, 2]) ** 2).mean() ** 0.5,
                                          ((dt[:, 1] - dt[:, 3]) ** 2).mean() ** 0.5))

    np.save(os.path.join(output_dir, 'test-output.npy'), dt)

def main(global_params, output_dir):
    training_input = open_memmap(global_params['training-x'])
    training_output = open_memmap(global_params['training-y'])

    model, pmodel = create_training_model(global_params, training_input.shape[1:],
                                          training_output.shape[1])

    if os.path.isdir(output_dir):
        if any((not f.startswith('.')) for f in os.listdir(output_dir)):
            print('Clearing {}...'.format(output_dir))
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)
    else:
        os.makedirs(output_dir)

    train_model(model, pmodel, global_params, training_input, training_output, output_dir)

    print('Adopting the weights to a new model for CPU')
    cpu_model = \
        convert_model_to_noncudnn(os.path.join(output_dir, 'final-cudamodel.hdf5'),
                                  os.path.join(output_dir, 'final-model.hdf5'),
                                  training_input.shape[1:], training_output.shape[1],
                                  global_params)

    del training_input, training_output

    print('Evaluation of the CUDA model')
    evaluate_model(pmodel, global_params, output_dir)

    print('Evaluation of the CPU model')
    evaluate_model(cpu_model, global_params, output_dir)


if __name__ == '__main__':
    random.seed(922)

    global_params = {
        'ngpu': 2,
        'optimizer': 'adam',
        'earlystopping_min_delta': 0.0005,
        'earlystopping_patience': 30,
        'batchsize_train': 3072,
        'batchsize_test': 2048,
        'epochs': 300,
        'validation_split': 0.25,
        'transform-y': 'dataarrays/scaling-transform.txt',
        'training-x': 'dataarrays/signals-training.npy',
        'training-y': 'dataarrays/scaling-training.npy',
        'testing-x': 'dataarrays/signals-testing.npy',
        'testing-y': 'dataarrays/scaling-testing.npy',
    }

    main(global_params, 'trained-model')

