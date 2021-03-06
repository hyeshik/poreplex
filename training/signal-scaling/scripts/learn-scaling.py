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
    def lstm_layer(units, return_sequences):
        if cudnn:
            return tf.keras.layers.LSTM(units=units, return_sequences=return_sequences)
        else:
            return tf.keras.layers.RNN(tf.keras.layers.LSTMCell(units=units),
                                       return_sequences=return_sequences)

    model = tf.keras.models.Sequential()

    model.add(tf.keras.layers.GaussianNoise(1.5, input_shape=input_shape))

    model.add(tf.keras.layers.Bidirectional(lstm_layer(units=48, return_sequences=True)))
    model.add(tf.keras.layers.Dropout(0.1))

    model.add(lstm_layer(units=48, return_sequences=False))
    model.add(tf.keras.layers.Dropout(0.2))

    model.add(tf.keras.layers.BatchNormalization())

    model.add(tf.keras.layers.Dense(output_variables, activation='linear'))
    return model

def create_training_model(params, input_shape, num_outputs):
    print('Compiling...')
    strategy = tf.distribute.MirroredStrategy()
    with strategy.scope():
        model = build_layers(input_shape, num_outputs, cudnn=True)
        model.compile(loss='mean_squared_error',
                      optimizer=params['optimizer'],
                      metrics=['mae'])

    print('  - done.')

    return model

def convert_model_to_noncudnn(cudamodel, global_params, training_input,
                              training_output, output_dir):
    num_outputs = training_output.shape[1]
    input_shape = training_input.shape[1:]

    with tf.device('/cpu:0'):
        cpu_model = build_layers(input_shape, num_outputs, cudnn=False)
        cpu_model.set_weights(cudamodel.get_weights())

        cpu_model.compile(loss='mean_squared_error',
                          optimizer=global_params['optimizer'],
                          metrics=['mae'])

        cpu_model.fit(training_input, training_output,
                      batch_size=global_params['batchsize_train'],
                      validation_split=global_params['validation_split'],
                      epochs=1, verbose=1)

    cpu_model.save(os.path.join(output_dir, 'final-model.hdf5'))
    return cpu_model

def train_model(model, global_params, training_input, training_output, output_dir):
    callbacks = [
        # Logging
        tf.keras.callbacks.CSVLogger(os.path.join(output_dir, 'training-log.csv')),
        # Save best models
        tf.keras.callbacks.ModelCheckpoint(os.path.join(output_dir, 'checkpoint-best.hdf5'),
                save_best_only=True, monitor='val_loss', save_weights_only=True),
        # Early stop at stalled learning
        tf.keras.callbacks.EarlyStopping(monitor='val_loss',
            min_delta=global_params['earlystopping_min_delta'],
            patience=global_params['earlystopping_patience'],
            verbose=1)
    ]

    hist = model.fit(training_input, training_output,
                     batch_size=global_params['batchsize_train'],
                     epochs=global_params['epochs'],
                     validation_split=global_params['validation_split'],
                     verbose=1, callbacks=callbacks)

    model.save(os.path.join(output_dir, 'final-cudamodel.hdf5'))

def evaluate_model(model, global_params, output_dir):
    testing_input = open_memmap(global_params['testing-x'])
    testing_output = open_memmap(global_params['testing-y'])

    xfrm_params = eval(open(global_params['transform-y']).read())

    predmtx = model.predict(testing_input, global_params['batchsize_test'])
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

    model = create_training_model(global_params, training_input.shape[1:],
                                  training_output.shape[1])

    if os.path.isdir(output_dir):
        if any((not f.startswith('.')) for f in os.listdir(output_dir)):
            print('Clearing {}...'.format(output_dir))
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)
    else:
        os.makedirs(output_dir)

    train_model(model, global_params, training_input, training_output, output_dir)

    print('Adopting the weights to a new model for CPU')
    cpu_model = \
        convert_model_to_noncudnn(model, global_params, training_input,
                                  training_output, output_dir)

    del training_input, training_output

    print('Evaluation of the CUDA model')
    evaluate_model(model, global_params, output_dir)

    print('Evaluation of the CPU model')
    evaluate_model(cpu_model, global_params, output_dir)


if __name__ == '__main__':
    random.seed(922)

    global_params = {
        'optimizer': 'adam',
        'earlystopping_min_delta': 0.001,
        'earlystopping_patience': 30,
        'batchsize_train': 1024,
        'batchsize_test': 1024,
        'epochs': 500,
        'validation_split': 0.25,
        'transform-y': 'dataarrays/scaling-transform.txt',
        'training-x': 'dataarrays/signals-training.npy',
        'training-y': 'dataarrays/scaling-training.npy',
        'testing-x': 'dataarrays/signals-testing.npy',
        'testing-y': 'dataarrays/scaling-testing.npy',
    }

    main(global_params, 'trained-model')

