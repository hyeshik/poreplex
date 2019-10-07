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

import h5py
from weighted_metrics import WeightedCategoricalCrossentropy, WeightedCategoricalAccuracy
from tensorflow.keras.layers import (
    Dense, Dropout, GaussianNoise, Bidirectional)
from tensorflow.keras.models import Sequential
from tensorflow.keras.callbacks import (
    Callback, EarlyStopping, CSVLogger, ModelCheckpoint, TensorBoard)
import tensorflow as tf
import shutil
from time import time
import pandas as pd
import numpy as np
import os
#import pickle


class CustomModelCheckpoint(Callback):

    def __init__(self, model, periodic, besttraining, bestvalidation, period=10):
        super().__init__()
        self.model_for_saving = model
        self.periodic_output = periodic
        self.besttraining_output = besttraining
        self.bestvalidation_output = bestvalidation
        self.period = period

        self.best_loss = self.best_val_loss = None

    def on_epoch_end(self, epoch, logs=None):
        if 'loss' not in logs or 'val_loss' not in logs: # maybe in termination
            return

        if self.best_loss is None or logs['loss'] < self.best_loss:
            self.model_for_saving.save_weights(self.besttraining_output, overwrite=True)
            self.best_loss = logs['loss']

        if self.best_val_loss is None or logs['val_loss'] < self.best_val_loss:
            self.model_for_saving.save_weights(self.bestvalidation_output, overwrite=True)
            self.best_val_loss = logs['val_loss']

        if (epoch + 1) % self.period == 0:
            modelpath = self.periodic_output.format(epoch=epoch)
            self.model_for_saving.save_weights(modelpath, overwrite=True)


def load_data(datapath, datatype):
    with h5py.File(datapath, 'r') as h5:
        print('Loading {} data from {}'.format(datatype, datapath))

        currentgroup = {}
        group = h5[datatype]
        for datasetname in ['signals', 'labels', 'weights', 'onehot', 'readid']:
            currentgroup[datasetname] = group[datasetname][:]

        print('  - done.')
        return currentgroup

def build_layers_LSTM1(input_shape, num_classes, cudnn=False):
    def lstm_layer(units, return_sequences):
        if cudnn:
            return tf.keras.layers.LSTM(units=units, return_sequences=return_sequences)
        else:
            return tf.keras.layers.RNN(tf.keras.layers.LSTMCell(units=units),
                                       return_sequences=return_sequences)

    model = Sequential()

    model.add(GaussianNoise(0.05, input_shape=input_shape))

    model.add(Bidirectional(lstm_layer(units=48, return_sequences=True)))
    model.add(Dropout(0.2))

    model.add(lstm_layer(units=64, return_sequences=False))
    model.add(Dropout(0.3))

    model.add(Dense(num_classes, activation='softmax'))

    return model

def build_layers_GRU1(input_shape, num_classes, cudnn=False):
    def gru_layer(units, return_sequences):
        if cudnn:
            return tf.keras.layers.GRU(units=units, return_sequences=return_sequences)
        else:
            return tf.keras.layers.RNN(tf.keras.layers.GRUCell(units=units),
                                       return_sequences=return_sequences)

    model = Sequential()

    model.add(GaussianNoise(0.05, input_shape=input_shape))

    model.add(Bidirectional(gru_layer(units=48, return_sequences=True)))
    model.add(Dropout(0.2))

    model.add(gru_layer(units=64, return_sequences=False))
    model.add(Dropout(0.3))

    model.add(Dense(num_classes, activation='softmax'))

    return model

def create_metric_functions(params, num_classes):
    cost_matrix = np.ones([num_classes, num_classes], dtype=np.float64)
    cost_matrix[1:, 1:] *= params['metrics']['cross_contamination_penalty']
    WeightedCategoricalAccuracy.name = 'w_acc'
    lossfun = WeightedCategoricalCrossentropy(cost_matrix)
    accfun = WeightedCategoricalAccuracy(cost_matrix)
    return lossfun, accfun

def create_model(params, layerdef, input_shape, num_classes):
    lossfun, accfun = create_metric_functions(params, num_classes)

    print('Creating model...')
    with tf.device('/cpu:0'):
        model = globals()['build_layers_' + layerdef](input_shape, num_classes, cudnn=True)

    print('Compiling...')
    strategy = tf.distribute.MirroredStrategy()
    with strategy.scope():
        pmodel = globals()['build_layers_' + layerdef](input_shape, num_classes, cudnn=True)
        pmodel.compile(loss=lossfun,
                       optimizer=params['optimizer'],
                       metrics=['accuracy', accfun])

    model.compile(loss=lossfun,
                  optimizer=params['optimizer'],
                  metrics=['accuracy', accfun])

    print('  - done.')

    return model, pmodel

def convert_model_to_noncudnn(layerdef, cudamodel, global_params, training_data, output_dir):
    num_classes = training_data['weights'].shape[0]
    input_shape = training_data['signals'].shape[1:]
    lossfun, accfun = create_metric_functions(global_params, num_classes)

    with tf.device('/cpu:0'):
        build_layers = globals()['build_layers_' + layerdef]
        cpu_model = build_layers(input_shape, num_classes, cudnn=False)
        cpu_model.set_weights(cudamodel.get_weights())

        cpu_model.compile(loss=lossfun,
                          optimizer=global_params['optimizer'],
                          metrics=['accuracy', accfun])

        cpu_model.fit(training_data['signals'],
                      training_data['onehot'],
                      batch_size=global_params['batchsize_train'],
                      validation_split=global_params['validation_split'],
                      epochs=1, verbose=1,
                      class_weight=training_data['weights'])

    cpu_model.save(os.path.join(output_dir, 'final-model.hdf5'))
    return cpu_model

def train_model(model, pmodel, global_params, training_data, output_dir):
    callbacks = [
        # Checkpoint saver
        CustomModelCheckpoint(model,
            periodic=output_dir + '/checkpoints-epoch{epoch:03d}.hdf5',
            besttraining=output_dir + '/bestmodel-training.hdf5',
            bestvalidation=output_dir + '/bestmodel-validation.hdf5',
            period=global_params['model_checkpoint_period']),
        # Logging
        CSVLogger(output_dir + '/training-log.csv'),
        # Early stop at stalled learning
        EarlyStopping(monitor='val_loss',
            min_delta=global_params['earlystopping_min_delta'],
            patience=global_params['earlystopping_patience'],
            verbose=1),
        TensorBoard(log_dir=output_dir + '/tensorboard')
    ]

    hist = pmodel.fit(training_data['signals'],
                      training_data['onehot'],
                      batch_size=global_params['batchsize_train'],
                      epochs=global_params['epochs'],
                      validation_split=global_params['validation_split'],
                      verbose=1, callbacks=callbacks,
                      class_weight=training_data['weights'])

    model.save(os.path.join(output_dir, 'final-cudamodel.hdf5'))


def evaluate_model(pmodel, global_params, testdata, output_dir):
    print('Evaluating the model performance...')
    loss, acc, w_acc = pmodel.evaluate(testdata['signals'], testdata['onehot'],
                                batch_size=global_params['batchsize_eval'],
                                verbose=False)
    print(loss, acc, w_acc, file=open(output_dir + '/evaluation.txt', 'w'), sep='\t')
    print('== Evaluation result ==')
    print(' * Test loss:', loss)
    print(' * Test accuracy:', acc)


def predict_classifications(pmodel, global_params, testdata, output_prefix):
    print('Making predictions for testing data...')
    predmtx = pmodel.predict(testdata['signals'],
                             batch_size=global_params['batchsize_test'],
                             verbose=False)
    np.save(output_prefix + '.npy', predmtx)

    predlabels = np.argmax(predmtx, axis=1)
    predlabel_probs = np.amax(predmtx, axis=1)

    summary_table = pd.DataFrame({
        'read_id': list(map(bytes.decode, testdata['readid'])),
        'prior_label': testdata['labels'],
        'pred_label': predlabels,
        'pred_prob': predlabel_probs,
    }).set_index('read_id')

    summary_table.to_csv(output_prefix + '.txt', sep='\t')


def main(global_params, layer_def, dataset_file, output_dir):
    training_data = load_data(dataset_file, 'training')

    model, pmodel = create_model(global_params, layer_def,
                                 training_data['signals'].shape[1:],
                                 training_data['weights'].shape[0])

    if os.path.isdir(output_dir):
        if any((not f.startswith('.')) for f in os.listdir(output_dir)):
            print('Clearing {}...'.format(output_dir))
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)
    else:
        os.makedirs(output_dir)

    train_model(model, pmodel, global_params, training_data, output_dir)

    print('Adopting the weights to a new model for CPU')
    cpu_model = \
        convert_model_to_noncudnn(layer_def, pmodel, global_params,
                                  training_data, output_dir)

    print('Evaluation of the CUDA model')
    predict_classifications(pmodel, global_params, training_data,
                            os.path.join(output_dir, 'training-predictions'))
    del training_data
    test_data = load_data(dataset_file, 'testing')
    predict_classifications(pmodel, global_params, test_data,
                            os.path.join(output_dir, 'test-predictions'))
    evaluate_model(pmodel, global_params, test_data, output_dir)

    print('Evaluation of the CPU model')
    evaluate_model(cpu_model, global_params, test_data, output_dir)


if __name__ == '__main__':
    import sys

    global_params, layer_def, dataset_file, output_dir = sys.argv[1:]

    import json

    global_params = json.loads(global_params)

    main(global_params, layer_def, dataset_file, output_dir)

