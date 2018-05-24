#!/usr/bin/env python3
#
# Copyright (c) 2018 Hyeshik Chang
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
from keras.layers import Dense, Dropout, GRU, LSTM
from keras.models import Sequential
from keras.callbacks import Callback
from keras.utils.training_utils import multi_gpu_model
import tensorflow as tf
import shutil
from time import time
import pandas as pd
import numpy as np
import os
#import pickle


class CustomModelCheckpoint(Callback):

    def __init__(self, model, modelpath, logpath, period=10):
        super().__init__()
        self.modelpath = modelpath
        self.model_for_saving = model
        self.period = period
        self.logfile = open(logpath, 'w')

    def on_epoch_end(self, epoch, logs=None):
        loss = logs['val_loss']
        accuracy = logs['acc']

        print(time(), epoch, logs['loss'], logs['acc'],
              logs['val_loss'], logs['val_acc'], sep='\t', file=self.logfile)
        self.logfile.flush()
        if epoch % self.period == 0 and epoch > 0:
            modelpath = self.modelpath.format(epoch=epoch)
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


def create_model(params, layerdef, input_shape, num_classes):
    print('Creating model...')

    with tf.device('/cpu:0'):
        model = Sequential()

        for i, defrow in enumerate(layerdef):
            celltype = eval(defrow[0])
            param1 = defrow[1]
            kwargs = {}

            if i == 0:
                kwargs['input_shape'] = input_shape
                kwargs['return_sequences'] = True

            model.add(celltype(param1, **kwargs))

        model.add(Dense(num_classes, activation='softmax'))

    print('Compiling...')

    if params['ngpu'] > 1:
        pmodel = multi_gpu_model(model, gpus=params['ngpu'])
        pmodel.compile(loss='categorical_crossentropy',
                       optimizer=params['optimizer'],
                       metrics=['accuracy'])
    else:
        pmodel = model

    model.compile(loss='categorical_crossentropy',
                  optimizer=params['optimizer'],
                  metrics=['accuracy'])

    print('  - done.')

    return model, pmodel

def train_model(model, pmodel, global_params, training_data, output_dir):
    modelcheckpoint = CustomModelCheckpoint(model,
        output_dir + '/checkpoints-epoch{epoch:03d}.hdf5',
        output_dir + '/training-log.txt',
        global_params['model_checkpoint_period'])

    hist = pmodel.fit(training_data['signals'],
                      training_data['onehot'],
                      batch_size=global_params['batchsize_train'],
                      epochs=global_params['epochs'],
                      validation_split=global_params['validation_split'],
                      verbose=1,
                      callbacks=[modelcheckpoint],
                      class_weight=training_data['weights'])

    model.save(output_dir + '/model.hdf5')


def evaluate_model(pmodel, global_params, testdata, output_dir):
    print('Evaluating the model performance...')
    score, acc = pmodel.evaluate(testdata['signals'], testdata['onehot'],
                                 batch_size=global_params['batchsize_eval'])
    print(score, acc, file=open(output_dir + '/evaluation.txt', 'w'), sep='\t')
    print('== Evaluation result ==')
    print(' * Test score:', score)
    print(' * Test accuracy:', acc)


def predict_test_classifications(pmodel, global_params, testdata, output_dir):
    print('Making predictions for testing data...')
    predmtx = pmodel.predict(testdata['signals'],
                             batch_size=global_params['batchsize_test'],
                             verbose=1)
    np.save(output_dir + '/test-prediction-output.npy', predmtx)

    predlabels = np.argmax(predmtx, axis=1)
    predlabel_probs = np.amax(predmtx, axis=1)

    summary_table = pd.DataFrame({
        'read_id': list(map(bytes.decode, testdata['readid'])),
        'prior_label': testdata['labels'],
        'pred_label': predlabels,
        'pred_prob': predlabel_probs,
    }).set_index('read_id')

    summary_table.to_csv(output_dir + '/test-predictions.txt', sep='\t')


def main(global_params, layer_def, dataset_file, output_dir):
    training_data = load_data(dataset_file, 'training')

    model, pmodel = create_model(global_params, layer_def,
                                 training_data['signals'].shape[1:],
                                 training_data['weights'].shape[0])

    if os.path.isdir(output_dir):
        if any((not f.startswith('.')) for f in os.listdir(output_dir)):
            answer = input('Output directory {} exists. Clear it? (y/N) '.format(output_dir))
            if answer.strip().lower().startswith('y'):
                print('Clearing {}...'.format(output_dir))
                shutil.rmtree(output_dir)
                os.mkdir(output_dir)
            else:
                return
    else:
        os.makedirs(output_dir)

    train_model(model, pmodel, global_params, training_data, output_dir)

    del training_data


    testdata = load_data(dataset_file, 'testing')
    predict_test_classifications(pmodel, global_params, testdata, output_dir)
    evaluate_model(pmodel, global_params, testdata, output_dir)


if __name__ == '__main__':
    #global_params = '{"ngpu": 2, "epochs": 12, "validation_split": 0.1, "batchsize_train": 512, "batchsize_eval": 64, "batchsize_test": 64, "optimizer": "adam", "model_checkpoint_period": 10}'
    #layer_def = '[["GRU", 64], ["Dropout", 0.5], ["GRU", 64], ["Dropout", 0.5]]'
    #dataset_file = '../traindata/signals-MXG3.1-s2000-t500.hdf5'
    #output_dir = '../models/MXG3.1-s2000-t500-lGRU_64_64'
    import sys

    global_params, layer_def, dataset_file, output_dir = sys.argv[1:]

    import json

    global_params = json.loads(global_params)
    layer_def = json.loads(layer_def)

    main(global_params, layer_def, dataset_file, output_dir)
