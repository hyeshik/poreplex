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

import os
import sys
from io import StringIO

__all__ = ['tensorflow', 'keras', 'WeightedCategoricalCrossentropy',
           'WeightedCategoricalAccuracy']

try: # Suppress warnings and informative messages from keras and tf.
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    sys.stderr, saved_stderr = StringIO(), sys.stderr
    from tensorflow import keras
    import tensorflow
    tensorflow.get_logger().setLevel('ERROR')
finally:
    sys.stderr = saved_stderr


# Weighted versions of Crossentropy and Accuracy metrics from eliadi:
#  https://github.com/keras-team/keras/issues/2115#issuecomment-530762739

import tensorflow.keras.backend as K
from tensorflow.keras.losses import CategoricalCrossentropy
from tensorflow.keras.metrics import CategoricalAccuracy

class WeightedCategoricalCrossentropy(CategoricalCrossentropy):
    def __init__(self, cost_mat, name='weighted_categorical_crossentropy', **kwargs):
        assert(cost_mat.ndim == 2)
        assert(cost_mat.shape[0] == cost_mat.shape[1])

        super().__init__(name=name, **kwargs)
        self.cost_mat = K.cast_to_floatx(cost_mat)

    def __call__(self, y_true, y_pred):
        return super().__call__(
            y_true=y_true,
            y_pred=y_pred,
            sample_weight=get_sample_weights(y_true, y_pred, self.cost_mat),
        )

def get_sample_weights(y_true, y_pred, cost_m):
    num_classes = len(cost_m)

    #y_pred.shape.assert_has_rank(2)
    #y_pred.shape[1].assert_is_compatible_with(num_classes)
    #y_pred.shape.assert_is_compatible_with(y_true.shape)

    y_pred = K.one_hot(K.argmax(y_pred), num_classes)

    y_true_nk1 = K.expand_dims(y_true, 2)
    y_pred_n1k = K.expand_dims(y_pred, 1)
    cost_m_1kk = K.expand_dims(cost_m, 0)

    sample_weights_nkk = cost_m_1kk * y_true_nk1 * y_pred_n1k
    sample_weights_n = K.sum(sample_weights_nkk, axis=[1, 2])

    return sample_weights_n

class WeightedCategoricalAccuracy(CategoricalAccuracy):
    def __init__(self, cost_mat, name='weighted_categorical_accuracy', **kwargs):
        assert(cost_mat.ndim == 2)
        assert(cost_mat.shape[0] == cost_mat.shape[1])

        super().__init__(name=name, **kwargs)
        self.cost_mat = K.cast_to_floatx(cost_mat)

    def update_state(self, y_true, y_pred, sample_weight=None):
        return super().update_state(
            y_true=y_true,
            y_pred=y_pred,
            sample_weight=get_sample_weights(y_true, y_pred, self.cost_mat),
        )

