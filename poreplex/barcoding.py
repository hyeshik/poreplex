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

from functools import partial
from bisect import bisect_right
import numpy as np
import h5py
import os
from .keras_wrap import keras, WeightedCategoricalAccuracy, WeightedCategoricalCrossentropy

class BarcodeDemultiplexer:

    PAD_FILLER = -1000.

    def __init__(self, config, qualitythreshold):
        self.config = config
        self.calibration_table = [-1]
        self.model = self.load_model()
        self.signals = []
        self.signal_assoc_read = []

        if len(self.calibration_table) - 1 < qualitythreshold:
            raise ValueError('The current demultiplexer does not support calibrated score '
                             'of {}. Consider lowering --barcoding-quality-filter value.'
                             .format(qualitythreshold))
        self.score_threshold = self.calibration_table[qualitythreshold]

    def clear(self):
        del self.signals[:]
        del self.signal_assoc_read[:]

    def load_model(self):
        model_file = os.path.join(os.path.dirname(__file__),
                        'presets', self.config['demux_model'])

        with h5py.File(model_file, 'r') as modf:
            calibtable = modf['poreplex_params/calibration']
            if np.any(calibtable['phred'][:] != np.arange(len(calibtable))):
                raise RuntimeError('Calibration table in {} is not continuous.'.format(
                                        model_file))
            self.calibration_table = list(calibtable['pred_score'])

            cost_mtx = modf['poreplex_params/loss_weights'][:]
            keras.utils.get_custom_objects().update({
                'WeightedCategoricalAccuracy': partial(WeightedCategoricalAccuracy,
                                                       cost_mat=cost_mtx),
                'WeightedCategoricalCrossentropy': partial(WeightedCategoricalCrossentropy,
                                                           cost_mat=cost_mtx),
            })

        return keras.models.load_model(model_file)

    def lookup_calibrated_phred_score(self, score):
        if score <= 0.:
            return 0
        return bisect_right(self.calibration_table, score)

    @staticmethod
    def normalize_signal(sig):
        med = np.median(sig)
        mad = np.median(np.abs(sig - med))
        return (sig - med) / max(0.01, (mad * 1.4826))

    def push(self, npread, signal):
        minlen = self.config['minimum_dna_length']
        maxlen = self.config['maximum_dna_length']

        if not minlen <= len(signal) <= maxlen:
            return

        trimlength = self.config['signal_trim_length']
        if len(signal) > trimlength:
            signal = self.normalize_signal(signal[-trimlength:])
        elif len(signal) < trimlength:
            signal = np.pad(self.normalize_signal(signal),
                            (trimlength - len(signal), 0), 'constant',
                            constant_values=self.PAD_FILLER)
        else:
            signal = self.normalize_signal(signal)

        self.signals.append(signal)
        self.signal_assoc_read.append(npread)

    def predict(self):
        if len(self.signals) > 0:
            signals = np.array(self.signals)[:, :, np.newaxis]
            predweights = self.model.predict(signals,
                batch_size=self.config['maximum_batch_size'], verbose=0)
            predlabels = np.argmax(predweights, axis=1) - self.config['number_of_decoy_labels']
            predscores = np.amax(predweights, axis=1)

            for npread, bcid, score in zip(self.signal_assoc_read, predlabels, predscores):
                if bcid >= 0 and score >= self.score_threshold:
                    effective_bcid = int(bcid)
                else:
                    effective_bcid = None

                calib_score = self.lookup_calibrated_phred_score(score)
                npread.set_barcode(effective_bcid, int(bcid), calib_score)


if __name__ == '__main__':
    import yaml
    config = yaml.load(open('poreplex/presets/rna-r941.cfg'))
    demux = BarcodeDemultiplexer(config)

    import h5py
    with h5py.File('training/traindata/signals-MXG3.1-s2000-t300.hdf5', 'r') as h5:
        for i, (readid, signal) in enumerate(zip(h5['testing/readid'], h5['testing/signals'])):
            signal = signal[:, 0]
            demux.push(readid.decode(), signal)
            if i > 150:
                break

    print(demux.predict())
