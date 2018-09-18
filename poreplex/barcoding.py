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

import numpy as np
import os
from .keras_wrap import keras

class BarcodeDemultiplexer:

    def __init__(self, config):
        self.config = config
        self.model = self.load_model()
        self.signals = []
        self.signal_assoc_read = []

    def clear(self):
        del self.signals[:]
        del self.signal_assoc_read[:]

    def load_model(self):
        model_file = os.path.join(os.path.dirname(__file__),
                        'presets', self.config['demux_model'])
        return keras.models.load_model(model_file)

    def push(self, npread, signal):
        minlen = self.config['minimum_dna_length']
        maxlen = self.config['maximum_dna_length']

        if not minlen <= len(signal) <= maxlen:
            return

        trimlength = self.config['signal_trim_length']
        if len(signal) > trimlength:
            signal = signal[-trimlength:]
        elif len(signal) < trimlength:
            signal = np.pad(signal, (trimlength - len(signal), 0), 'constant')

        self.signals.append(signal)
        self.signal_assoc_read.append(npread)

    def predict(self):
        if len(self.signals) > 0:
            signals = np.array(self.signals)[:, :, np.newaxis]
            predweights = self.model.predict(signals,
                batch_size=self.config['maximum_batch_size'], verbose=0)
            predlabels = np.argmax(predweights, axis=1) - self.config['number_of_decoy_labels']
            predlabel_probs = np.amax(predweights, axis=1)

            predvalid = ((predlabels >= 0) &
                         (predlabel_probs >= self.config['pred_weight_cutoff']))
            for npread, is_valid, bcid in zip(self.signal_assoc_read, predvalid, predlabels):
                if is_valid:
                    npread.set_barcode(int(bcid))

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
