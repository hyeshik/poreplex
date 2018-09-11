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
import numpy as np
import pandas as pd
import os
from functools import partial
from .keras_wrap import keras
from ont_fast5_api.fast5_file import Fast5File

class SignalLoader:

    def __init__(self, config, fast5prefix):
        self.config = config
        self.fast5prefix = fast5prefix
        self.scaler_model, self.scaler_cfg = self.load_scaler_model()
        self.head_signals = []
        self.head_signal_readids = []

    def clear(self):
        del self.head_signals[:]
        del self.head_signal_readids[:]

    def load_scaler_model(self):
        model_file = os.path.join(os.path.dirname(__file__),
                        'presets', self.config['scaler_model'])

        with h5py.File(model_file, 'r') as modelh5:
            mw = modelh5['model_weights'].attrs
            scaler_cfg = eval(mw['input_defs'].decode())
            scaler_cfg['model_version'] = mw['model_version'].decode()

            xfrm = eval(mw['output_transform'].decode())
            scaler_cfg['xfrm_scale'] = np.poly1d([xfrm['scale_std'], xfrm['scale_mean']])
            scaler_cfg['xfrm_shift'] = np.poly1d([xfrm['shift_std'], xfrm['shift_mean']])

        return keras.models.load_model(model_file), scaler_cfg

    def prepare_loading(self, filename):
        fast5path = os.path.join(self.fast5prefix, filename)

        fast5 = Fast5File(fast5path, 'r')
        if len(fast5.status.read_info) != 1:
            return {
                'filename': filename,
                'status': 'irregular_fast5',
            }

        read_info = fast5.status.read_info[0]
        channel_info = fast5.get_channel_info()
        tracking_id = fast5.get_tracking_id()

        sampling_rate = channel_info['sampling_rate']

        metainfo = {
            'filename': filename,
            'read_id': read_info.read_id,
            'channel': channel_info['channel_number'],
            'start_time': round(read_info.start_time / sampling_rate, 3),
            'run_id': tracking_id['run_id'],
            'sample_id': tracking_id['sample_id'],
            'duration': read_info.duration,
            'num_events': read_info.event_data_count,
            'sequence_length': 0,
            'mean_qscore': 0.,
        }

        context = {
            'read_info': read_info,
            'channel_info': channel_info,
            'tracking_id': tracking_id,
            'fast5': fast5,
            'sampling_rate': sampling_rate,
            'meta': metainfo,
            'status': '',
            'full_raw_signal': None,
        }

        error_status = self.load_padded_signal_head(fast5, read_info)
        if error_status:
            context['status'] = error_status

        return context

    def load_padded_signal_head(self, fast5, read_info):
        length_limit = self.scaler_cfg['length']
        stride = self.scaler_cfg['stride']

        sigload_length = min(length_limit, read_info.duration)
        sigload_length = sigload_length - sigload_length % stride

        signal = fast5.get_raw_data(end=sigload_length, scale=True)
        if len(signal) % stride > 0:
            signal = signal[:-(len(signal) % stride)]

        if len(signal) < self.scaler_cfg['min_length']:
            return 'scaler_signal_too_short'

        signal_means = signal.reshape([len(signal) // stride, stride]
                                      ).mean(axis=1, dtype=np.float32)
        length_limit //= stride
        if len(signal_means) < length_limit:
            signal_means = np.pad(signal_means, [length_limit - len(signal_means), 0],
                                  'constant')

        # Add the processed signal and read id to the queue
        self.head_signals.append(signal_means)
        self.head_signal_readids.append(read_info.read_id)

    def fit_scalers(self, contexts):
        batch_size = self.config['scaler_maximum_batch_size']
        scfg = self.scaler_cfg

        predmtx = self.scaler_model.predict(
                        np.array(self.head_signals)[:, :, np.newaxis], batch_size)
        scaling_params = np.transpose([
            scfg['xfrm_scale'](predmtx[:, 0]), scfg['xfrm_shift'](predmtx[:, 1])])

        for read_id, paramset in zip(self.head_signal_readids, scaling_params):
            context = contexts[read_id]
            paramset = np.array(paramset, dtype=scfg['dtype'])

            f_load = partial(load_signal, context, paramset, scfg['stride'])

            context['signal_scaling'] = paramset
            context['load_signal'] = f_load


def load_signal(context, scaling, stride=None, end=None, pool=True, pad=False):
    # Load from the cache if available
    if context['full_raw_signal'] is not None:
        sig = context['full_raw_signal']
    else:
        sig = context['fast5'].get_raw_data(end=end, scale=True)
        if end is None:
            context['full_raw_signal'] = sig

    if pool:
        cutend = len(sig) - len(sig) % stride
        sig = sig[:cutend].reshape([len(sig) // stride, stride]
                                   ).mean(axis=1, dtype=np.float32)
    else:
        stride = 1

    if end is not None:
        expected_size = end // stride
        if len(sig) > expected_size:
            sig = sig[-expected_size:]
        elif pad and len(sig) < expected_size:
            sig = np.pad(sig, [expected_size - len(sig), 0], 'constant')

    return np.poly1d(scaling)(sig), stride


if __name__ == '__main__':
    import yaml
    config = yaml.load(open('poreplex/presets/rna-r941.cfg'))['signal_processing']
    scaler = SignalLoader(config, 'testinput/live')

    from ont_fast5_api.fast5_file import Fast5File
    import os
    inpfiles = [f for f in os.listdir('testinput/live') if f.endswith('.fast5')]

    for inpfile in inpfiles:
        print('Loading', inpfile)
        scaler.prepare_loading(inpfile)

    scaler.determine_scaling()
    print(next(iter(scaler.contexts.values())))

