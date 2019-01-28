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
from scipy.stats import norm
from functools import partial
from .keras_wrap import keras
from .signal_analyzer import SignalAnalysisError
from .fast5_file import Fast5Reader

__all__ = ['SignalLoader', 'NanoporeRead']

class SignalLoader:

    def __init__(self, config, fast5prefix):
        self.config = config
        self.fast5prefix = fast5prefix
        self.scaler_model, self.scaler_cfg = self.load_scaler_model()
        self.head_signals = []
        self.head_signal_assoc_reads = []

    def clear(self):
        del self.head_signals[:]
        del self.head_signal_assoc_reads[:]

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

            # Compute the bounds of predicted outcomes for passing QC. Since the training
            # data has 1.8 times boosted stdev of the original, the chance of falling
            # out of this filter should be lower than the setted threshold.
            qc_threshold_level = self.config['scaler_qc_threshold']
            qc_threshold_level = [qc_threshold_level, 1 - qc_threshold_level]
            qc_scale_range = norm.ppf(qc_threshold_level, xfrm['scale_mean'], xfrm['scale_std'])
            qc_shift_range = norm.ppf(qc_threshold_level, xfrm['shift_mean'], xfrm['shift_std'])

            scaler_cfg['qc_scale'] = (
                lambda x, range=qc_scale_range: (x >= range[0]) & (x <= range[1]))
            scaler_cfg['qc_shift'] = (
                lambda x, range=qc_shift_range: (x >= range[0]) & (x <= range[1]))

        return keras.models.load_model(model_file), scaler_cfg

    def prepare_loading(self, filename, read_id):
        npread = NanoporeRead(filename, self.fast5prefix, read_id)
        signal_means = npread.load_padded_signal_head(
                self.scaler_cfg['length'], self.scaler_cfg['stride'],
                self.scaler_cfg['min_length'])

        if signal_means is not None:
            self.head_signals.append(signal_means)
            self.head_signal_assoc_reads.append(npread)

        return npread

    def fit_scalers(self):
        if len(self.head_signals) <= 0:
            return

        batch_size = self.config['scaler_maximum_batch_size']
        scfg = self.scaler_cfg

        predmtx = self.scaler_model.predict(
                        np.array(self.head_signals)[:, :, np.newaxis], batch_size)
        scaling_params = np.transpose([
            scfg['xfrm_scale'](predmtx[:, 0]), scfg['xfrm_shift'](predmtx[:, 1])])
        qc_pass_scale = scfg['qc_scale'](scaling_params[:, 0])
        qc_pass_shift = scfg['qc_shift'](scaling_params[:, 1])
        qc_pass = qc_pass_scale & qc_pass_shift

        for npread, paramset, ok in zip(self.head_signal_assoc_reads, scaling_params, qc_pass):
            if ok:
                paramset = np.array(paramset, dtype=scfg['dtype'])
                npread.set_scaling_params(paramset)
            else:
                npread.set_status('scaling_qc_fail', stop=True)


class NanoporeRead:

    fast5 = full_raw_signal = error_message = None
    sequence_length = mean_qscore = num_events = 0
    sequence = scaling_params = label = barcode = polya = None

    def __init__(self, filename, srcdir, read_id):
        self.fullpath = os.path.join(srcdir, filename)
        self.filename = filename
        self.read_id = read_id
        self.status = 'okay'
        self.stopped = False
        self.load()

    def __del__(self):
        self.close()

    def set_status(self, newstatus, stop=False):
        self.status = newstatus
        self.stopped = self.stopped or stop

    def set_error(self, status, error_message):
        self.status = status
        self.error_message = error_message

    def set_scaling_params(self, params):
        self.scaling_params = params

    def set_label(self, newlabel):
        self.label = newlabel

    def set_barcode(self, newbarcode):
        self.barcode = newbarcode

    def set_adapter_trimming_length(self, newlength):
        if self.sequence is None:
            raise Exception('Sequence is not set.')
        self.sequence = self.sequence[:2] + (newlength,)

    def set_polya_tail(self, polya_info):
        self.polya = polya_info

    def is_stopped(self):
        return self.stopped

    def close(self):
        self.full_raw_signal = None
        if self.fast5 is not None:
            self.fast5.close()

    def report(self):
        rep = {'filename': self.filename, 'read_id': self.read_id,
               'status': self.status}

        if self.fast5 is not None:
            rep.update({
                'channel': self.fast5.channel_number,
                'start_time': round(self.fast5.start_time / self.fast5.sampling_rate, 3),
                'run_id': self.fast5.run_id,
                'sample_id': self.fast5.sample_id,
                'duration': self.fast5.duration,
                'num_events': self.num_events,
                'sequence_length': self.sequence_length,
                'mean_qscore': self.mean_qscore,
            })

        if self.sequence is not None:
            rep['sequence'] = self.sequence

        if self.error_message:
            rep['error_message'] = self.error_message

        if self.label is not None:
            rep['label'] = self.label

        if self.barcode is not None:
            rep['barcode'] = self.barcode

        if self.polya is not None:
            rep['polya'] = self.polya

        return rep

    def load(self):
        try:
            fast5 = Fast5Reader(self.fullpath, self.read_id)
        except:
            import traceback
            traceback.print_exc()
            self.set_status('irregular_fast5', stop=True)
            return

        self.fast5 = fast5
        self.sampling_rate = fast5.sampling_rate

    def load_padded_signal_head(self, length_limit, stride, min_length):
        sigload_length = min(length_limit, self.fast5.duration)
        sigload_length = sigload_length - sigload_length % stride

        signal = self.fast5.get_raw_data(end=sigload_length)
        if len(signal) % stride > 0:
            signal = signal[:-(len(signal) % stride)]

        if len(signal) < min_length:
            self.set_status('scaler_signal_too_short', stop=True)
            return

        signal_means = signal.reshape([len(signal) // stride, stride]
                                      ).mean(axis=1, dtype=np.float32)
        length_limit //= stride
        if len(signal_means) < length_limit:
            signal_means = np.pad(signal_means, [length_limit - len(signal_means), 0],
                                  'constant')

        return signal_means

    def load_signal(self, end=None, pool=None, pad=False, scale=True):
        # Load from the cache if available
        if self.full_raw_signal is not None:
            sig = self.full_raw_signal
        else:
            if self.fast5 is None:
                raise Exception('Fast5 must be open for getting signals.')
            sig = self.fast5.get_raw_data(end=end)
            if end is None:
                self.full_raw_signal = sig

        if pool is not None and pool > 1:
            cutend = len(sig) - len(sig) % pool
            sig = sig[:cutend].reshape([len(sig) // pool, pool]
                                       ).mean(axis=1, dtype=np.float32)
        else:
            pool = 1

        if end is not None:
            expected_size = end // pool
            if len(sig) > expected_size:
                sig = sig[-expected_size:]
            elif pad and len(sig) < expected_size:
                sig = np.pad(sig, [expected_size - len(sig), 0], 'constant')

        if scale:
            if self.scaling_params is None:
                raise Exception('Scaling parameters were not set for this read.')

            return np.poly1d(self.scaling_params)(sig)
        else:
            return sig

    def load_fast5_events(self):
        if self.fast5 is None:
            raise Exception('Fast5 must be open for getting events.')

        bcall = self.fast5.get_basecall()
        if bcall is None:
            raise SignalAnalysisError('not_basecalled')

        self.sequence_length = bcall['sequence_length']
        self.mean_qscore = bcall['mean_qscore']
        self.num_events = bcall['num_events']
        self.sequence = bcall['sequence'], bcall['qstring'], 0

        return bcall['events']

    def call_albacore(self, albacore):
        rawdata = self.load_signal(pool=False, scale=False)
        bcall = albacore.basecall(
                    rawdata, self.fast5,
                    os.path.basename(self.filename).rsplit('.', 1)[0])
        if bcall is None:
            raise SignalAnalysisError('not_basecalled')

        self.sequence_length = bcall['sequence_length']
        self.mean_qscore = bcall['mean_qscore']
        self.num_events = bcall['called_events']
        self.sequence = bcall['sequence'], bcall['qstring'], 0

        return bcall['events']

