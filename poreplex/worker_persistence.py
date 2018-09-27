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

import sys
import pandas as pd
import numpy as np
import importlib.util as imputil
from pomegranate import (
    HiddenMarkovModel, GeneralMixtureModel, State, NormalDistribution)


__all__ = ['WorkerPersistenceStorage']


class WorkerPersistenceStorage:

    STORAGE_NAME = '__poreplex_persistence'
    MODCACHE_SPACE = sys.modules
    VARIABLES = [
        'segmodel', 'unsplitmodel', 'kmermodel', 'kmersize', 'loader',
        'demuxer', 'albacore', 'polyaanalyzer']

    def __init__(self, config):
        self.config = config

    def retrieve_objects(self, target):
        if self.STORAGE_NAME not in self.MODCACHE_SPACE:
            storage = self.init_persistence_objects(self.config)
        else:
            storage = self.MODCACHE_SPACE[self.STORAGE_NAME].storage

        for varname in self.VARIABLES:
            if varname in storage:
                setattr(target, varname, storage[varname])

        for varname in ['loader', 'demuxer']:
            if varname in storage:
                storage[varname].clear()

    def init_persistence_objects(self, config):
        storage = {
            'segmodel': load_segmentation_model(config['segmentation_model']),
            'unsplitmodel': load_segmentation_model(config['unsplit_read_detection_model']),
            'kmermodel': pd.read_table(config['kmer_model'], header=0, index_col=0),
        }
        storage['kmersize'] = len(storage['kmermodel'].index[0])

        if config['barcoding']:
            from .barcoding import BarcodeDemultiplexer
            storage['demuxer'] = BarcodeDemultiplexer(config['demultiplexing'])

        if config['measure_polya']:
            from .polya import PolyASignalAnalyzer
            storage['polyaanalyzer'] = PolyASignalAnalyzer(config['polya_dwell'])

        if config['albacore_onthefly']:
            from .basecall_albacore import AlbacoreBroker
            storage['albacore'] = AlbacoreBroker(config['albacore_configuration'],
                                                 storage['kmersize'])

        from .signal_loader import SignalLoader
        storage['loader'] = SignalLoader(config['signal_processing'], config['inputdir'])

        fakespec = imputil.spec_from_file_location('spam', 'egg.py')
        persmod = imputil.module_from_spec(fakespec)
        self.MODCACHE_SPACE[self.STORAGE_NAME] = persmod
        persmod.storage = storage

        return storage


# Internal serialization implementation pomegranate to json does not accurately
# recover the original. Use a custom format here.
def load_segmentation_model(modeldata):
    model = HiddenMarkovModel('model')

    states = {}
    for s in modeldata:
        if len(s['emission']) == 1:
            emission = NormalDistribution(*s['emission'][0][:2])
        else:
            weights = np.array([w for _, _, w in s['emission']])
            dists = [NormalDistribution(mu, sigma)
                     for mu, sigma, _ in s['emission']]
            emission = GeneralMixtureModel(dists, weights=weights)
        state = State(emission, name=s['name'])

        states[s['name']] = state
        model.add_state(state)
        if 'start_prob' in s:
            model.add_transition(model.start, state, s['start_prob'])

    for s in modeldata:
        current = states[s['name']]
        for nextstate, prob in s['transition']:
            model.add_transition(current, states[nextstate], prob)

    model.bake()

    return model

