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

__all__ = ['check_albacore', 'AlbacoreBroker']

# All albacore modules are loaded from the insides of functions to enable
# graceful detection of ImportErrors in check_albacore().
from itertools import product
from collections import OrderedDict
import pandas as pd
import os

REQUIRED_ALBACORE_VERSION = '2.3.0'

def prepare_albacore(configpath, flowcell, kit, version_min=REQUIRED_ALBACORE_VERSION):
    from distutils.version import LooseVersion
    try:
        from albacore import __version__
    except ImportError:
        raise Exception("Could not load `albacore.' Check PYTHONPATH.")

    if LooseVersion(__version__) < LooseVersion(version_min):
        raise Exception('Albacore {} is detected although {} is required.'.format(
            __version__, version_min))

    # Generate configuration file from template
    from albacore.path_utils import get_default_path
    from albacore.config_selector import choose_config
    from configparser import ConfigParser

    try:
        albacore_data_path = get_default_path('.')
        template = choose_config(albacore_data_path, flowcell, kit)
    except:
        raise Exception('Unsupported flowcell or kit from albacore.')

    config = ConfigParser(interpolation=None)
    config.read_file(open(template[0]))
    config['basecaller']['model_path'] = albacore_data_path
    config['calib_detector']['method'] = ''
    config.write(open(configpath, 'w'))

    return __version__


class AlbacoreBroker:

    LAYOUT_FILE = 'layout_raw_basecall_1d.jsn'
    WORKERS_SINGLE_PROCESS = 0

    def __init__(self, configpath, kmer_size):
        from albacore.pipeline_core import PipelineCore
        from albacore.path_utils import get_default_path

        self.albacore_data_path = get_default_path('.')
        self.descr_file = os.path.join(self.albacore_data_path, self.LAYOUT_FILE)
        self.core = PipelineCore(self.descr_file, configpath, self.WORKERS_SINGLE_PROCESS)

        self.kmer_decode_table = list(map(''.join, product('ACGT', repeat=kmer_size)))

    def adopt_basecalled_table(self, result):
        field_names = 'mean start stdv length model_state move p_model_state weights'.split()
        return pd.DataFrame.from_dict(OrderedDict(
            [(field, result[field] if field != 'model_state'
              else list(map(self.kmer_decode_table.__getitem__, result[field])))
             for field in field_names]))

    def basecall(self, rawdata, fast5, filename):
        read_id = fast5.read_id

        self.core.pass_data({
            'sampling_rate': fast5.sampling_rate,
            'start_time': fast5.start_time,
            'read_id': read_id,
            'data_id': filename,
            'raw': rawdata,
        })
        try: self.core.finish_all_jobs()
        except RuntimeError: pass
        results = self.core.get_results()

        if read_id not in results or 'basecall_1d_callback' not in results[read_id]:
            return None

        bc1dresult = results[read_id]['basecall_1d_callback']
        return {
            'events': self.adopt_basecalled_table(bc1dresult),
            'sequence': bc1dresult['sequence'][::-1].replace('T', 'U'),
            'qstring': bc1dresult['qstring'][::-1],
            'sequence_length': bc1dresult['sequence_length'],
            'mean_qscore': round(bc1dresult['mean_qscore'], 3),
            'called_events': bc1dresult['called_events'],
        }


if __name__ == '__main__':
    import sys
    try:
        version = prepare_albacore(*sys.argv[1:])
    except Exception as exc:
        print(str(exc))
    else:
        print('okay', version)
