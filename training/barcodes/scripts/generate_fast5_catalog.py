#!/usr/bin/env python3
#
# Copyright (c) 2017 Institute for Basic Science
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

import pandas as pd
import numpy as np
import feather
import os
import sys
import glob

def load_catalog(datadir, runname):
    summary = pd.read_table(os.path.join(datadir, 'sequencing_summary.txt'),
                            low_memory=False)

    summary['run_name'] = runname
    return summary

def process_all(allrundirs):
    allcats = []

    for rundir in allrundirs:
        runname = rundir.split('/')[0]
        print('Loading {} ... '.format(runname), end='')
        sys.stdout.flush()
        allcats.append(load_catalog(rundir, runname))
        print('{:,}'.format(len(allcats[-1])))

    print('Loading finished. Saving merged table... ', end='')
    sys.stdout.flush()
    allcats = pd.concat(allcats, axis=0).sort_values(
        by=['run_name', 'start_time', 'channel']).reset_index(drop=True)
    print('({:,} reads total)'.format(len(allcats)))

    feather.write_dataframe(allcats, 'sequencing_summary.feather')

    print('Done.')

if __name__ == '__main__':
    input_dirs = sys.argv[1:]
    process_all(input_dirs)
