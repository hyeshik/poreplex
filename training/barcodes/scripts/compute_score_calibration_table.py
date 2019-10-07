#!/usr/bin/env python3
#
# Copyright (c) 2019 Seoul National University
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

from rpy2 import robjects as ro
from scipy.interpolate import interp1d
from scipy.optimize import bisect
import pandas as pd
import numpy as np


def loess_fit(x, y, px=None, model=None, alpha=0.5):
    if model is None:
        model = ro.r('y ~ x')

    if px is None:
        px = np.linspace(min(x), max(y), 22)[1:-1]

    fitframe = ro.DataFrame({'x': ro.FloatVector(x), 'y': ro.FloatVector(y)})
    loessmodel = ro.r.loess(model, fitframe, span=alpha)

    predframe = ro.DataFrame({'x': ro.FloatVector(px)})
    predy = ro.r.predict(loessmodel, predframe)
    preddata = [(x, predy[i]) for i, x in enumerate(px)]

    return np.array(preddata).transpose()


class ScoreCalibrationTableGenerator:

    SCORING_BINNING_PARAMS = [ # [window_size, minimum_size, interval]
        [10000, 2500, 3300],
        [2000, 1000, 1000],
        [1000, 500, 500],
    ]
    SCORING_STDEV_THRESHOLD = 0.02
    EXTRAPOLATION_SUPPORT_POINTS = 3
    INTERPOLATION_LOESS_ALPHA = 0.3

    def __init__(self, predtbl):
        self.preds = pd.read_table(predtbl).sort_values(by='pred_prob', ascending=False)

        # Remove predictions from/to decoys
        self.preds = self.preds[(self.preds['prior_label'] > 0) &
                                (self.preds['pred_label'] > 0)].copy()
        self.preds['correct'] = (self.preds['prior_label'] == self.preds['pred_label'])

    def calc_error_rate(self, start, end):
        rows = self.preds[start:end]
        rows = rows[rows['pred_label'] != 0]
        return (~rows['correct']).mean()

    def scan_error_rates(self, window_size, min_width, interval):
        totalpoints = len(self.preds)
        windows = []

        # Prepare binning ranges to count errors
        for start in range(0, totalpoints - min_width + 1, interval):
            end = min(totalpoints, start + window_size)
            start = max(0, start)
            width = end - start

            windows.append([start, end, width])

        # Calculate error rates
        erate_bins = pd.DataFrame(windows, columns=['start', 'end', 'width'])
        erate_bins['error_rate'] = erate_bins.apply(lambda row:
            self.calc_error_rate(int(row['start']), int(row['end'])),
            axis=1)

        predscoreloc = self.preds['pred_prob'].iloc

        # Calculate score metrics
        erate_bins['score_high'] = erate_bins['start'].apply(predscoreloc.__getitem__)
        erate_bins['score_low'] = (erate_bins['end']-1).apply(predscoreloc.__getitem__)
        erate_bins['score_mean'] = erate_bins.apply(lambda x:
            predscoreloc[int(x['start']):int(x['end'])].mean(), axis=1)
        erate_bins['score_std'] = erate_bins.apply(lambda x:
            predscoreloc[int(x['start']):int(x['end'])].std(), axis=1)

        return erate_bins

    def build_multiscale_error_table(self):
        errortbl = None

        for binparams in self.SCORING_BINNING_PARAMS:
            errorstat = self.scan_error_rates(*binparams)
            if errortbl is None:
                errortbl = errorstat
            else:
                prevtbl = errortbl[errortbl['score_std'] < self.SCORING_STDEV_THRESHOLD]
                prevminscore = prevtbl['score_mean'].min()
                newtbl = errorstat[errorstat['score_mean'] < prevminscore]
                errortbl = pd.concat([prevtbl, newtbl], axis=0)

        return errortbl

    def compute_low_extrapolated_table(self, errortbl):
        probscan_lowlimit = errortbl.iloc[-1]['score_mean']

        # Extrapolate predicted score -> error rate mappings using internal supports
        extrapol_supports = errortbl.iloc[-self.EXTRAPOLATION_SUPPORT_POINTS:]
        coeffs = np.polyfit(extrapol_supports['score_mean'],
                            extrapol_supports['error_rate'], 1)
        extrapol_range = np.array([0, probscan_lowlimit])
        extrapol_y = np.poly1d(coeffs)(extrapol_range)
        extrapol_inv = np.poly1d([1/coeffs[0], - coeffs[1]/coeffs[0]])

        extrapol_phred_max = int(-np.log10(probscan_lowlimit) * 10)
        extrapol_phred_scores = np.arange(extrapol_phred_max + 1)
        extrapol_score_pred = extrapol_inv(10 ** (-extrapol_phred_scores / 10))

        ret = []
        print('Extrapolated map for score range [0, {:.4})'.format(probscan_lowlimit))
        for scorepred, phred in zip(extrapol_score_pred, extrapol_phred_scores):
            print('{:2}  {:.6}'.format(phred, scorepred))
            ret.append([phred, scorepred])

        # Fix table so that scores in range are non-negative.
        assert ret[0][0] == 0
        ret[0][1] = 0.0

        return ret

    def compute_interpolated_table(self, errortbl, extrapol_phred_max):
        probscan_lowlimit = errortbl.iloc[-1]['score_mean']
        probscan_hilimit = errortbl.iloc[0]['score_mean']
        interpol_phred_max = int(-np.log10(errortbl.iloc[0]['error_rate']) * 10)

        print('\nInterpolated map for score range [{:.4}, {:.6}]'.format(
            probscan_lowlimit, probscan_hilimit))

        # Interpolate
        pred_x = np.array(sorted(np.hstack([
            np.linspace(probscan_lowlimit, probscan_hilimit, 100),
            errortbl['score_mean']])))
        pred_x, pred_y = loess_fit(errortbl['score_mean'], errortbl['error_rate'],
                                   pred_x, alpha=self.INTERPOLATION_LOESS_ALPHA)

        ret = []
        for phred in range(extrapol_phred_max + 1, interpol_phred_max + 1):
            exp_error_rate = 10 ** (-phred / 10)
            ifun = interp1d(pred_x, pred_y - exp_error_rate)
            scorepred = bisect(ifun, probscan_lowlimit, probscan_hilimit)

            print('{:2}  {:.6}'.format(phred, scorepred))
            ret.append([phred, scorepred])
        return ret

    def build_calibration_table(self):
        errortbl = self.build_multiscale_error_table()
        lowtbl = self.compute_low_extrapolated_table(errortbl)
        hightbl = self.compute_interpolated_table(errortbl, lowtbl[-1][0])
        return pd.DataFrame(lowtbl + hightbl, columns=['phred', 'pred_score'])


if __name__ == '__main__':
    import sys
    score_pred_file, output_file = sys.argv[1:]
    calgen = ScoreCalibrationTableGenerator(score_pred_file)
    caltbl = calgen.build_calibration_table()
    caltbl.to_csv(output_file, sep='\t', index=False)

