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

__all__ = ['union_intervals', 'errx', 'errprint', 'ensure_dir_exists']

import sys
import os

def union_intervals(iset):
    merged = []

    for begin, end in sorted(iset):
        if merged:
            if merged[-1][-1] >= begin:
                if merged[-1][-1] < end:
                    merged[-1][-1] = end
                continue
        merged.append([begin, end])

    return merged

def errx(msg):
    errprint(msg)
    sys.exit(254)

def errprint(msg, end=None):
    #import traceback, sys
    #traceback.print_stack(file=sys.stderr)
    print(msg, file=sys.stderr, end=end)

def ensure_dir_exists(filepath):
    dirname = os.path.dirname(filepath)
    if not os.path.isdir(dirname):
        try:
            os.makedirs(dirname)
        except FileExistsError:
            pass
