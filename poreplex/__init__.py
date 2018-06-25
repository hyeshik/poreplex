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

__all__ = [
    '__version__',
    'OUTPUT_NAME_PASSED', 'OUTPUT_NAME_FAILED',
    'OUTPUT_NAME_ARTIFACT', 'OUTPUT_NAME_BARCODES',
]

__version__ = '0.1'

OUTPUT_NAME_PASSED = 'pass'
OUTPUT_NAME_FAILED = 'fail'
OUTPUT_NAME_ARTIFACT = 'artifact'
OUTPUT_NAME_BARCODES = 'BC{n}'

