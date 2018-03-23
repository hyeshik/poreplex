#
# Copyright (c) 2018 Hyeshik Chang
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

include("./Fast5.jl")

mutable struct Fast5SignalProcessor
    fast5::Fast5Reader
    filter_size::Int64

    function Fast5SignalProcessor(filename, peek_duration=3.0::Float64,
                                  filter_size=5::Int64)
        fast5 = Fast5Reader(filename)
        new(fast5, filter_size)
    end
end

function load_filtered_signal(proc::Fast5SignalProcessor)
    
end

# main
filename = "/guava/ourseq/2018/20180309-Nanopore-Mux_v3/20180308_0547_Mux_v3/fast5/pass/0/driller_20180307_FC21_MN17871_mux_scan_Mux_v3_62788_read_99_ch_490_strand.fast5"
proc = Fast5SignalProcessor(filename, 5)

sig = load_scaled_raw_signals(f5file, 99, 0.1)
println(typeof(sig))
println(sig)

