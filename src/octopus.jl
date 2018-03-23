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

struct Fast5SignalProcessor
    fast5::Fast5Reader
    readnumber::Int64
    peekduration::Float64
    filtersize::Int64
    firstindex::Int64

    function Fast5SignalProcessor(filename, readnumber,
                                  peekduration=3.0::Float64,
                                  filtersize=5::Int64)
        fast5 = Fast5Reader(filename)

        if filtersize ≤ 0 || filtersize % 2 != 1
            throw(ArgumentError("filtersize must be a positive odd number."))
        end

        firstindex::Int64 = (filtersize - 1) / 2
        new(fast5, readnumber, peekduration, filtersize, firstindex)
    end
end


function load_filtered_signal(proc::Fast5SignalProcessor)::Array{Float64}
    raw = load_scaled_raw_signal(proc.fast5, proc.readnumber, proc.peekduration)

    [median(raw[left:(left + proc.filtersize - 1)])
     for left = 1:(length(raw) - proc.filtersize + 1)]
end



# main
filename = "/guava/ourseq/2018/20180309-Nanopore-Mux_v3/20180308_0547_Mux_v3/fast5/pass/0/driller_20180307_FC21_MN17871_mux_scan_Mux_v3_62788_read_99_ch_490_strand.fast5"
proc = Fast5SignalProcessor(filename, 99, 0.1, 5)
sig = load_filtered_signal(proc)
println(sig)

