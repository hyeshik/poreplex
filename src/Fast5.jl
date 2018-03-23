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

using HDF5

function read_fast5_attributes(file::HDF5File)
    channel_attrs = attrs(file["UniqueGlobalKey/channel_id"])

    (read(channel_attrs["channel_number"]),
     read(channel_attrs["range"]),
     read(channel_attrs["sampling_rate"]),
     read(channel_attrs["digitisation"]),
     read(channel_attrs["offset"]))
end

mutable struct Fast5Reader
    filename::String
    h5::HDF5File

    # Signal sampling-related attributes
    channel_number::String
    signal_range::Float64
    sampling_rate::Float64
    digitization::Float64
    offset::Float64

    function Fast5Reader(filename)
        h5 = h5open(filename, "r")
        finalizer(h5, close)

        channel_number, range, sampling_rate, digitization, offset =
            read_fast5_attributes(h5)

        return new(filename, h5, channel_number, range, sampling_rate,
                   digitization, offset)
    end
end

function load_raw_signals(reader::Fast5Reader, readnumber::Int64,
                          duration::Float64)::Array{Int16}
    sigdataset = reader.h5["Raw/Reads/Read_$readnumber/Signal"]
    sampled_length = length(sigdataset)
    requested_length = Int(round(reader.sampling_rate * duration))
    
    ret_length = min(sampled_length, requested_length)
    sigdataset[1:ret_length]
end

function load_scaled_raw_signals(reader::Fast5Reader, readnumber::Int64,
                                 duration::Float64)::Array{Float64}
    (load_raw_signals(reader, readnumber, duration) + reader.offset) * (
        reader.signal_range / reader.digitization)
end
