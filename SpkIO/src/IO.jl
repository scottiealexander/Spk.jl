
export load_ts_file, load_wf_file, spk_read_channel, spk_read_header,
    spk_channel_labels, spk_read_par

import JSON

# ============================================================================ #
# SpkHeader type for internal use only
# ============================================================================ #
mutable struct SpkHeader
    nchan::UInt16
    offset::Vector{UInt64}
    file::AbstractString
end
# ============================================================================ #
function load_ts_file(ifile::AbstractString)
    ext = splitext(ifile)[2]
    if ext == ".txt"
        return float(open(readlines, ifile))
    elseif ext != ".ts"
        warn("Input file does not have the proper extension, results may vary")
    end

    open(ifile, "r") do io
        seekend(io)
        nbyte = position(io)
        nel = floor(Int64, nbyte / sizeof(Float64))
        seekstart(io)

        return read(io, Float64, nel)
    end
end
# ============================================================================ #
function load_wf_file(ifile::AbstractString)
    open(ifile, "r") do io
        row = read(io, UInt64)
        col = read(io, UInt64)
        x = read(io, Float64, row, col)

        return x
    end
end
# ============================================================================ #
function spk_read_channel(ifile::AbstractString, label::AbstractString)
    hdr = spk_read_header(ifile)
    labels = spk_channel_labels(hdr)

    label in labels || error("Given channel does not exist")

    kchan = findfirst(isequal(label), labels)

    open(ifile, "r") do io

        seek(io, hdr.offset[kchan] + 64)

        nts = read(io, UInt64)
        wf_pts = read(io, UInt64)
        nmrk = read(io, UInt64)

        data = Dict{Symbol, Any}()

        data[:ts] = Vector{Float64}(undef, nts)
        read!(io, data[:ts])

        if wf_pts > 0
            data[:wf] = Matrix{Float64}(undef, wf_pts, nts)
            read!(io, data[:wf])
        else
            data[:wf] = Matrix{Float64}(undef, 0, 0)
        end
        data[:mrk] = Vector{Float64}(undef, nmrk)
        read!(io, data[:mrk])

        return data
    end
end
# ============================================================================ #
function spk_read_header(ifile::AbstractString)
    isfile(ifile) || error("File \"$(ifile)\" does not exist")

    open(ifile, "r") do io
        nchan = read(io, UInt16)
        offset = Vector{UInt64}(undef, Int64(nchan))
        read!(io, offset)

        return SpkHeader(nchan, offset, ifile)
    end
end
# ============================================================================ #
function spk_channel_labels(ifile::AbstractString)
    return spk_channel_labels(spk_read_header(ifile))
end
# ---------------------------------------------------------------------------- #
function spk_channel_labels(hdr::SpkHeader)
    open(hdr.file, "r") do io
        labels = Vector{AbstractString}(undef, hdr.nchan)
        for k = 1:hdr.nchan
            seek(io, hdr.offset[k])
            tmp = Vector{UInt8}(undef, 64)
            read!(io, tmp)
            labels[k] = strip(String(tmp), '\0')
        end

        return labels
    end
end
# ============================================================================ #
function spk_read_par(ifile::String)
    file, ext = splitext(ifile)
    parfile = file * ".json"
    if !isfile(parfile)
        error("Failed to locate PAR/JSON file: \"$(parfile)\"")
    end
    return JSON.parse(read(parfile, String))
end
# ============================================================================ #
