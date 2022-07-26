import Base: show
using StaticArrays
using Dates

"""
    DcdFile

Handle to a dcd File. Opening in read mode reads the header 
and stores the offsets to the frames for fast random access.     

# Useful fields, all initialized on opening in read mode.

- natoms::Int
- nframes::Int
- time::Vector{Float32}(nframes)

# Constructor

    DcdFile(name::AbstractString, mode::AbstractString)

mode is the usual string r,w,a
"""
mutable struct DcdFile
    file::IOStream
    filename::AbstractString
    mode::AbstractString
    natoms::Int
    nframes::Int
    startstep::Int
    stepsperframe::Int
    offsets::Vector{UInt64}
    time::Vector{Float32}
    nfixedatoms::Int
    title::Vector{String}
    dt::Float64
    has64bitblocksize::Bool
    has4d::Bool
    hasbox::Bool
    swapendian::Bool
    ischarmm::Bool

    function DcdFile(name::AbstractString, mode::AbstractString; kwargs...)
        if 'w' in mode            
            # ensure that title has at least 2 lines and that
            # lines have at most 79 charachters.
            if ! haskey(kwargs, :title)
                title = String["MDTrajectoryFiles.jl",
                                "Created $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))."]
            elseif typeof(kwargs[:title] == String)
                title = String[kwargs[:title][1:79], 
                                "Created $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))."]
            elseif typeof(kwargs[:title] == Vector{String})
                title = [str[1:79] for str in kwargs[:title]]
            else
                error("title should be either String or Vector{String}.")
            end
            if haskey(kwargs, :natoms)
                natoms = Int(kwargs[:natoms])
                if natoms <= 0
                    error("natoms must be > 0")
                end
            else
                error("missing natoms parameter")
            end
            if haskey(kwargs, :startstep)
                startstep = Int(kwargs[:startstep])
                if startstep <= 0
                    error("startstep must be > 0")
                end
            else
                error("missing startstep parameter")
            end
            if haskey(kwargs, :stepsperframe)
                stepsperframe = Int(kwargs[:stepsperframe])
                if stepsperframe <= 0
                    error("stepsperframe must be > 0")
                end
            else
                error("missing stepsperframe parameter")
            end
            if haskey(kwargs, :dt)
                dt = Float64(kwargs[:dt])
                if dt <= 0.0
                    error("dt must be > 0")
                end
            else
                error("missing dt parameter")
            end 
            file = open(name, "w+")
            write_dcd_header(file, natoms, 0, startstep, stepsperframe, dt, title)
            DcdF = new(file, name, mode, natoms, 0, startstep, stepsperframe, 
                       [], [], 0, title, dt, false,  false, true, false, true)
        elseif 'r' in mode || 'a' in mode
            file = open(name, "a+")
            stuff = read_dcd_header(file)
            DcdF = new(file, name, mode, stuff.natoms, stuff.nframes, 
                       stuff.startstep, stuff.stepsperframe, stuff.offsets, stuff.time,
                       stuff.nfixedatoms, stuff.title, stuff.deltat, stuff.has64bitblocksize,
                       stuff.has4d, stuff.hasbox, stuff.swapendian, stuff.ischarmm)
        else
            error("mode should be 'r', 'w' or 'a'.")
        end
        closeme(dcdf) = close(dcdf.file)
        finalizer(closeme, DcdF)
    end
end

function show(io::IO, dcdf::DcdFile)
    result = "DcdFile(\"$(dcdf.filename)\" open for \"$(dcdf.mode)\" "
    result *= "$(dcdf.nframes) frames with $(dcdf.natoms) atoms each)"
    print(io, result)
end

function readblocksize(file, blocksizetype, swapendian)
    blocksize = read(file, blocksizetype)
    if swapendian
        blocksize = bswap(blocksize)
    end
    return blocksize
end

function write_dcd_header(file, natoms, nframes, startstep, stepsperframe, dt, title)
    titlebuf = zeros(UInt8, 80, length(title)) 
    for (j, str) in enumerate(title)
        for (i, c) in enumerate(str)
            titlebuf[i, j] = UInt8(c)
        end
    end
    blocksize = Int32[84]
    inctrl = zeros(Int32, 20)
    deltat = Float32(dt / 48.88821f-3 ) # time in AKMA units
    inctrl[1] = Int32(nframes)
    inctrl[2] = Int32(startstep)
    inctrl[3] = Int32(stepsperframe)
    inctrl[4] = inctrl[1] * inctrl[3]
    inctrl[10] = reinterpret(Int32, deltat)
    inctrl[11] = one(Int32)
    inctrl[20] = Int32(24) # pretend to be charmm 24
    magic = "CORD"
    write(file, blocksize)
    write(file, magic)
    write(file, inctrl)
    write(file, blocksize)
    blocksize[1] = Int32(length(titlebuf) + 4)
    write(file, blocksize)
    write(file, Int32(size(titlebuf, 2)))
    write(file, titlebuf)
    write(file, blocksize)
    blocksize[1] = Int32(4)
    write(file, blocksize)
    write(file, Int32(natoms))
    write(file, blocksize)
end

function read_dcd_header(file)
    has4d = false
    ischarmm = false
    has64bitblocksize = false
    blocksizetype = Int32
    swapendian = false
    seekend(file)
    filesize = position(file)
    seekstart(file)
    # Guess endianness and size of block read_xtc_headers
    blocksize64 = read(file, Int64)
    tmp32 = reinterpret(reshape, Int32, [blocksize64])
    blocksize32 = tmp32[1]
    blocksize = blocksize32
    magic = String(vec(reinterpret(reshape, UInt8, [tmp32[2]])))
    if blocksize64 == 84 || blocksize64 == 6052837899185946624 
        blocksizetype = Int64
        has64bitblocksize = true
        if blocksize == 6052837899185946624
            swapendian = true
        end
        tmp8 = Vector{UInt8}(undef, 4)
        read!(file, tmp8)
        magic = String(tmp8)
        blocksize = blocksize64
    elseif blocksize32 == 1409286144
        swapendian = true
    elseif blocksize32 != 84
        error("Wrong first block size in dcd file.")
    end
    if swapendian
        blocksize = bswap(blocksize)
    end
    if magic == "VELD"
        error("Velocity file not supported, yet.")
    elseif magic != "CORD"
        error("Unknown header magic string.")
    end
    # finish reading first block
    inctrl = Vector{Int32}(undef, 20)
    read!(file, inctrl)
    doublets = Vector{Int32}(undef,2)
    doublets[1] = inctrl[10]
    doublets[2] = inctrl[11]
    if swapendian
        inctrl .= bswap.(inctrl)
    end
    nframes = Int(inctrl[1])
    startstep = Int(inctrl[2])
    stepsperframe = Int(inctrl[3])
    nfixedatoms = Int(inctrl[9])
    if nfixedatoms > 0 
        error("Fixed atoms not supported, yet")
    end
    ischarmm = inctrl[20] == 0 ? false : true
    if ischarmm
        deltat = reinterpret(Float32, inctrl[10])
        hasbox = inctrl[11] > 0
        has4d = inctrl[12] > 0
        if inctrl[13] > 0
            error("Fluctuating charges not supported.")
        end
    else
        deltat64 = reinterpret(reshape, Float64, doublets)
        if swapendian
            deltat64 = bswap(deltat64)
        end
        deltat = Float32(deltat64)
    end
    deltat *= 48.88821f-3 # we use picoseconds, not AKMA units
    if blocksize != readblocksize(file, blocksizetype, swapendian) # this closes first block
        error("Block size mismatch.")
    end
    blocksize = readblocksize(file, blocksizetype, swapendian) # open second block
    ntitle = read(file, Int32)
    if swapendian
        ntitle = bswap(ntitle)
    end
    if blocksize != 4 + 80 * ntitle
        error("Block size should be 4 + 80 * NTITL")
    end  
    title = String[]
    for i in 1:ntitle
        chars = Vector{UInt8}(undef, 80)
        read!(file, chars)
        push!(title, strip(String(chars), '\0'))
    end
    if blocksize !=  readblocksize(file, blocksizetype, swapendian) # this closes second block
        error("Block size mismatch.")
    end  
    blocksize = readblocksize(file, blocksizetype, swapendian) # open third block
    if blocksize != 4 
        error("Third block should be 4 bytes.")
    end
    nat = read(file, Int32)
    if swapendian
        nat = bswap(nat)
    end
    natoms = Int(nat)
    if blocksize !=  readblocksize(file, blocksizetype, swapendian) # this closes third block
        error("Block size mismatch.")
    end  
    datastart = position(file)  # end of header
    framesize = (3 + Int(inctrl[12])) * (sizeof(blocksizetype) * 2 + natoms * 4)
    if hasbox 
        framesize += sizeof(blocksizetype) * 2 + 6 * 8
    end
    if nframes == 0
        nframes = (filesize - datastart) ÷ framesize
    end 
    if nframes * framesize != filesize - datastart || mod(filesize - datastart, framesize) != 0

        error("Wrong file size.")
    end

    offsets = UInt64[datastart + (i - 1) * framesize for i in 1:nframes]

    time = Float64[(startstep + (i - 1) * stepsperframe) * deltat for i in 1:nframes]

    return (; natoms, nframes, startstep, stepsperframe, deltat, offsets, time,nfixedatoms,
            title, has64bitblocksize, has4d, hasbox, swapendian, ischarmm)
end


"""
    read_dcd_atoms(file::DcdFile, frame::Integer, coords::AbstractMatrix{T}) where {T <: Real}

Reads atom coordinates for a frame.
    
# Parameters
- file: must be already opened for "r"
- frame: must be in 1:nframes
- coords: must be preallocated with dimension (3, natoms)
"""
function read_dcd_atoms(file::DcdFile, frame::Integer, coords::AbstractMatrix{T}) where {T <: Real}
    buffer = Vector{Float32}(undef, file.natoms)
    offset = file.offsets[frame]
    blocksizetype = file.has64bitblocksize ? Int64 : Int32
    if file.hasbox
        offset += sizeof(blocksizetype) * 2 + 6 * 8
    end
    seek(file.file, offset)
    for i in 1:3
        blocksize = readblocksize(file.file, blocksizetype, file.swapendian)
        if blocksize != file.natoms * 4
            error("wrong number of coordinates.")
        end
        read!(file.file, buffer)
        if file.swapendian
            coords[i, :] .= T.(bswap.(buffer))
        else
            coords[i, :] .= T.(buffer)
        end
        if blocksize != readblocksize(file.file, blocksizetype, file.swapendian)
            error("Block size mismatch.")
        end
    end
    return nothing
end 


"""
    read_dcd_box(file::DcdFile, frame::Integer, box::AbstractVector{T}) where {T <: Real}

# Parameters
    - file: must be already opened for "r"
    - frame: must be in 1:nframes
    - box: must be preallocated with dimension (6,)
"""
function read_dcd_box(file::DcdFile, frame::Integer, box::AbstractVector{T}) where {T <: Real}
    seek(file.file, file.offsets[frame])
    blocksizetype = file.has64bitblocksize ? Int64 : Int32
    blocksize = readblocksize(file.file, blocksizetype, file.swapendian)
    bsize = 6 * 8
    if blocksize != bsize 
        error("Extrablock size should be $bsize.")
    end
    rawbox = MVector{6, Float64}(undef)
    read!(file.file, rawbox)
    if file.swapendian
        rawbox .= bswap.(rawbox)
    end
    box .= T.(rawbox)
    if blocksize != readblocksize(file.file, blocksizetype, file.swapendian)
        error("Block size mismatch.")
    end
    return nothing
end

"""
    write_dcd_frame(file::DcdFile, box::AbstractVector{S}, 
                    coords::AbstractMatrix{T}) where {S, T <: Real}

# Parameters
    - file: must be already opened for "w" or "a"
    - box dimension (6, ) The interpretation of these
      parameter varies with the program and version
    - coords: must be preallocated with dimension (3, natoms)
"""
function write_dcd_frame(file::DcdFile, box::AbstractVector{S}, coords::AbstractMatrix{T}) where{S, T <: Real}
    buffer = Vector{Float32}(undef, file.natoms)
    boxbuffer = MVector{6, Float64}(undef)
    blocksize = Int32[6 * 8]
    inctrl = Vector{Int32}(undef, 4)
    seekend(file.file)
    write(file.file, blocksize)
    boxbuffer .= Float64.(box)
    write(file.file, boxbuffer)
    write(file.file, blocksize)
    blocksize[1] = Int32(file.natoms * 4)
    for i in 1:3
        write(file.file, blocksize)
        buffer .= Float32.(coords[i, :])
        write(file.file, buffer)
        write(file.file, blocksize)
    end
    seek(file.file, 8)
    #increment nframes
    read!(file.file, inctrl)
    inctrl[1] += 1
    inctrl[4] = inctrl[1] * inctrl[3]
    write(file.file, inctrl)
    return nothing
end

"""
    read_dcd_file(filename) -> Array{Float32}(3, natoms, nframes)

Reads all frames in a xtc trajectory files and returns coordinates.    
"""
function read_dcd_file(name)
    dcdfile = DcdFile(name, "r")
    coordinates = Array{Float32}(undef, 3, dcdfile.natoms, dcdfile.nframes)
    for frame in 1:dcdfile.nframes
        read_dcd_atoms(dcdfile, frame, view(coordinates, :, :, frame))
    end
    return coordinates
end