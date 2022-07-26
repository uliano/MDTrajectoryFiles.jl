# MDTrajectoryFiles

Aim of this package is to provide **low level** and **native** (in the sense of no foreign language dependency) Julia access to trajectory files. At the moment only XTC and DCD files are supported, other formats will eventually be implemented, Amber .nc and Desmond .dtr being next in line.

API is really unstable and may evolve quickly. For each `struct XXXFile` handle the following methods are provided:

```julia
XXXFile(filename, mode) # constructor
read_XXX_atoms(handle::XXXFile, frame, AbstractArray(3, natoms))
read_XXX_box(handle::XXXFile, frame, AbstractArray(XXX dependent))
write_XXX_frame(handle::XXXFile, frame, box, atoms)
```

Also, among different `struct XXXFile`s, a few common proprieties are provided:

```julia
XXXFile.natoms::Int
XXXFile.nframes::Int
XXXFile.time::Vector{Float32}(nframes)
```

## Units

- time: ps
- distance: Å

# XTC Files

Feature complete and blazingly fast (reads @ ~2x with respect to C implementation). 

### Coordinates

```julia
julia> using MDTrajectoryFiles

julia> spam=XtcFile("/home/uliano/spam.xtc","r")
XtcFile("/home/uliano/spam.xtc" open for "r" 11 frames with 5558 atoms each)

julia> coords=Array{Float32}(undef, 3, spam.natoms, spam.nframes);

julia> for frame in 1:spam.nframes
           read_xtc_atoms(spam, frame, view(coords, :, :, frame))
       end

julia> coords
3×5558×11 Array{Float32, 3}:
[:, :, 1] =
 22.37  22.21  23.12  22.78  21.19  …  16.76  16.45  17.44  17.73  16.79
 18.9   18.33  18.22  19.91  18.91     14.77  13.73   3.44   3.34   2.74
 18.93  17.99  19.39  18.73  19.76     29.57  30.63  33.15  32.24  33.27

[:, :, 2] =
 18.01  17.4   17.39  18.68  18.73  …  21.43  20.14  22.95  22.85  23.78
 19.94  20.44  19.49  20.72  18.77     11.64  10.88  33.61  32.99  33.37
 12.97  13.74  12.16  12.55  13.62      5.17   5.39  35.77  36.49  35.36

;;; … 
```

### Time vector

```julia
julia> spam.time
11-element Vector{Float32}:
    0.0
  100.0
  200.0
  300.0
  400.0
  500.0
  600.0
  700.0
  800.0
  900.0
 1000.0

julia> 

```

### Unit Cell

```julia
julia> periodicbox = Matrix{Float32}(undef, 3, 3);

julia> read_xtc_box(spam, 1, periodicbox)

julia> periodicbox
3×3 Matrix{Float32}:
 39.0   0.0   0.0
  0.0  39.0   0.0
  0.0   0.0  39.0

julia> 
```

# DCD Files

Reads and writes DCD files **without** fixed atoms (support may come in the future), 4th dimension is ignored. Due to the lack of standardization this needs more testing (at the moment only OpenMM and NAMD), you can help providing your datasets.

## Coordinates

```julia
julia> spam=DcdFile("/home/uliano/spam.dcd","r")
DcdFile("/home/uliano/spam.dcd" open for "r" 10 frames with 5558 atoms each)

julia> coords=Array{Float32}(undef, 3, spam.natoms, spam.nframes);

julia> for frame in 1:spam.nframes
           read_dcd_atoms(spam, frame, view(coords, :, :, frame))
       end

julia> coords
3×5558×10 Array{Float32, 3}:
[:, :, 1] =
 13.9547  14.6422  13.1881  14.6035  …   0.264754  19.7127    4.90919
 18.9011  18.4144  19.3824  19.6031      6.32946   23.6911   21.4199
 17.587   18.3115  18.2311  17.0208     13.831      8.78946  33.1622

[:, :, 2] =
 11.456   11.866   11.0812  12.3559  …   3.52806  20.7359   3.04038
 18.0378  17.0576  18.5171  18.6206      2.65504  15.5784  19.7867
 13.8721  14.197   14.8017  13.581      10.7042   13.2966  36.9427
```

## Time vector

```julia
julia> spam.time
10-element Vector{Float32}:
  300.0
  400.00003
  500.00003
  600.0
  700.00006
  800.00006
  900.00006
 1000.00006
 1100.0
 1200.0

julia> 
```

## Unit Cell

Due to the plethora of different implementation ([see here](https://github.com/MDAnalysis/mdanalysis/issues/187)), at the moment **no effort** is done to interpret these data but they are provided *as is*

```julia
julia> box = Vector(undef, 6);

julia> box = Vector{Float64}(undef, 6);

julia> read_dcd_box(spam, 1, box)

julia> box
6-element Vector{Float64}:
 38.25117249723212
  0.0
 38.25117249723212
  0.0
  0.0
 38.25117249723212

julia>
```










