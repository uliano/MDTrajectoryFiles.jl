# MDTrajectoryFiles

Aim of this package is to provide **low level** and **native** (in the sense of no foreign language dependency) Julia access to trajectory files. At the moment only XTC and DCD files are supported, other formats will eventually be implemented, Amber .nc and Desmond .dtr being next in line.

API unstable and may evolve quickly. For each `struct XXXFile` handle the following methods are provided:

```julia
XXXFile(filename, mode) # constructor
read_atoms(handle::XXXFile, frame, AbstractArray(3, natoms))
read_atoms(handle::XXXFile, frame, AbstractArray(3, natoms), selection::BitVector)
read_box(handle::XXXFile, frame, AbstractArray(XXX dependent))
write_frame(handle::XXXFile, frame, box, atoms)
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

## File Handles

`XtcFile` Feature complete and blazingly fast (reads @ ~2x with respect to C implementation).

`DcdFile` Reads and writes DCD files **without** fixed atoms (support may come in the future), 4th dimension is ignored. Due to the lack of standardization this needs more testing (at the moment only OpenMM and NAMD), you can help providing your datasets.

## Examples

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

### Atom Selection

```julia
julia> spam = XtcFile("/home/uliano/spam.xtc", "r")
XtcFile("/home/uliano/spam.xtc" open for "r" 11 frames with 5558 atoms each)

julia> selection = falses(spam.natoms);

julia> selection[10:20] .=true;

julia> selection[30:40] .=true;

julia> frames = falses(11);

julia> frames[1:2:11] .= true;

julia> coords = Array{Float64}(undef, 3, 22, 6);

julia> k = 1;

julia> for frame in 1:11
           if frames[frame]
               read_atoms(spam, frame, view(coords, :, :, k), selection)
               k += 1
           end
       end

julia> coords
3×22×6 Array{Float64, 3}:
[:, :, 1] =
 19.33  18.72  18.69  19.4   17.69  18.13  18.18  17.14  17.14  16.17  … 
 19.51  21.43  22.0   21.93  21.41  19.35  19.24  18.87  19.24  17.9     
 21.4   20.81  19.86  21.54  21.21  19.76  18.53  20.53  21.46  20.25    

[:, :, 2] =
 22.06  20.08  19.1   20.27  19.97  20.86  21.51  19.73  19.25  19.09  …
 22.98  22.23  22.23  21.16  22.98  24.1   25.11  24.12  23.25  25.25   
 16.96  16.74  17.26  16.5   15.92  18.18  17.97  18.89  18.94  19.53   

;;; … 
```









