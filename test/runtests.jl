using MDTrajectoryFiles
using Test
using Base.Filesystem
BASE_FOLDER = dirname(pathof(MDTrajectoryFiles))

@testset "Test XTC read and write" begin
    # read a file, write it and compare.
    tmp = tempname()
    testfile = joinpath(BASE_FOLDER, "data", "gromacs.xtc")
    xtcfile = XtcFile(testfile, "r")
    xtcfile2 = XtcFile(tmp, "w")
    nframes = xtcfile.nframes
    natoms = xtcfile.natoms
    coordinates = Matrix{Float32}(undef, 3, natoms)
    box = zeros(Float32, 3,3)
    for frame in 1:nframes
        read_box(xtcfile, frame, box)
        read_atoms(xtcfile, frame, coordinates)
        write_frame(xtcfile2, xtcfile.steps[frame], xtcfile.time[frame], box, coordinates)
    end
    close(xtcfile.file)
    close(xtcfile2.file) # this is required as XtcFile closes on destruction
    data = MDTrajectoryFiles.read_xtc_file(testfile)
    data2 = MDTrajectoryFiles.read_xtc_file(tmp)
    @test data == data2
end

@testset "Test XTC read and write, 5 atoms trajectory" begin
    # read a file, write it and compare.
    tmp = tempname()
    testfile = joinpath(BASE_FOLDER, "data", "fiveatoms.xtc")
    xtcfile = XtcFile(testfile, "r")
    xtcfile2 = XtcFile(tmp, "w")
    nframes = xtcfile.nframes
    natoms = xtcfile.natoms
    coordinates = Matrix{Float32}(undef, 3, natoms)
    box = zeros(Float32, 3,3)
    for frame in 1:nframes
        read_box(xtcfile, frame, box)
        read_atoms(xtcfile, frame, coordinates)
        write_frame(xtcfile2, xtcfile.steps[frame], xtcfile.time[frame], box, coordinates)
    end
    close(xtcfile.file)
    close(xtcfile2.file) # this is required as XtcFile closes on destruction
    data = MDTrajectoryFiles.read_xtc_file(testfile)
    data2 = MDTrajectoryFiles.read_xtc_file(tmp)
    @test data == data2
end

@testset "Test DCD read and write" begin
    # read a file, write it and compare.
    tmp = tempname()
    testfile = joinpath(BASE_FOLDER, "data", "openmm.dcd")
    dcdfile = DcdFile(testfile, "r")
    nframes = dcdfile.nframes
    natoms = dcdfile.natoms
    dt = dcdfile.dt
    startstep = dcdfile.startstep
    stepsperframe = dcdfile.stepsperframe
    dcdfile2 = DcdFile(tmp, "w"; natoms=natoms, stepsperframe=stepsperframe,
                       startstep=startstep, dt=dt)
    coordinates = Matrix{Float32}(undef, 3, natoms)
    box = zeros(6)
    for frame in 1:nframes
        read_box(dcdfile, frame, box)
        read_atoms(dcdfile, frame, coordinates)
        write_frame(dcdfile2, box, coordinates)
    end
    close(dcdfile.file)
    close(dcdfile2.file) # this is required as DcdFile closes on destruction
    data = MDTrajectoryFiles.read_dcd_file(testfile)
    data2 = MDTrajectoryFiles.read_dcd_file(tmp)

    @test data == data2
end

@testset "Test XTC read selection" begin
    testfile = joinpath(BASE_FOLDER, "data", "gromacs.xtc")
    xtcfile = XtcFile(testfile, "r")
    selection = falses(xtcfile.natoms)
    selection[10:20] .= true
    selection[30:40] .= true
    frames = falses(xtcfile.nframes)
    frames[2] = true
    frames[4] = true
    coords = Array{Float32}(undef, 3, 22, 2)
    read_atoms(xtcfile, 2, view(coords, :, :, 1), selection)
    read_atoms(xtcfile, 4, view(coords, :, :, 2), selection)
    data = MDTrajectoryFiles.read_xtc_file(testfile)

    @test coords == data[:, selection, frames]
end

@testset "Test XTC read selection 5 atoms" begin
    testfile = joinpath(BASE_FOLDER, "data", "fiveatoms.xtc")
    xtcfile = XtcFile(testfile, "r")
    selection = falses(xtcfile.natoms)
    selection[2] = true
    selection[4] = true
    frames = falses(xtcfile.nframes)
    frames[2] = true
    frames[4] = true
    coords = Array{Float32}(undef, 3, 2, 2)
    read_atoms(xtcfile, 2, view(coords, :, :, 1), selection)
    read_atoms(xtcfile, 4, view(coords, :, :, 2), selection)
    data = MDTrajectoryFiles.read_xtc_file(testfile)

    @test coords == data[:, selection, frames]
end

@testset "Test DCD read selection" begin
    testfile = joinpath(BASE_FOLDER, "data", "openmm.dcd")
    dcdfile = DcdFile(testfile, "r")
    selection = falses(dcdfile.natoms)
    selection[10:20] .= true
    selection[30:40] .= true
    frames = falses(dcdfile.nframes)
    frames[2] = true
    frames[4] = true
    coords = Array{Float32}(undef, 3, 22, 2)
    read_atoms(dcdfile, 2, view(coords, :, :, 1), selection)
    read_atoms(dcdfile, 4, view(coords, :, :, 2), selection)
    data = MDTrajectoryFiles.read_dcd_file(testfile)

    @test coords == data[:, selection, frames]
end
