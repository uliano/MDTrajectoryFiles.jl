using MDTrajectoryFiles
using Test
using Base.Filesystem
BASE_FOLDER = dirname(pathof(MDTrajectoryFiles))

@testset "Test XTC read and write" begin
    # read a file, write it and compare.
    tmp = tempname()
    testfile = joinpath(BASE_FOLDER, "data", "test.xtc")
    xtcfile = XtcFile(testfile, "r")
    xtcfile2 = XtcFile(tmp, "w")
    nframes = xtcfile.nframes
    natoms = xtcfile.natoms
    coordinates = Matrix{Float32}(undef, 3, natoms)
    box = zeros(Float32, 3,3)
    for frame in 1:nframes
        read_xtc_box(xtcfile, frame, box)
        read_xtc_atoms(xtcfile, frame, coordinates)
        write_xtc_frame(xtcfile2, xtcfile.steps[frame], xtcfile.time[frame], box, coordinates)
    end
    close(xtcfile.file)
    close(xtcfile2.file) # this is required as XtcFile closes on destruction
    data = read(testfile)
    data2 = read(tmp)
    @test data == data2
end
