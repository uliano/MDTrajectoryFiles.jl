module MDTrajectoryFiles

# Write your package code here.

include("xtcfiles.jl")
include("dcdfiles.jl")

export XtcFile, DcdFile, write_frame, read_atoms, read_box

end
