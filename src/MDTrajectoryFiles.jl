module MDTrajectoryFiles

# Write your package code here.

include("xtcfiles.jl")
export XtcFile, write_xtc_frame, read_xtc_atoms, read_xtc_box

include("dcdfiles.jl")
export DcdFile, write_dcd_frame, read_dcd_atoms, read_dcd_box

end
