# MT1D_DetectMagma
Repository for scripts related to Cordell et al. "On the detectability of melt-rich lenses in magmatic reservoirs using magnetotellurics" (submitted to Geophysical Research Letters December 2020)

MATLAB scripts required to re-create synthetic data and inversions for
1-D layered-magma bodies.

The primary scripts to run are:

MT1D_DetectMagma_VariableThick.m
MT1D_DetectMagma_VariableMelt.m

To run the scripts, open them in MATLAB and run each code block. See comments in scripts for more details.

Other underlying functions are also included but these do not need to be edited or opened by the user.
	Functions: calc_fwd_1d, logerrorbar, MAL, MAL_solve_phi, manual_legend, and OCCAM1d_simple
