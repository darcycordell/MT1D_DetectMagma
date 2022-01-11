# MT1D_DetectMagma

Repository for scripts related to Cordell et al. "Estimating melt fraction in silicic systems using Bayesian inversion of magnetotelluric data" 
(submitted to Journal of Volcanology and Geothermal Research)

Note that this is an active set of scripts which I update regularly so if you are looking for the JVGR-specific codes, see the "JVGR" release which is a time capsule of the scripts as they were when the JVGR article was published.

The primary scripts to run are:

Run_MCMC_JVGR.m: Code which reads model txt files, creates synthetic data, and runs MCMC inversion

Within the JVGR release is "Figures_for_Paper.m" which is a script which reproduces the JVGR figures. This is no longer in the main repository.

To run the scripts, open them in MATLAB and run each code block. See comments in scripts for more details.

Other underlying functions are also included but these do not need to be edited or opened by the user.

The scripts are under the very permissive MIT license. See LICENSE.txt file for details.
