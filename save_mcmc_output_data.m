function save_mcmc_output_data(models,rms,MTdata)

    save([MTdata.name,'_run01_MCMCoutput'],'models','rms');

end