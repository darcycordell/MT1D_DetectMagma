%MCMC JVGR

%% LOAD TXT FILE MODEL AND CREATE SYNTHETIC DATA
close all; clearvars;

FileName = 'geometry_3_err_005.txt';

    [~,name,ext] = fileparts(FileName);
    OutputFileName = [name,'.mat'];

    synthetic = 1;
    fid = fopen(FileName);
    val = fscanf(fid,'%f %f',[2 Inf])';

    model = val(1:end-1,:);
    err = val(end,1);
    clear val

    if ~exist('freq_array','var')
        freq_array = logspace(-3,3,30);
    end

    [fwd]=calc_fwd_1d(model(:,1),model(:,2),freq_array,err);
    fwd.rhoerr = (2*fwd.rho.*fwd.Zerr)./abs(fwd.Z);
    fwd.phierr = (180/pi)*(fwd.Zerr)./abs(fwd.Z);

    nskip = 0;
    mode = 0;
    is = 1;

    %These are the variables that matter to MCMC
    MTdata.Z = [real(fwd.Z) imag(fwd.Z)]';
    MTdata.Zerr = [fwd.Zerr fwd.Zerr]';
    MTdata.freq_array = freq_array;

    %These variables are just for housekeeping to remember what parameters were
    %used
    MTdata.name = name; %Filename loaded
    MTdata.synthetic = synthetic;
    MTdata.fwd = fwd; % Apparent resistivity, phase data
    MTdata.err = err; %For synthetic this is the error added, for real data this is the error floor applied
    MTdata.nskip = nskip; %Number of frequencies skipped (for synthetic this is always 0)
    MTdata.mode = mode; % TE = 1, TM = 2, Determinant = 3 (for synthetic this is always 0)
    MTdata.is = is; % Site index (for synthetic data this is always 1)

    %save(OutputFileName,'MTdata')

%% DESIGN LIKELIHOOD AND PRIORS
 % Likelihood---------------------------------
    lognormpdf=@(x,mu,sigma)(-0.5*((x-mu)./sigma).^2 - log(sqrt(2*pi).*sigma));

    P.nL = 4 ; % Number of layers (MCMC is fixed dimension) so you must edit this to be something which is 
                % logical for the data you are using
    P.nparams = 2*P.nL-1; %Number of model parameters

    % Sum the misfit to get the log likelihood probability (note the
    % normalization factor)
    L.logLike=@(m)sum(lognormpdf(MTdata.Z,forwardmodel_melt_z(MTdata.freq_array,m(1:P.nL),m(P.nL+1:end)),MTdata.Zerr));

    L.forward = @(m)forwardmodel_melt_z(MTdata.freq_array,m(1:P.nL),m(P.nL+1:end));

    % Define a function to compute the r.m.s.
    rms_fun = @(dobs,dpred,err)sqrt((1/length(dobs))*nansum(((dobs-dpred)./err).^2));
    L.rms_funm = @(m)rms_fun(MTdata.Z,forwardmodel_melt_z(MTdata.freq_array,m(1:P.nL),m(P.nL+1:end)),MTdata.Zerr);

    % Priors------------------------------------
    % Priors can be made very flexible. In this example, we use a simple upper
    % and lower bound on model parameters as well as a restriction on how thin
    % a layer can be (to avoid completely unresolvable layers).

    %User inputs-------------------------------------------------------
    % Must define upper and lower boundaries on each model parameter
    minrho = 0;
    maxrho = 1;
    mindepth = log10(50);
    maxdepth = log10(50000);

    %The thickness of each layer is restricted to be a function of depth
    % The ratio of each depth must be greater than this
    % value. e.g. if depth_ratio = 1.01, then at 10 km
    % depth, the thinnest allowable layer is 100 m (because
    % 10.1 divided by 10 = 1.01). At 5 km depth, the
    % thinnest allowable layer is 50 m..
    P.depth_ratio = 1.01;

    %----------------------------------------------------------------------



    P.priors = zeros(P.nparams,2);
    P.priors(1:P.nL,:) = repmat([minrho maxrho],P.nL,1);
    P.priors(P.nL+1:end,:) = repmat([mindepth maxdepth],P.nL-1,1);

    % First prior defines bounds on resistivity and depth
    P.logprior{1} = @(m)(all(m>P.priors(:,1))&all(m<P.priors(:,2)));

    % The second prior defines bounds on the thickness as a function of depth
    P.logprior{2} = @(m)(all(diff((m(P.nL+1:end)))>=log10(P.depth_ratio)));

%% MCMC SETTINGS
    S.SaveEvery = 10000000; %Save final data or not (0 = no saving; 1 = save data)
    S.verbose = true; %outputs to command window if verbose

    S.nsteps = 50000; %The length of the Markov Chain for each Walker. 
                     %Common range for good convergence: 100,000 to 1,000,000

    % This parameter is more problem-specific and depends on the number of
    % model paramters you are inverting for. It is recommended to use 2 times
    % the number of model parameters. This determines how many chains (or
    % walkers) you have in your MCMC ensemble.
    S.nWalkers = 20;

    S.stepsize = 2; %Must be >1. Common range: 1.1 to 5. If you are using adaptive MCMC then this is only
                    % the initial value. If you are not using adaptive MCMC,
                    % then this value is fixed.

    S.adaptive = true; %Adaptive MCMC changes step size to optimize at 77% acceptance rate


    %Probably don't need to change anything past here--------------------------

    S.burnin = 0; % Preferably set burnin to 0 and then determine burn-in manually after
    S.thin = 1; % Preferably do not thin chains, but you could maybe do thin = 2 to decrease autocorrelation.


    S.alpha = 1 ; %Set this to add a constant to the acceptance ratio. MCMC operates in ln space so ln(alpha) is what is added.
                  % For this reason, alpha = 1 is a zero. Basically, larger
                  % alpha (e.g. alpha = 100 or alpha = 1000) allows more models
                  % to be accepted with higher r.m.s.. Generally leave alpha = 1

    %Calculated automatically. Do not edit.
    S.totcount = S.nWalkers*S.nsteps; % This is the total number of forward simulations you will need
                                % And can be used to estimate the total time
                                % required.

%% INITIALIZE MODEL
count = 0;
while 1

    %Draw m0 from uniform distribution bounded by the upper and lower
    %bounds of the priors
    m0 = P.priors(:,1) + (P.priors(:,2)-P.priors(:,1)).*rand(P.nparams,S.nWalkers);
    
    % If you are inverting for depth, then you must also sort the depths
    m0(P.nL+1:end,:) = sort(m0(P.nL+1:end,:),1);


    % I also perform a check to ensure that the starting guesses satisfy
    % the priors. (By definition they should satisfy the first prior, but I
    % also check to make sure they satisfy the second prior thickness
    % constraint). I also compute the initial r.m.s. for each walker
    for i = 1:S.nWalkers
        for j = 1:length(P.logprior)
            val{j}(i) = P.logprior{j}(m0(:,i));
        end
        rms0(i) = L.rms_funm(m0(:,i));
    end

    %Check that the priors are satisfied and that the mean intiial r.m.s.
    %for all walkers is less than 50. The r.m.s check is done just to
    %ensure that some models aren't starting from a highly unlikely point.
    %Usually the while loop only needs to run once or twice.
    for j = 1:length(P.logprior)
        prior_satisfied(j) = all(val{j}==1);
    end
    
    if all(prior_satisfied==1) && mean(rms0)<50
        break
    end

    if count > 1000
        %Escape. In this case, you may still be able to run the MCMC sampler, but you
        %might have really high r.m.s. to begin with and might slow down
        %convergence. You may also have a crash if the priors are not
        %satisfied
        disp('Warning: Your priors are not satisifed and/or your initial model is *very* far from the true data, consider adjusting your priors')
        break
    else
        count = count+1;
    end

end

%%
[models, ~,rms] = gwmcmc_jvgr(m0,[P.logprior {L.logLike}],S.totcount ,L.rms_funm,S.alpha,S.adaptive,S.SaveEvery,S.verbose,'burnin',S.burnin,'stepsize',S.stepsize,'thinchain',S.thin);
