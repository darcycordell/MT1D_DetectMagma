function [models,logP,rms,reject]=gwmcmc(minit,logPfuns,mccount,stats_fun,alpha,adaptive,SaveEvery,verbose,varargin)
%% Cascaded affine invariant ensemble MCMC sampler. "The MCMC hammer"
%
% GWMCMC is an implementation of the Goodman and Weare 2010 Affine
% invariant ensemble Markov Chain Monte Carlo (MCMC) sampler. MCMC sampling
% enables bayesian inference. The problem with many traditional MCMC samplers
% is that they can have slow convergence for badly scaled problems, and that
% it is difficult to optimize the random walk for high-dimensional problems.
% This is where the GW-algorithm really excels as it is affine invariant. It
% can achieve much better convergence on badly scaled problems. It is much
% simpler to get to work straight out of the box, and for that reason it
% truly deserves to be called the MCMC hammer.
%
% (This code uses a cascaded variant of the Goodman and Weare algorithm).
%
% USAGE:
%  [models,logP]=gwmcmc(minit,logPfuns,mccount, Parameter,Value,Parameter,Value);
%
% INPUTS:
%     minit: an MxW matrix of initial values for each of the walkers in the
%            ensemble. (M:number of model params. W: number of walkers). W
%            should be atleast 2xM. (see e.g. mvnrnd).
%  logPfuns: a cell of function handles returning the log probality of a
%            proposed set of model parameters. Typically this cell will
%            contain two function handles: one to the logprior and another
%            to the loglikelihood. E.g. {@(m)logprior(m) @(m)loglike(m)}
%   mccount: What is the desired total number of monte carlo proposals.
%            This is the total number, -NOT the number per chain.
%
% Named Parameter-Value pairs:
%   'StepSize': unit-less stepsize (default=2.5).
%   'ThinChain': Thin all the chains by only storing every N'th step (default=10)
%   'ProgressBar': Show a text progress bar (default=true)
%   'Parallel': Run in ensemble of walkers in parallel. (default=false)
%   'BurnIn': fraction of the chain that should be removed. (default=0)
%
% OUTPUTS:
%    models: A MxWxT matrix with the thinned markov chains (with T samples
%            per walker). T=~mccount/p.ThinChain/W.
%      logP: A PxWxT matrix of log probabilities for each model in the
%            models. here P is the number of functions in logPfuns.
%
% Note on cascaded evaluation of log probabilities:
% The logPfuns-argument can be specifed as a cell-array to allow a cascaded
% evaluation of the probabilities. The computationally cheapest function should be
% placed first in the cell (this will typically the prior). This allows the
% routine to avoid calculating the likelihood, if the proposed model can be
% rejected based on the prior alone.
% logPfuns={logprior loglike} is faster but equivalent to
% logPfuns={@(m)logprior(m)+loglike(m)}
%
% TIP: if you aim to analyze the entire set of ensemble members as a single
% sample from the distribution then you may collapse output models-matrix
% thus: models=models(:,:); This will reshape the MxWxT matrix into a
% Mx(W*T)-matrix while preserving the order.
%
%
% EXAMPLE: Here we sample a multivariate normal distribution.
%
% %define problem:
% mu = [5;-3;6];
% C = [.5 -.4 0;-.4 .5 0; 0 0 1];
% iC=pinv(C);
% logPfuns={@(m)-0.5*sum((m-mu)'*iC*(m-mu))}
%
% %make a set of starting points for the entire ensemble of walkers
% minit=randn(length(mu),length(mu)*2);
%
% %Apply the MCMC hammer
% [models,logP]=gwmcmc(minit,logPfuns,100000);
% models(:,:,1:floor(size(models,3)*.2))=[]; %remove 20% as burn-in
% models=models(:,:)'; %reshape matrix to collapse the ensemble member dimension
% scatter(models(:,1),models(:,2))
% prctile(models,[5 50 95])
%
%
% References:
% Goodman & Weare (2010), Ensemble Samplers With Affine Invariance, Comm. App. Math. Comp. Sci., Vol. 5, No. 1, 65–80
% Foreman-Mackey, Hogg, Lang, Goodman (2013), emcee: The MCMC Hammer, arXiv:1202.3665
%
% WebPage: https://github.com/grinsted/gwmcmc
%
% -Aslak Grinsted 2015


    
persistent isoctave;  
if isempty(isoctave)
	isoctave = (exist ('OCTAVE_VERSION', 'builtin') > 0);
end

if nargin<7
    error('GWMCMC:toofewinputs','GWMCMC requires atleast 5 inputs.')
end
M=size(minit,1);
if size(minit,2)==1
    minit=bsxfun(@plus,minit,randn(M,M*5));
end

if ~exist('verbose','var')
    verbose = true;
end


p=inputParser;
if isoctave
    p=p.addParamValue('StepSize',2,@isnumeric); %addParamValue is chosen for compatibility with octave. Still Untested.
    p=p.addParamValue('ThinChain',10,@isnumeric);
    p=p.addParamValue('ProgressBar',true,@islogical);
    p=p.addParamValue('Parallel',false,@islogical);
    p=p.addParamValue('BurnIn',0,@(x)(x>=0)&&(x<1));
    p=p.parse(varargin{:});
else
    p.addParameter('StepSize',2,@isnumeric); %addParamValue is chose for compatibility with octave. Still Untested.
    p.addParameter('ThinChain',10,@isnumeric);
    p.addParameter('ProgressBar',true,@islogical);
    p.addParameter('Parallel',false,@islogical);
    p.addParameter('BurnIn',0,@(x)(x>=0)&&(x<1));
    p.parse(varargin{:});
end
p=p.Results;


Nwalkers=size(minit,2);

if size(minit,1)*2>size(minit,2)
    warning('GWMCMC:minitdimensions','Check minit dimensions.\nIt is recommended that there be atleast twice as many walkers in the ensemble as there are model dimension.')
end

if p.ProgressBar
    progress=@textprogress;
else
    progress=@noaction;
end


Nkeep=ceil(mccount/p.ThinChain/Nwalkers); %number of samples drawn from each walker
mccount=(Nkeep-1)*p.ThinChain+1;

models=nan(M,Nwalkers,Nkeep); %pre-allocate output matrix
rms = nan(Nwalkers,Nkeep);

models(:,:,1)=minit;


if ~iscell(logPfuns)
    logPfuns={logPfuns};
end

NPfun=numel(logPfuns);

%calculate logP state initial pos of walkers
logP=nan(NPfun,Nwalkers,Nkeep);
for wix=1:Nwalkers
    for fix=1:NPfun
        v=logPfuns{fix}(minit(:,wix));
        if islogical(v) %reformulate function so that false=-inf for logical constraints.
            
            %Original Download from Github had the following line. It seems
            %like a bug to me because it overwrites logPfuns{fix} as a new
            %function and you lose the information related to the original logPfuns{fix}
            %v=-1/v;logPfuns{fix}=@(m)-1/logPfuns{fix}(m); %experimental implementation of experimental feature
            
            %Fixed to have logPfuns{fix} held as a function within a
            %function
            funholder{fix} = logPfuns{fix};
            %v=-1/v;
            logPfuns{fix}=@(m)log(double(funholder{fix}(m)));
            v = logPfuns{fix}(minit(:,wix));
        end
        logP(fix,wix,1)=v;
        
    end
end

if ~all(all(isfinite(logP(:,:,1))))
    error('Starting points for all walkers must have finite logP')
end

%%
reject=zeros(Nwalkers,1);
%reject_count = zeros(Nwalkers,1);
%bigStep = zeros(1,Nwalkers);
reject_first = zeros(size(reject));
reject_second = zeros(NPfun,Nwalkers);

% min_step = p.StepSize;
% max_step = 5;
% max_rms = 5;
% p.StepSize = max_step*ones(1,Nwalkers);
% rms(:,1) = 10^6*ones(Nwalkers,1);
p.StepSize = p.StepSize*ones(1,Nwalkers);
curm=models(:,:,1);
curlogP=logP(:,:,1);
if verbose
    progress(0,0,0,0,0,0);
end
rejectrate = ones(1,Nwalkers);
timestep = 0; timeavg = 0;
totcount=Nwalkers; %ind = [1:Nwalkers];
for row=1:Nkeep
    time= tic;
    for jj=1:p.ThinChain
  
        %generate proposals for all walkers
        %(done outside walker loop, in order to be compatible with parfor - some penalty for memory):
        %-Note it appears to give a slight performance boost for non-parallel.
        %rix=mod((1:Nwalkers)+floor(rand*(Nwalkers-1)),Nwalkers)+1; %pick a random partner
        
        %bigStep(reject_count>10) = 10;%p.StepSize;
        
        rix=randpermfull(Nwalkers);
        zz=(((p.StepSize - 1).*rand(1,Nwalkers) + 1).^2./(p.StepSize));
        %zz=(((bigStep+p.StepSize - 1).*rand(1,Nwalkers) + 1).^2./(bigStep+p.StepSize));
        %zz = (((p.StepSize^2 - 1)*rand(1,Nwalkers) + 1)/p.StepSize).^2; %variation
        proposedm=curm(:,rix) - bsxfun(@times,(curm(:,rix)-curm),zz);
        logrand=log(rand(NPfun+1,Nwalkers)); %moved outside because rand is slow inside parfor

        for wix=1:Nwalkers
            acceptfullstep=true;
            proposedlogP=nan(NPfun,1);
            if logrand(1,wix)<(numel(proposedm(:,wix))-1)*log(zz(wix))
                for fix=1:NPfun
                    proposedlogP(fix)=logPfuns{fix}(proposedm(:,wix));
                    if logrand(fix+1,wix)>log(alpha)+(proposedlogP(fix)-curlogP(fix,wix)) || ~isreal(proposedlogP(fix)) || isnan(proposedlogP(fix))
                        acceptfullstep=false;
                        reject_second(fix,wix) = reject_second(fix,wix)+1;
                        break
                    end
                end
            else
                reject_first(wix) = reject_first(wix)+1;
                acceptfullstep=false;
            end
            if acceptfullstep
                curm(:,wix)=proposedm(:,wix); curlogP(:,wix)=proposedlogP;
                %reject_count(wix) = 0;
                %bigStep(wix) = 0;
            else
                %reject_count(wix) = reject_count(wix)+1;
                reject(wix)=reject(wix)+1;
            end
            rms(wix,row) = stats_fun(curm(:,wix));
            
            rejectrate(wix) = reject(wix)./(totcount/Nwalkers);
            
            if adaptive
                if rejectrate(wix)>0.77
                    p.StepSize(wix) = p.StepSize(wix)*0.9;
                     if p.StepSize(wix)<=1
                        p.StepSize(wix) = 1.01;
                    end
                else
                    p.StepSize(wix) = p.StepSize(wix)*1.1;
                end

                if rand>0.9
                    max_step = 5; min_step = 1.01;
                    max_rms = 4;
                    k = log(max_step/min_step)/(max_rms-1);
                    A = min_step/exp(k);

                    p.StepSize(wix) = min(A*exp(k*rms(wix,row)),max_step);
                end
            end
%                 if row>1 && acceptfullstep
%                     if abs(rms(wix,row)-rms(wix,row-1))<0.001
%                         s = 1;
%                     end
%                 end
        end
        

        
        
        %length(find(bigStep~=0))
%         k = log(max_step/min_step)/(max_rms-1);
%         A = min_step/exp(k);
%         
%         p.StepSize = A*exp(k*rms(:,row)');
%         
%         p.StepSize(rms(:,row)>max_rms)=max_step;
        
        %p.StepSize'
        
        totcount=totcount+Nwalkers;

        if verbose
            progress((row-1+jj/p.ThinChain)/Nkeep,curm,sum(reject)/totcount,rms(:,row),p.StepSize,timeavg)
        end
    end
    models(:,:,row)=curm;
    logP(:,:,row)=curlogP;
    %ind = find(rms(:,row)==min(rms(:,row)));

    %progress bar
    timestep = toc(time)+timestep;

    timeavg = timestep/totcount*Nwalkers;
    
    
    if mod(row,SaveEvery)==0
        
        %save('MCMC_Output.mat','models','rms')
        
    end


end
if verbose
    progress(1,0,0,0);
end
if p.BurnIn>0
    crop=ceil(Nkeep*p.BurnIn);
    models(:,:,1:crop)=[]; %TODO: never allocate space for them ?
    logP(:,:,1:crop)=[];
end

save('MCMC_Output.mat','models','rms')

pause(0.1)

end


% TODO: make standard diagnostics to give warnings...
% TODO: make some diagnostic plots if nargout==0;



function textprogress(pct,curm,rejectpct,rms,reject,timeavg)
persistent lastNchar lasttime starttime
if isempty(lastNchar)||pct==0
    lasttime=cputime-10;starttime=cputime;lastNchar=0;
    pct=1e-16;
end
if pct==1
    fprintf('%s',repmat(char(8),1,lastNchar));lastNchar=0;
    return
end
if (cputime-lasttime>0.1)

    ETA=datestr((cputime-starttime)*(1-pct)/(pct*60*60*24),13);
    progressmsg=[183-uint8((1:40)<=(pct*40)).*(183-'*') ''];
    %progressmsg=['-'-uint8((1:40)<=(pct*40)).*('-'-'•') ''];
    %progressmsg=[uint8((1:40)<=(pct*40)).*'#' ''];
    %curmtxt=sprintf('% 9.3g\n',[curm(1:min(end,20),1); rms(1)]);
    pval(1:2:2*length(rms)) = rms;
    pval(2:2:2*length(rms)) = reject;
    curmtxt=sprintf('% 9.3g %9.3g\n',pval(1:length(rms)),pval(length(rms)+1:end));
    
    
    %curmtxt=mat2str(curm);
    progressmsg=sprintf('\nGWMCMC %5.1f%% [%s] %s\n%3.0f%% rejected. Average time per iteration = %5.6f \n%s\n',pct*100,progressmsg,ETA,rejectpct*100,timeavg,curmtxt);

    fprintf('%s%s',repmat(char(8),1,lastNchar),progressmsg);
    drawnow;lasttime=cputime;
    lastNchar=length(progressmsg);
end
end

function noaction(varargin)

% Acknowledgements: I became aware of the GW algorithm via a student report
% which was using emcee for python. Great stuff.
end

