clearvars
geometrynum = 4;
% Load appropriate model
fid = fopen(['geometry_',num2str(geometrynum),'_err_005.txt']);
load(['geometry_',num2str(geometrynum),'_err_005_FRun02.mat']);
load(['geometry_',num2str(geometrynum),'_FRun02_MCMC_Output.mat']);

if geometrynum ~=4
    load(['geometry_',num2str(geometrynum),'_FRun02_Occam.mat']);
end


val = fscanf(fid,'%f %f',[2 Inf])';
orig_model = val(1:end-1,:);
fclose(fid);

[occ.phi,occ.numiter] = solve_phi_MAL(min(occ.model(occ.depth>2000)),1000,0.61,1.5);
doc2000 = nearestpoint(2000,occ.depth);
[~, doccid] = min(occ.model(doc2000:end));%Minimum res
doccid = doc2000+doccid-1;

burnin = 200000;

nanind = all(isnan(rms),1);
rms(:,nanind) =[];
%Find "good" walkers which converged to a well-mixed low rms
target_rms = 1.5; %Play with this number
                  % Start large (e.g. 100) to see all walkers.
ind = find(rms(:,200000)<=target_rms);

rhob = 1./MAL(1/1000,1/0.61,1.5,models(1:P.nL,:,:));
tic
rhoa = zeros(length(freq_array),length(1:2:length(ind)),length(burnin:1000:S.nsteps));
phia = zeros(size(rhoa));
for i = 1:2:length(ind)
    for j = 1:1000:S.nsteps-burnin+1
        [fwd]=calc_fwd_1d([0; 10.^models(P.nL+1:end,ind(i),j+burnin-1)],rhob(:,ind(i),j+burnin-1),freq_array,0);
        rhoa(:,i,j) = fwd.rho;
        phia(:,i,j) = fwd.phi;
    end
    disp(['Walker #',num2str(i)])
end
toc

%% SUPPLEMENTARY FIGURE: RMS CONVERGENCE PLOT
%Usually the first iterations have very high rms so I plot two plots: the
%first from 1:iterplot and the second from iterplot:end. The burnin period
%(i.e. the iterations it takes to converge to low rms)
%can be anywhere from 10% - 40% of the full run.
iterplot = 200000;
for i = 1:length(ind)
    subplot(3,1,1)
    plot(1:iterplot,squeeze(rms(ind(i),1:iterplot))','-'); hold on
    title(['r.m.s. For First ',num2str(iterplot),' Iterations'])
    ylabel('r.m.s. misfit')
    
    subplot(3,1,2)
    plot(iterplot:length(rms),squeeze(rms(ind(i),iterplot:end))','-'); hold on
    %axis([iterplot length(rms) 0.8 1.1])
    xlabel('Iteration #')
    ylabel('r.m.s. misfit')
    title(['r.m.s. From ',num2str(iterplot),' Iterations Onward'])
end

% RMS HISTOGRAM PLOT
subplot(3,1,3)
histogram(rms(ind,burnin:end))
xlabel('r.m.s. misfit');
ylabel('Count')

%% SUPPLEMENTARY FIGURE: PLOT RAW HISTOGRAMS OF MODEL PARAMETERS
figure(1)
for i = 1:P.nL
    subplot(2,ceil(P.nparams/2),i)
    h = histogram(models(i,ind,burnin:end),'FaceColor','k','FaceAlpha',1); hold on
    %plot([log10(model_res(i)) log10(model_res(i))],[0 max(h.Values)],'-r','LineWidth',2)
    title(['Layer #',num2str(i)])
    xlabel(['Melt Fraction']); ylabel('Count')
    %axis([-2 5 0 10000])
end

for i = P.nL+1:2*P.nL-1
    subplot(2,ceil(P.nparams/2),i+1)
    h = histogram(10.^models(i,ind,burnin:end),'FaceColor','k','FaceAlpha',1); hold on
    %plot([log10(model_res(i)) log10(model_res(i))],[0
    %max(h.Values)],'-r','LineWidth',2)
    xlabel('Depth'); ylabel('Count')
    %axis([-2 5 0 10000])
end

%% TABLE 1: SORT OUT STATISTICS FOR MODEL PARAMETERS

[mler,mle_plusr,mle_minusr,~,~] = mle_calc(log10(rhob(2,ind,burnin:end)),0.1,1);

[mled,mle_plusd,mle_minusd,~,~] = mle_calc(10.^(models(2*P.nL-2,ind,burnin:end)),0.1,1);
[mlet,mle_plust,mle_minust,~,~] = mle_calc(10.^(models(2*P.nL-1,ind,burnin:end))-10.^(models(2*P.nL-2,ind,burnin:end)),0.1,25);

[mlem,mle_plusm,mle_minusm,~,~] = mle_calc((models(P.nL-1,ind,burnin:end)),0.1,1);

excel_table = [10^mler, 10^mle_minusr, 10^mle_plusr; mled, mle_minusd, mle_plusd; mlet, mle_minust, mle_plust; mlem, mle_minusm, mle_plusm];

excel_table(1,4) = occ.model(doccid);
excel_table(2,4) = occ.depth(doccid);
excel_table(3,4) = NaN;
excel_table(4,4) = occ.phi;


%% FIGURE 2: GRID RESULTS AND PLOT MODEL AND RESPONSES
mp = log10(rhob(:,ind,burnin:end));
mp=mp(:,:)'; %flatten the chain


md = 10.^(models(P.nL+1:end,ind,burnin:end));
md = md(:,:)';
md = [zeros(length(md),1) md];

md = md(:);
mp = mp(:);

minrho = P.priors(1,1); maxrho = P.priors(1,2);
maxdepth = P.priors(P.nL+1,2);

minrho = -1; maxrho = 4;
     
rhogrid = linspace(minrho,maxrho,200); %Gridded log resistivity
dgrid = linspace(-100,10^maxdepth,200); %Gridded depths (meters)

[count_interfaces,Xedges,Yedges] = histcounts2(mp,md,rhogrid,dgrid);


% PLOT COLORMAP OF LOG PROBABILITY
axlim = [-1 4 -100 10000];

uc = ceil(max(max(log10(count_interfaces./(length(mp)))))); %Upper colorbar limit (log(p))
lc = floor(log10(1./length(mp))); %Lower colorbar limit (log(p))

set_figure_size(2);
subplot(2,2,[2 4])
pcolor(rhogrid(1:end-1),dgrid(1:end-1),log10(count_interfaces'./(length(mp)))); hold on
axis(axlim); shading flat;
axis ij

stairs(log10(occ.model),occ.depth,'-m','LineWidth',1);
stairs(log10([orig_model(1,2); orig_model(:,2); orig_model(end,2)]),[orig_model(:,1); orig_model(end,1); 10^6],'-k','LineWidth',1); hold on
manual_legend('Original Model','-k','Occam Model','-m');
     
title(['Geometry #',num2str(geometrynum)]);

hcb = colorbar;
hcb.Label.String = 'Log10(Probability Density)';
colormap(flipud(pink(uc-lc))); caxis([lc uc]);
grid on
ylabel('Depth (m)'); 
xlabel('Log10 Resistivity')
set(gca,'layer','top')
%set(gca,'XScale','log')


subplot(2,2,1); 
for i = 1:size(rhoa,2)
    for j = 1:10000:size(rhoa,3)
        loglog(1./freq_array,rhoa(:,i,j),'-','Color',[0.8 0.8 0.8]);hold on
    end
end
loglog(1./freq_array,MTdata.fwd.rho,'.k','MarkerSize',12);
loglog(1./freq_array,occ.ra_mod,'-b','LineWidth',1)
xlabel('Period (s)')
ylabel('Resistivity (\Omega m)')
manual_legend('Synthetic Data 5% Noise','.k','Occam Respose','-b','Bayesian Responses','-k');
grid on

subplot(2,2,3)
for i = 1:size(rhoa,2)
    for j = 1:10000:size(rhoa,3)
        loglog(1./freq_array,phia(:,i,j),'-','Color',[0.8 0.8 0.8]);hold on
    end
end
semilogx(1./freq_array,MTdata.fwd.phi,'.k','MarkerSize',12); hold on
semilogx(1./freq_array,occ.ph_mod,'-b','LineWidth',1)
xlabel('Period (s)');
ylabel('Phase (deg)');
grid on



%% FIGURE 3: PROBABILITY MAP FOR MELT FRACTION AND DEPTH
% FIGURE 4: PROBABLITY MAP OF MELT FRACTION AND THICKNESS

true_depth = 5100;
thick_flag = 1;

if geometrynum == 1 || geometrynum == 3
    true_thick = 1000;
else
    true_thick = 200;
end


tic


nL = P.nL;

if nL == 3
    mod_ind = 2;
elseif nL == 4
    mod_ind = 3;
end

mp = models(mod_ind,ind,burnin:end);
mp=mp(:,:)'; %flatten the chain




if thick_flag
    md = 10.^(models(mod_ind+nL,ind,burnin:end))-10.^(models(mod_ind+nL-1,ind,burnin:end));
else
    md = 10.^(models(mod_ind+nL-1,ind,burnin:end));
end
md = md(:,:)';
%md = [zeros(length(md),1) md];

md = md(:);
mp = mp(:);

minrho = P.priors(1,1); maxrho = P.priors(1,2);
maxdepth = P.priors(nL+1,2);

if thick_flag
    maxdepth = log10(4000);
end
     
rhogrid = linspace(minrho,maxrho,200); %Gridded log resistivity
dgrid = linspace(0,10^maxdepth,400); %Gridded depths (meters)

md(md>10.^maxdepth) = 10.^maxdepth;

[count_interfaces,Xedges,Yedges] = histcounts2(mp,md,rhogrid,dgrid);

 toc
% PLOT COLORMAP OF LOG PROBABILITY
synthetic = 0;
if thick_flag
    axlim = [0 1 0 4000];
else
    axlim = [0 1 3000 7000];
end

uc = ceil(max(max(log10(count_interfaces./(length(mp)))))); %Upper colorbar limit (log(p))
lc = floor(log10(1./length(mp))); %Lower colorbar limit (log(p))

figure(10);
pcolor(rhogrid(1:end-1),dgrid(1:end-1),log10(count_interfaces'./(length(mp)))); hold on
title(['Geometry #',num2str(geometrynum)]);
if synthetic
    stairs(log10([model_res(1); model_res'; model_res(end)]),[10^-6 depth(2:end) depth(end) 10^6],'-k','LineWidth',2); hold on
end
axis(axlim); shading flat; 

if thick_flag
    plot(0.3,true_thick,'oy','MarkerSize',10,'MarkerFaceColor','y','MarkerEdgeColor','k')
    plot([occ.phi occ.phi],[axlim(3) axlim(4)],'--k')
else
    plot(0.3,true_depth,'oy','MarkerSize',10,'MarkerFaceColor','y','MarkerEdgeColor','k')
    plot(occ.phi,(occ.depth(doccid)+occ.depth(doccid+1))/2,'om','MarkerSize',10,'MarkerFaceColor','m','MarkerEdgeColor','k')
end

if ~thick_flag
    axis ij
end
hcb = colorbar;
hcb.Label.String = 'Log10(Probability Density)';
colormap(flipud(pink(uc-lc))); caxis([lc uc]);
grid on
if thick_flag
    ylabel('Thickness (m)'); 
else
    ylabel('Depth (m)'); 
end
xlabel('Melt Fraction')
set(gca,'layer','top')

%s071c = load('res_curve_LDMs071c.mat');

%plot(log10(s071c.res(:,4)),s071c.m.cz+2160,'-b','LineWidth',2);

%% FIGURE 5: PLOT VOLUMES
clearvars
for imodel = 1:4
    geometrynum = imodel;
    % Load appropriate model
    fid = fopen(['geometry_',num2str(geometrynum),'_err_005.txt']);
    load(['geometry_',num2str(geometrynum),'_err_005_FRun02.mat']);
    load(['geometry_',num2str(geometrynum),'_FRun02_MCMC_Output.mat']);

    if geometrynum ~=4
        load(['geometry_',num2str(geometrynum),'_FRun02_Occam.mat']);
    end

    val = fscanf(fid,'%f %f',[2 Inf])';
    orig_model = val(1:end-1,:);
    fclose(fid);

    [occ.phi,occ.numiter] = solve_phi_MAL(min(occ.model(occ.depth>2000)),1000,0.61,1.5);
    doc2000 = nearestpoint(2000,occ.depth);
    [~, doccid] = min(occ.model(doc2000:end));%Minimum res
    doccid = doc2000+doccid-1;

    burnin = 200000;

    nanind = all(isnan(rms),1);
    rms(:,nanind) =[];
    %Find "good" walkers which converged to a well-mixed low rms
    target_rms = 1.5; %Play with this number
                      % Start large (e.g. 100) to see all walkers.
    ind = find(rms(:,200000)<=target_rms);

    true_depth = 5100;
    thick_flag = 1;

    if geometrynum == 1 || geometrynum == 3
        true_thick = 1000;
    else
        true_thick = 200;
    end

    nL = P.nL;

    if nL == 3
        mod_ind = 2;
    elseif nL == 4
        mod_ind = 3;
    end

    mp = models(mod_ind,ind,burnin:end);
    mp=mp(:,:)'; %flatten the chain

    if thick_flag
        md = 10.^(models(mod_ind+nL,ind,burnin:end))-10.^(models(mod_ind+nL-1,ind,burnin:end));
    else
        md = 10.^(models(mod_ind+nL-1,ind,burnin:end));
    end
    md = md(:,:)';
    %md = [zeros(length(md),1) md];

    md = md(:);
    mp = mp(:);

    %plot(log10(md),log10(mp),'.k')

    v = md.*mp;

    %v = md.*mp;
    if imodel == 4
        v(v>316)=NaN;
    end

    subplot(2,2,imodel)
    histogram(v);
    xlabel('Volume of Melt (m^3)')
    ylabel('Count')
    title(['Geometry #',num2str(imodel)])
    grid on
end



