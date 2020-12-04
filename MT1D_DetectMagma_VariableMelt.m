%% 1-D SYNTHETIC TESTS FOR DETECTING MAGMA BODIES: VARYING MELT FRACTION
%
% Reference:
%
%   Cordell, D., Hill, G.J., Bachmann, O., Moorkamp, M., Huber, C., in review.
%   Detectability of melt-rich lenses in magmatic reservoirs using
%   magnetotellurics, Geophysical Research Letters
%
% This script replicates the models used by Rasht-Behesht et al. (2020)
% using the MT method. Synthetic data are produced using Wait's recursion
% and then noise are added prior to inversion using Occam algorithm of
% Constable et al. (1987).
%
% This script uses Model Geometry #1 and #2 (thick and thin reservoir,
% respectively) but varies the melt fraction between 0.1 and 0.9
%
% These scripts allow the computation of the data as well as the creation 
% figure 4 in the paper
%
clear all

%Some user defaults
reslims = [1 5000]; 
depthlims = [0 10000];

%Set up the four models
% The depths/thicknesses are taken directly from Rasht-Behesht
input_model{1} = [0 5000 5100 6100 6600];
input_model{2} = [0 5000 5100 5300 5800];

err_to_add = NaN; %NaN triggers frequency-dependent error

%Melt fraction is varied from 0.1 to 0.9
phi = 0.1:0.2:0.9;
% The resistivities are calculated using MAL (m = 1.5, 1000 Ohm m matrix,
% 0.61 Ohm m melt). The melt resistivity is calculated using Guo et al.
% (2016) at 800 C, 100 MPa and 4 wt% H2O.
rho_layer = 1./MAL(phi,1/1000,1/0.61,1.5);

%Loop through the melt levels
for i = 1:length(rho_layer)
     
    filename = ['melt_',sprintf('%02d',i)];
    
    resistivity{1} = [1000 51.8 rho_layer(i) 18.9 1000];
    resistivity{2} = resistivity{1};

    halfspace_value = 1000; %Halfspace in Ohm m
    rms_goal = 1; %Desired rms
   
    disp(['Melt Fraction: ',num2str(phi(i))]);

    tic
    for j = 1:length(input_model) %Loop over the 4 models
        %This function computes the forward data and performs the OCCAM
        %inversion
        %
        % This function outputs a "mod" structure which contains all
        % the information about the inversion
        mod(j,1) = OCCAM1D_simple(input_model{j},resistivity{j},err_to_add,halfspace_value,rms_goal);

    end

    toc
    %Automatically save the results as a matfile
    save(filename);

end

%% PLOT RESULTS
clear all
load('melt_01'); 
eidx = 1;
for E = 1:length(phi)
    
    load(['melt_',sprintf('%02d',E)]); 
    reslims = [0.1 5000]; 
    for   M = 1:2 %Loop through the 2 models 
      
        [~,id2] = min(abs(2500-mod(M,1).depth)); %Find the index at 2 km
        [~,id10] = min(abs(10000-mod(M,1).depth)); %Index at 10 km
        [~,id5] = min(abs(5000-mod(M,1).depth)); %Index at 5 km

        %Get all the 100 models into a single matrix (nlayers x nruns)
        allmod = log10(mod(M,1).model);
        [allminmod,~] = nanmin(allmod(id2:id10,:)); %Find the minimum bulk resistivity between 2 and 10 km depth
    
        %Calculate the mean and std of all 100 runs
        minmod = 10.^nanmean(allminmod);
        stdmod = 10.^nanstd(allminmod);
      

        %Calculate the estimated melt fraction for all 100 runs
        %To solve for phi using MAL, it must be done numerically. Here
        %the function uses Newton-Raphson
        melt_frac = MAL_solve_phi(10.^allminmod,1000,0.61,1.5);
     
        fullmeltfrac(M,E) = melt_frac;
    end
end

layer = [1000 200];
for M = 1:2
    subplot(1,2,M)
    plot(phi,fullmeltfrac(M,:),'.k','MarkerSize',8); hold on
    xlabel('True Melt Fraction');
    ylabel('Estimated Melt Fraction');
    title(['Model Geometry #',num2str(M),': ',num2str(layer(M)),' m thick'])
    plot([0 1],[0 1],'--r'); axis equal
    axis([0 1 0 1])
    grid on
end



