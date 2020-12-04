%% 1-D SYNTHETIC TESTS FOR DETECTING MAGMA BODIES
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
% There are 4 different noise levels used (2.5%, 5%, 10%, and
% frequency-dependent noise) as well as 4 different models (thin and thick
% melt-rich layer; with and without overlying sediment fill). Thus there
% are a total of 16 different (model,noise) combinations. 100 inversions
% are performed for each combination with different sampling of noise.
%
% These scripts allow the computation of the data as well as the creation 
% figures 2 and 3 in the paper as well as all the 15 supplementary figures
%
clear all

%Some user defaults
reslims = [1 5000]; 
depthlims = [0 10000];
numruns = 5; %Each error level is run 100 times (1600 total inversions calculations)


%Set up the four models
% The depths/thicknesses are taken directly from Rasht-Behesht
% The resistivities are calculated using MAL (m = 1.5, 1000 Ohm m matrix,
% 0.61 Ohm m melt). The melt resistivity is calculated using Guo et al.
% (2016) at 800 C, 100 MPa and 4 wt% H2O.
input_model{1} = [0 5000 5100 6100 6600]; %1000 m thick, no cover
input_model{2} = [0 5000 5100 5300 5800]; %200 m thick, no cover
input_model{3} = [0 500 5000 5100 6100 6600]; %1000 m thick, with cover
input_model{4} = [0 500 5000 5100 5300 5800]; %200 m thick with cover

resistivity{1} = [1000 51.8 3.7 18.9 1000];
resistivity{2} = resistivity{1};
resistivity{3} = [10 1000 51.8 3.7 18.9 1000];
resistivity{4} = resistivity{3};

%There are four error levels to be added to the data.
filename = {'01_low','02_medium','03_high','04_variable'};
err = [0.025, 0.05, 0.1, NaN]; %If err = NaN then this triggers the frequency-dependent noise.

%Loop through error levels
for i = 1:length(err)

    err_to_add = err(i); %Gaussian error to add
    halfspace_value = 1000; %Halfspace in Ohm m
    rms_goal = 1; %Desired rms
   
    disp(['Error level: ',num2str(err(i))]);
    for j = 1:numruns %Create 100 different random data sets with the prescribed noise level
        tic
        for k = 1:length(input_model) %Loop over the 4 models
            %This function computes the forward data and performs the OCCAM
            %inversion. Outputs a "mod" data structure which contains all
            % the information about the inversion
            mod(k,j) = OCCAM1D_simple(input_model{k},resistivity{k},err_to_add,halfspace_value,rms_goal);

        end
        disp(['RUN #',num2str(j)])
        toc
    end
    
    %Automatically save the results as a matfile
    save(filename{i});

end

%% PLOT RESULTS FOR FIGURE 2 AND 3
clear all
plotfig =true; %Set this as true to recreate Figure 2 and all the supplementary figures (S1 to S15)
filename = {'01_low','02_medium','03_high','04_variable'}; %The matfile names which must be in the current directory
colorchart = [0 0 0; 1 0 0; 0 1 0; 0 0 1]; %Colors for each error level
LT = {'o','s','*','v'}; %Line type for each error level
MS = [5,5,5,5]; %Marker size for each error level
count = 0;
for i = [3 2 4 1] %Loop through the 4 error cases from highest noise (3) to lowest noise (1). The variable error (4) lies between 1 and 2

    load(filename{i}); %File must be in current directory

    elevel = err(i);
    if isnan(elevel)
        elevel = 'Variable';
    end
    
    for   k = 1:length(input_model) %Loop through the 4 models
        
        if plotfig %If plotting
            %Set figure size:
            screensize=get(groot,'Screensize');
            fig=figure(1); clf
            set(fig,'Position',[0.1*screensize(3) 0.1*screensize(4) 0.8*screensize(3) 0.8*screensize(4)])
            
            gray = [150 150 150]/255;

            for j = 1:numruns %Loop through the 100 runs and plot each inversion model

                m = mod(k,j);

                %PLOT RESULTS--------------------------------------------------------------
                subplot(2,2,1)
                loglog(m.T,m.ra_mod,'Color',gray,'LineWidth',1); hold on

                subplot(2,2,3)
                semilogx(m.T,m.ph_mod,'Color',gray,'LineWidth',1); hold on

                subplot(2,2,[2 4])
                stairs(m.model,m.depth/1000,'Color',gray,'LineWidth',1); hold on

            end

            %Plot apparent resistivity and phase data for one of the
            %datasets
            subplot(2,2,1)
            [h, g] = logerrorbar(m.T,m.rhoa',m.rhoerr,'.r','-r'); hold on
            set(h,'MarkerSize',12)
            set(g,'LineWidth',2)
            loglog(m.T,halfspace_value*ones(length(m.T),1),'--r','LineWidth',2);
            xlabel('Period (s)')
            ylabel('App Res (\Omega m)')
            axis([min(m.T) max(m.T) reslims])
            grid on
            title([num2str(length(mod(1,:))),' Runs with ',num2str(err_to_add),' Noise Added'])

            subplot(2,2,3)
            errorbar(m.T,m.pha',m.phaerr,'.r','MarkerSize',12,'LineWidth',1); hold on
            semilogx(m.T,45*ones(length(m.T),1),'--r','LineWidth',2)
            set(gca,'XScale','log')
            xlabel('Period (s)')
            ylabel('Phase (deg)')
            axis([min(m.T) max(m.T) 0 90])
            grid on

            %Plot the true resistivity model
            subplot(2,2,[2 4])
            stairs([resistivity{k}(1) resistivity{k} resistivity{k}(end)],[input_model{k} input_model{k}(end) 10^6]/1000,'-r','LineWidth',3); hold on
            axis([reslims depthlims/1000])
                        set(gca,'XScale','log')
            xlabel('Resistivity (\Omega m)')
            ylabel('Depth (km)')
            axis ij
            grid on

        end

       
        %Calculate the mean and stddev predicted bulk resistivity for all
        %the inverse models
        meanmod = nanmean(log10(reshape([mod(k,:).model],length(mod(k,1).depth),length(mod(k,:)))),2);
        stdmod = nanstd(log10(reshape([mod(k,:).model],length(mod(k,1).depth),length(mod(k,:)))),0,2);

        if plotfig
            subplot(2,2,[2 4])%Plot the mean model from all 100 runs
            stairs(10.^meanmod,m.depth/1000,'-k','LineWidth',3); hold on
            manual_legend('True Model','-r','Average Predicted Model','-k');
        end
      
        [~,id2] = min(abs(2500-mod(k,1).depth)); %Find the index at 2 km
        [~,id10] = min(abs(10000-mod(k,1).depth)); %Index at 10 km
        [~,id5] = min(abs(5000-mod(k,1).depth)); %Index at 5 km

        %Get all the 100 models into a single matrix (nlayers x nruns)
        allmod = log10(reshape([mod(k,:).model],length(mod(k,1).depth),length(mod(k,:))));
        [allminmod,allid] = nanmin(allmod(id2:id10,:)); %Find the minimum bulk resistivity between 2 and 10 km depth
        allid = allid+id2-1; %The index of the minimum bulk resistivity
        for L = 1:length(allid) %Loop through all the model runs to find the depth to the melt-rich layer
            %Note that the depth to the layer is necessarily defined at a
            %model edge (e.g. 5 km). So the model cell above and below the
            %model edge both have equal validity in determining the
            %"depth". Really, any depth between 4.5 km and 5.5 km is
            %"detecting" the depth of the layer within the mesh dimensions.
            alldepth(L) = mod(k,1).depth(allid(L)-1);
            if alldepth(L)<2000 %In some (high-noise) cases, the layer is not detected at all
                %In this case, set the depth to NaN
                alldepth(L) = NaN;
            end
        end

        for L = 1:length(mod(k,1).model) %Loop through model layers and find how many model runs found the minimum
            %at a given depth.
            nmin_layer(L) = sum((alldepth == mod(k,1).depth(L)));
        end

        %Fill a 16 x 16 matrix with the number of models which found the
        %correct depth within 4.5 km and 5.5 km depth.
        layers(k,i) = nmin_layer(id5) + nmin_layer(id5-1);

        %Calculate the mean and std of all 100 runs
        minmod = 10.^nanmean(allminmod);
        stdmod = 10.^nanstd(allminmod);
        if plotfig
            title(['Estimated Depth: ',num2str(nanmean(alldepth)/1000),'km. Estimated Rho: ',num2str(minmod),' \Omega m'])
            if strcmp(err_to_add,'Variable')
                elevel = 'Variable';
            end
            printfile = ['Model_',num2str(k),'_',num2str(elevel),'_error_',num2str(length(mod(1,:))),'runs_1000hs'];
            print('-deps',[printfile,'.eps'])

        end

        %Calculate the estimated melt fraction for all 100 runs
        for L = 1:length(allminmod)
            %To solve for phi using MAL, it must be done numerically. Here
            %the function uses Newton-Raphson
            melt_frac(L) = MAL_solve_phi(10.^allminmod(L),1000,0.61,1.5);
            
            if melt_frac(L)>0.3
                count = count+1;
            end
        end 

        %Plot the results of each of the 100 runs with estimated depth on
        %y-axis and estimated melt fraction on x-axis
        figure(2)
        subplot(2,2,k)
        plot(melt_frac,alldepth/1000,'LineStyle','none','Marker',LT{i},'Color',colorchart(i,:),'MarkerSize',MS(i),'MarkerFaceColor',colorchart(i,:)); hold on
        plot(0.3,5,'.m','MarkerSize',40)
        xlabel('Melt Fraction');
        ylabel('Depth'); axis ij
        axis([0 0.45 2.5 7])
        grid on

        %Collect all the relevant variables into 16 x 16 matrices.
        fullmeanmod(k,i) = minmod;
        fullstdmod(k,i) = stdmod;
        mindepth(k,i) = min(alldepth)/1000;
        maxdepth(k,i) = max(alldepth)/1000;
        minphi(k,i) = min(melt_frac);
        maxphi(k,i) = max(melt_frac);
        fullmeandepth(k,i) = mean(alldepth)/1000;
        fullstddepth(k,i) = std(alldepth)/1000;
        fullmeanphi(k,i) = mean(melt_frac);
        fullstdphi(k,i) = std(melt_frac);
    end
end

for k = 1:4 %Loop through models and error levels to plot mean depth and melt fraction
    subplot(2,2,k)
    for i = [3 2 4 1]
        %c = rand(3,1);
        c = colorchart(i,:);
        plot(fullmeanphi(k,i),fullmeandepth(k,i),'LineStyle','none','Marker',LT{i},'Color',c,'MarkerSize',MS(i)+6,'MarkerFaceColor',c); hold on
        plot([fullmeanphi(k,i)-fullstdphi(k,i) fullmeanphi(k,i)+fullstdphi(k,i)],[fullmeandepth(k,i) fullmeandepth(k,i)],'-k','Color',c,'LineWidth',3);
        plot([fullmeanphi(k,i) fullmeanphi(k,i)],[fullmeandepth(k,i)-fullstddepth(k,i) fullmeandepth(k,i)+fullstddepth(k,i)],'-k','Color',c,'LineWidth',3);
        color(i,:) = c;
    end
    
    grid on; axis ij
    axis([0 0.45 3.5 6])
    title(['Model Geometry #',num2str(k)])
    if k == 4
        manual_legend('Low Noise (0.025)',color(1,:),'Medium Noise (0.05)',color(2,:),'High Noise (0.1)',color(3,:),'Frequency-Dependent Noise',color(4,:));
    end
    
end


