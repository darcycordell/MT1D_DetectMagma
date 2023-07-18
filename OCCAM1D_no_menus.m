function result = OCCAM1D_no_menus(MTdata,plot_fig)
% This code solves the OCCAM 1D Solution for the MT problem. This was
% primarily written by Ersan but re-worked and re-written by Darcy with
% comments. Ersan's GUI was also removed to make for an easier and more
% usable code
%
% User is able to load EDI data and choose to invert either the TE, TM or
% determinant impedance
%
% User is also able to load a text file containing columns of frequency,
% apparent resistivity, phase, apparent resistivity error and phase error
%
% Finally, user is able to make synthetic data within the program and
% invert it.
%
% For more info on OCCAM inverison see:
%   Constable, S., Parker, R. L., & Constable, C. (1987). Occam’s inversion: 
%   A practical algorithm for generating smooth models from electromagnetic 
%   sounding data. Geophysics, 52(3), 289–300. 
%
% (C) 2020 Unsworth Research Group (University of Alberta, Edmonton, Canada)

%close all
% Plotting options
reslims = [1 10000]; %resistivity limits 
depthlims = [0 10000]; %depth limits for plotting

rms_goal = 1;
halfspace_value = 100;
error_floor = 0;


track = 0; trackcount = 1;
%while 1

%main_menu = menu('',get_data,mesh_design,inv_params,inv_run,'Plot Results','Save Results','Quit');

main_menu = 1;
if main_menu == 1
    track(trackcount) = main_menu;
    % DATA INPUT OPTIONS  ----------------------------------------------------
    %data_menu = menu('','Load EDI','Load Text File','Make Synthetic Data (default)');


        %error_floor = MTdata.err;
        f = MTdata.freq_array';
        nd = length(f);
        T = 1./f;
        lat = 0;
        lon = 0;
        elev = 0;
        rot = 0;


        Z = MTdata.Z(1:nd)+1i*MTdata.Z(nd+1:end);
        dZ = MTdata.Zerr(1:nd);
        floorZ = abs(Z);
        
        [rhoa,pha] = calc_rho_pha(Z,dZ,T);
        rhoa = rhoa';
        pha = pha';

    %Set constants and frequencies
    w=2*pi*f';mua=4*pi*10^-7;
end



main_menu = 2;
if main_menu == 2 && max(track) >= 1
    track(trackcount) = main_menu;
    % MESH DESIGN ------------------------------------------------------------
    %The mesh can be designed manually using skin depth, but the way the code
    %is currently implemented has a mesh generated automatically using Bostick
    %inversion. (Bostick is automated and works well)
    %mesh_menu = menu('Mesh Design Using:','Bostick Automated (default)','Manually Set','Load Text File');
    mesh_menu = 2;
    if mesh_menu == 2 %Manually Set Thicknesses
         meshtype = 'Manually Entered';
%         prompt={'Number of Layers','First Thickness in meters (default is 1/3 min skin depth)','Maximum Depth in meters (default is 3*(max skin depth))','Geometric Factor'};
%         dlg_title='Mesh Parameters';
%         %def={num2str(80),num2str(503*sqrt(min(T)*100)/3),num2str(503*sqrt(max(T)*100)*3),num2str(1.2)};
%         def={num2str(137),num2str(10),num2str(200000),num2str(1.1)};
%         num_lines=1;
%         dinp = inputdlg(prompt,dlg_title,num_lines,def);

%         nl = str2double(dinp{1});
%         first_layer = str2double(dinp{2});
%         maxdepth = str2double(dinp{3});
%         geo_factor = str2double(dinp{4});

        first_layer = 3;
        maxdepth = 1000;
        geo_factor = 1.1;

        thick=[]; thick(1) = first_layer; i =2;
        while sum(thick)<=maxdepth
            thick(i)=thick(i-1)*geo_factor;
            i=i+1;
        end
        
        %thick = [thick 200 300 408 500*ones(1,16)];
        thick = [thick 100*ones(1,80)];
        
        maxdepth = 500000;
        geo_factor = 1.2; i = length(thick)+1;
        while sum(thick)<=maxdepth
            thick(i)=thick(i-1)*geo_factor;
            i=i+1;
        end
        
    end

    nl = length(thick)+1; %Number of layers
    depth=[0 cumsum(thick)]; %Depth of layer tops
end

main_menu = 3;
if main_menu == 3 && max(track) >= 2
    % INVERSION PARAMETERS ---------------------------------------------------
    
    dinp = {'200',num2str(error_floor),num2str(rms_goal),num2str(halfspace_value),'0'};
    
    if ~isempty(dinp)
        track(trackcount) = main_menu;
        itmax = str2double(dinp{1}); %Maximum number of iterations (usually 10 - 30 is ok)
        err_flr = str2double(dinp{2}); %The data error floor (usually 0.01 to 0.1 is ok). This is a true error floor as a direct
                %percentage of the apparent resistivity and phase, independently.
                %There is no scaling or error propagation carried out.
                %The choice of error floor has a significant impact on the
                %inversion fit and needs to be tuned accordingly.

        %taumax = 100; %Maximum conductance (for plotting purposes)
        rmstreshold = str2double(dinp{3}); %RMS misfit threshold (usually 1 is good for synthetic)
        inres = str2double(dinp{4}); %Starting model resistivity in Ohm m.
        discontinuity = str2double(dinp{5});
        
        if discontinuity
            prompt={'Depth of Disontinuity or Discontinuities (km)'};
            dlg_title='Discontinuity Parameters';
            def={num2str([1 5])};
            num_lines=1;
            disc_inp = inputdlg(prompt,dlg_title,num_lines,def);
            
            if isempty(disc_inp)
                discontinuity = 0;
            else
                disc_depth = str2num(disc_inp{1});
                for i = 1:length(disc_depth)
                    disc_dist = abs(depth - disc_depth(i)*1000);
                    [~, disc_idx(i)] = min(disc_dist);
                end
            end
                    
        end
        
    end

end

main_menu = 4;
if main_menu == 4 && max(track) >= 3
    track(trackcount) = main_menu;
    inv_run = '(4) Run Inversion (DONE)';
    
    %Remove NaN from data
    ind = isnan(Z);
    Z(ind) = [];
    dZ(ind) = [];
    f(ind) = [];
    floorZ(ind) = [];
    pha(ind) = [];
    rhoa(ind) = [];
    T(ind) = [];
    w(ind) = [];
    nd = length(f);
    
    % OCCAM INVERSION----------------------------------------------------------

    %Apply error floor-------------------------------
    [dZ] = apply_errorfloor(dZ,err_flr,floorZ);

    rhoerr = (2*rhoa'.*dZ)./abs(Z);
    phaerr = (180/pi)*(dZ)./abs(Z);

    W = diag(1./[dZ; dZ]); %Data error weighting matrix used in inversion
    d=[real(Z);imag(Z)]'; %Data vector used in inversion

    %-------------starting model-------------------
    clear m J F
    m(1,:)=ones(1,nl)*inres;

    %FORWARD CALCULATION ---------------------------------------------
    %"F" is the forward of real and imaginary impedance values (across columns). Each
    %inversion iteration is saved in a new row.
    [F(1,:),ra,ph]=mt_fwd_occ(m(1,:),nl,w,thick); 
    %-----------------------------------------------------------------

    if plot_fig
        %Plot initial data and starting model--------------------------------------
        set_figure_size(1);
        subplot(2,2,1)
        logerrorbar(T,rhoa',rhoerr,'.k','-k'); hold on
        loglog(T,ra','--k');
        xlabel('Period (s)')
        ylabel('App Res (\Omega m)')
        axis([min(T) max(T) reslims])
        grid on

        subplot(2,2,3)
        errorbar(T,pha',phaerr,'.k'); hold on
        semilogx(T,ph','--k')
        set(gca,'XScale','log')
        xlabel('Period (s)')
        ylabel('Phase (deg)')
        axis([min(T) max(T) 0 90])
        grid on

        subplot(2,2,[2 4])
        stairs(m(1,:),depth/1000,'-k'); hold on;
    end

    %--------------------------------------------------------------------------

    %Begin inversion
    X2 = []; rms = []; R1= []; D1 = [];
    rms(1) = norm(W*(F(1,:))'-d')/sqrt(2*nd); %Initial rms of starting model
    rgh1=diag([0;ones(nl-1,1)]) + diag(-ones(nl-1,1),-1);         %delta matrix
    
    if discontinuity
        rgh1(disc_idx,:) = 0;
    end
    
    rgh2=rgh1'*rgh1;                                              %delta' x delta  
    mumax = []; mumax(1)=10000;
    iter = 2; 

    %Iterate to solve the inversion problem using Occam algorithm
    for iter=2:itmax

        if rms(iter-1)>2
            rmsdes=rms(iter-1)/1.5;
        else
            rmsdes=rmstreshold;
        end

        %Build the Jacobian matrix --------------------------------------------
        %I am not sure why Ersan chose to do this in logarithmic space but it
        %seems strange to me that he takes logarithms of both apparent
        %resistivity AND phase (since phase is a linear quantity). I modified
        %the code so that it works by inverted real and imaginary impedance
        %data instead of rho and pha.
        parameter = log10(m(iter-1,:));
        apt = (F(iter-1,:));
        for I=1:length(parameter);
            parameter(I)=parameter(I)+0.005;
            [cpt,~,~] = mt_fwd_occ(10.^parameter,nl,w,thick);      
            turev=(cpt-apt)/0.005;
            J(:,I)=turev';
            parameter(I)=parameter(I)-0.005;
        end

        %----------------------------------------------------------------------

        %Inverse algorithm (from Constable et al., 1987)
        son=(W*J)'*W*J;
        b=(W*J)'*W*(d-(F(iter-1,:))+(J*log10(m(iter-1,:))')')';
        %The solution is found using a "golden section search" algorithm to
        %find the minimum
        [mumax(iter), rms(iter), m(iter,:), F(iter,:),x3,f4]=golden_section(son,b,W,rgh2,d,rmsdes,0.0001,1,mumax(1), 0.01,2*nd,w,thick,nl);
        X3(:,iter-1)=x3';  F4(:,iter-1)=f4';
        R1(iter)=((rgh1*log10(m(iter,:))')'*rgh1*log10(m(iter,:))');                 
        DD(iter)=(log10(m(iter,:))-log10(m(iter-1,:)))*(log10(m(iter,:))-log10(m(iter-1,:)))';   
        if DD(iter)<0.01
            break
        end

        if abs(rms(iter)-rms(iter-1)) < 0.01
            break
        end

    end

    X2=rms.^2*nd;
    [~,ra_mod,ph_mod]=mt_fwd_occ(m(end,:),nl,w,thick); 

    if plot_fig
        %PLOT RESULTS--------------------------------------------------------------
        figure(1)
        subplot(2,2,1)
        h(1) = loglog(T,ra_mod,'-m','LineWidth',1);
        title(['Total Iterations = ',num2str(iter),'. Final Occam RMS = ',num2str(rms(end))])

        figure(1)
        subplot(2,2,3)
        h(2) = semilogx(T,ph_mod,'-m','LineWidth',1);

        figure(1)
        subplot(2,2,[2 4])
        h(3) = semilogx(m(iter,:),depth/1000,'-m','LineWidth',1);
        set(gca,'XScale','log')
        axis([reslims depthlims/1000])
        xlabel('Resistivity (\Omega m)')
        ylabel('Depth (km)')
        axis ij
        grid on
    end
    
end
% 
% if main_menu == 5 && max(track) >= 4
%     track(trackcount) = main_menu;
%     %PLOT ALL ITERATIONS--------------------------------------------------
%     im = 1;
%     while 1
% 
%         [~,ra_view,ph_view]=mt_fwd_occ(m(im,:),nl,w,thick); 
% 
%         h(1).YData = ra_view;
%         subplot(2,2,1)
%         title(['Total Iterations = ',num2str(im),'. Final Occam RMS = ',num2str(rms(im))])
%         h(2).YData = ph_view;
% 
%         figure(1)
%         subplot(2,2,[2 4])
%         stairs(m(im,:),depth/1000,'-m','LineWidth',1); hold off
%         set(gca,'XScale','log')
%         axis([reslims depthlims/1000])
%         xlabel('Resistivity (\Omega m)')
%         ylabel('Depth (km)')
%         axis ij
%         grid on
% 
%         next_menu = menu('','Next Iteration','Quit');
% 
%         if next_menu == 1
% 
%             im = im+1;
%             if im > iter
%                 break
%             end
% 
%         else
%             break
%         end
% 
%     end
% 
% end
% 
rms_mod = rms(end);
model = m(end,:)';
starting_model = m(1,:)';

result.T = T;
result.ra_mod = ra_mod;
result.ph_mod = ph_mod; result.rms_mod = rms_mod;
result.model = model; result.thick = thick; result.depth = depth;
result.iter = iter; result.lat = lat; result.lon = lon; result.elev = elev;
result.rot = rot; result.datatype = 'Bayesian'; result.meshtype = meshtype;
result.err_flr = err_flr; result.starting_model = starting_model;
result.rhoerr = rhoerr; result.phaerr = phaerr; result.rhoa = rhoa;
result.pha = pha; result.Z = Z; result.dZ = dZ;
% 
% if main_menu == 6 && max(track) >= 4
%     track(trackcount) = main_menu;
%     %SAVE RESULTS
% %     prompt={'Filename to Save'};
% %     dlg_title='Filename to Save';
% %     def={'filename_occam1d'};
% %     num_lines=1;
% %     dinp = inputdlg(prompt,dlg_title,num_lines,def);
%     
%     %filename = dinp{1};
%     
%     rms_mod = rms(end);
%     model = m(end,:)';
%     starting_model = m(1,:)';
%     save([filename,'.mat'],'T','ra_mod','ph_mod','rms_mod','model','thick','depth','iter','lat','lon','elev','rot','datatype','meshtype','err_flr','err_to_add','starting_model');
%     save([filename,'_data_input.mat'],'rhoerr', 'phaerr', 'rhoa', 'pha', 'Z', 'dZ');
%     
%     fid = fopen([filename,'_ef_',num2str(err_flr*100),'.mod'],'w');
%     A = [thick' model(1:end-1)];
%     fprintf(fid,'%12.3f %e\r\n',A');
%     fclose(fid);
%     
%     
%     fid = fopen([filename,'_ef_',num2str(err_flr*100),'.resp'],'w');
%     write_error = 10^-9*ones(1,nd);
%     A = [f ra_mod' ph_mod' write_error' write_error'];
%     fprintf(fid,'%e %e %e %e %e\r\n',A');
%     fclose(fid);
%     
% end
% 
% if main_menu > 6
%     %break
% end
% 
% trackcount = trackcount+1;
%     
% end
    
    
    





