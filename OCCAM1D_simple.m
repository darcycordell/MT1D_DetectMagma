function result = OCCAM1D_simple(input_model,resistivity,err,halfspace_value,rms_goal)
% This code solves the OCCAM 1D Solution for the MT problem using synthetic data.
%
% OCCAM1D_simple takes inputs and has no user GUI or menus. The model mesh is 
% generated automatically using set defaults. Number of frequencies is set to default (80).
% (See OCCAM1D.m for a more user-friendly interface with more options.)
%
% Usage: result = OCCAM1D_simple(input_model,resistivity,err_to_add,halfspace_value,error_floor,rms_goal)
%
% Inputs:
%       input_model: a vector of depths (in meters) to layer boundaries for the true model from which synthetic 
%                   data will be calculated
%       resistivity: a vector of resistivity values (Ohm m) for each layer
%       err: The error fraction (e.g. 0.05) of the modulus of the impedance to add to the data. 
%                   Note if NaN then "frequency-dependent" noise is added
%       halfspace_value: Starting halfspace resistivity (e.g. 1000 Ohm m)
%       rms_goal: Target r.m.s. (usually 1)
%
% Outputs:
%       result: A structure containing all the relevant variables including
%           Periods, modelled rho/pha, final r.m.s., final model resistivity
%           values, layer depths, and layer thicknesses, # of iterations
%           for convergence, starting model, synthetic input rho/pha/Z, and
%           errors used on Z.
%
%
% For more info on OCCAM inverison see:
%   Constable, S., Parker, R. L., & Constable, C. (1987). Occam’s inversion: 
%   A practical algorithm for generating smooth models from electromagnetic 
%   sounding data. Geophysics, 52(3), 289–300. 
%
% (C) 2020 Unsworth Research Group (University of Alberta, Edmonton, Canada)

close all

%DATA PARAMETERS-----------------------------------------------------------
datatype = 'SYNTHETIC';

num_freq = 80;
min_freq = 0.001;
max_freq = 1000;

model_depth = input_model; %Depth to top of each model layer
model_res = resistivity; %Resistivity of each layer

f = 10.^linspace(log10(min_freq),log10(max_freq),num_freq).';
nd = length(f);
T = 1./f;

[fwd]=calc_fwd_1d(model_depth,model_res,f',err);

Z = fwd.Z.';
rhoa = fwd.rho;
pha = fwd.phi;

dZ = fwd.Zerr';
floorZ = abs(Z);

%Set constants and frequencies
w=2*pi*f';mua=4*pi*10^-7; iwm=1i*mua*w;  
ndR=2*nd;


%MESH PARAMETERS-----------------------------------------------------------
meshtype = 'Bostick';

first_layer = 3;
maxdepth = 1000;
geo_factor = 1.1;

thick=[]; thick(1) = first_layer; i =2;
while sum(thick)<=maxdepth
    thick(i)=thick(i-1)*geo_factor;
    i=i+1;
end

thick = [thick 200 300 408 500*ones(1,16)];

maxdepth = 200000;
geo_factor = 1.2; i = length(thick)+1;
while sum(thick)<=maxdepth
    thick(i)=thick(i-1)*geo_factor;
    i=i+1;
end

nl = length(thick)+1; %Number of layers
depth=[0 cumsum(thick)]; %Depth of layer tops

%INVERSION PARAMETERS ---------------------------------------------------
itmax = 200; %Maximum number of iterations (usually 10 - 30 is ok)
rmstreshold = rms_goal; %RMS misfit threshold (usually 1 is good for synthetic)
inres = halfspace_value; %Starting model resistivity in Ohm m.

  
%RUN INVERSION-------------------------------------------------------------

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
rhoerr = (2*rhoa'.*dZ)./abs(Z);
phaerr = (180/pi)*(dZ)./abs(Z);

W = diag(1./[dZ; dZ]); %Data error weighting matrix used in inversion
d=[real(Z);imag(Z)]'; %Data vector used in inversion

%-------------starting model-------------------
clear m J F
m(1,:)=ones(1,nl)*inres;

%forward calculation------------------------
%"F" is the forward of real and imaginary impedance values (across columns). Each
%inversion iteration is saved in a new row.
[F(1,:),~,~]=mt_fwd_occ(m(1,:),nl,w,thick); 

%-------------------------------------------
%Begin inversion
X2 = []; rms = []; R1= []; D1 = [];
rms(1) = norm(W*(F(1,:))'-d')/sqrt(2*nd); %Initial rms of starting model
rgh1=diag([0;ones(nl-1,1)]) + diag(-ones(nl-1,1),-1);         %delta matrix

rgh2=rgh1'*rgh1;                                              %delta' x delta  
mumax = []; mumax(1)=10000;

%Iterate to solve the inversion problem using Occam algorithm
for iter=2:itmax

    if rms(iter-1)>2
        rmsdes=rms(iter-1)/1.5;
    else
        rmsdes=rmstreshold;
    end

    %Build the Jacobian matrix --------------------------------------------
    parameter = log10(m(iter-1,:));
    apt = (F(iter-1,:));
    for I=1:length(parameter)
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


rms_mod = rms(end);
model = m(end,:)';
starting_model = m(1,:)';

%Collect variables into data structure
result.T = T;
result.ra_mod = ra_mod;
result.ph_mod = ph_mod; result.rms_mod = rms_mod;
result.model = model; result.thick = thick; result.depth = depth;
result.iter = iter;
result.datatype = datatype; result.meshtype = meshtype;
result.starting_model = starting_model;
result.rhoerr = rhoerr; result.phaerr = phaerr; result.rhoa = rhoa;
result.pha = pha; result.Z = Z; result.dZ = dZ;

end %END MAIN


function [xmin, fmmin, m, F,x3,f4] = golden_section(son,b,W,rgh2,d,rmsdes,ax, bx, cx, tol,nd,w,thick,nl)
% Function which computes golden section algorithm to find minimum
%
% Written by Ersan Turkoglu, circa early 2000s
%
% Updated to work for MT impedances by Darcy Cordell, 2019
b=b';

ndG=nd;

C = (3-sqrt(5))/2;  R = 1-C;
x0 = ax;  x3 = cx;
if (abs(cx-bx) > abs(bx-ax))
  x1 = bx;  x2 = bx + C*(cx-bx);
else
  x2 = bx;  x1 = bx - C*(bx-ax);
end
A1=x1*rgh2+son;
m1=b/A1;
[F1,~,~]=mt_fwd_occ(10.^m1,nl,w,thick);
f1 = norm(W*F1'-W*d')/sqrt(ndG);
A2=x2*rgh2+son;
m2=b/A2;
[F2,~,~]=mt_fwd_occ(10.^m2,nl,w,thick);
f2 = norm(W*F2'-W*d')/sqrt(ndG);
k = 1;
while abs(x3-x0) > tol*(abs(x1)+abs(x2))
if f2 < f1
    x0 = x1;
    x1 = x2;
    x2 = R*x1 + C*x3;   
    f1 = f2;
    A2=x2*rgh2+son;
    m2=b/A2;
    [F2,~,~]=mt_fwd_occ(10.^m2,nl,w,thick);
    f2 = norm(W*F2'-W*d')/sqrt(ndG);
else
    x3 = x2;
    x2 = x1;
    x1 = R*x2 + C*x0;   
    f2 = f1;
    A1=x1*rgh2+son;
    m1=b/A1;
    [F1,~,~]=mt_fwd_occ(10.^m1,nl,w,thick);
    f1 = norm(W*F1'-W*d')/sqrt(ndG);
end
%xx1(k)=x1; ff1(k)=f1; xx2(k)=x2; ff2(k)=f2;  
k = k+1;
end
if f1 < f2
  xmin = x1;
  fmmin = f1;
  F=F1;
  m=10.^m1;
else
  xmin = x2;
  fmmin = f2;
  F=F2;
  m=10.^m2;
end
if fmmin<rmsdes
    xB=logspace(log10(xmin),log10(cx),20);
    for B=1:length(xB)
         AB=xB(B)*rgh2+son;
         mB=b/AB;
         [FB,~,~]=mt_fwd_occ(10.^mB,nl,w,thick);
         fB(B) = norm(W*FB'-W*d')/sqrt(ndG);
     end
    xx=logspace(log10(xmin),log10(cx),500);
    yy=spline(xB,fB,xx);

    xinds=max(find(yy<rmsdes));
    [~, kk]=size(xinds);
    if kk==0
        op=find(min(yy));
        xmin=xx(op);
        A=xmin*rgh2+son;
        m=b/A;
        [F,~,~]=mt_fwd_occ(10.^m,nl,w,thick);
        fmmin = norm(W*F'-W*d')/sqrt(ndG);
        m=10.^m;
    else
        xmin=xx(xinds);
        A=xmin*rgh2+son;
        m=b/A;
        [F,~,~]=mt_fwd_occ(10.^m,nl,w,thick);
        fmmin = norm(W*F'-W*d')/sqrt(ndG);
        m=10.^m;
    end
end

 x3=logspace(-5,4,30);
 for ii=1:30
     A3=x3(ii)*rgh2+son;
     m3=b/A3;
     [F3,~,~]=mt_fwd_occ(10.^m3,nl,w,thick);
     f4(ii) = norm(W*F3'-W*d')/sqrt(ndG);
 end

end %END golden_section

function [F,ra,ph]=mt_fwd_occ(pr,nl,w,thick)
%This function computes the MT response of HLE for Occam's inversion
%Ersan Turkoglu, 2004------eturk@phys.ualberta.ca
%f=frequency, pr=[res1...resn,t1...tn-1(meters)], ra=computed apparent resistivity, ph=computed phase.
%# of layers angular frequency permittivity and iwm are global variables
%nl=number of layers
%w=2*pi*f;                                 %Angular frequency
mua=4*pi*10^-7;                           %Magnetic permeability
iwm=1i*mua*w; 
pr=([1./pr(1:nl),thick]);
Q(nl,:)=ones(1,length(w));
for n=nl-1:-1:1 
    Q(n,:)=(Q(n+1,:)+sqrt(pr(n+1)/pr(n))*tanh(sqrt(iwm*pr(n))*pr(n+nl)))./(sqrt(pr(n+1)/pr(n))+Q(n+1,:).*tanh(sqrt(iwm*pr(n))*pr(n+nl)));
end
Z=sqrt(iwm/pr(1)).*Q(1,:);                %Empedance
ra=1./(w*mua).*abs(Z).^2;                 %Apparent Resistivity
ph=atan(imag(Z(1,:))./real(Z(1,:)))*180/pi;        %phase
% [f',ra',ph']

F=[real(Z)';imag(Z)']';

end %END mt_fwd_occ
    
    
    





