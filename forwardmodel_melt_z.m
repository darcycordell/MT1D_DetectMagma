function dpred = forwardmodel_melt_z(freq_array,phi,z)
%Forward model MT data given log10 resistivity and log10 depth
%
% Usage: dpred = forwardmodel_melt_z(freq_array,phi,z)
%
% Inputs: freq_array is a (1 x nf) vector ascending frequencies
%         melt is a (nL x 1) vector of melt fractions between 0 and 1 (nL is number of
%            layers in the model.
%         z is a (nL-1 x 1) vector of log10 depths (the top depth is pinned
%            at zero and is not included as a parameter
%
% Outputs: dpred is a (2*nf x 1) vector of MT impdeances [real(Z) imag(Z)];

m=1.5; %Connectivity parameter
sigf = 1/0.61; %fluid resistivity
sigm = 1/2000; %matrix resistivity

p = (log10(1-phi.^m)./(log10(1-phi)));
res =  1./(sigf.*(phi).^m+ sigm*(1-phi).^p);

z = [0; 10.^z];

fwd = calc_fwd_1d(z,res,freq_array,0);

dpred = [real(fwd.Z) imag(fwd.Z)]';

end