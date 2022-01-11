function [sig] = MAL(sigm,sigf,m,phi)
% Solve Modified Archie's Law (Glover et al., 2000)
%
% Usage: sig = archie(phi, sigm, sigf, m)
%
% Inputs:
%       phi = volume fraction of more conductive material
%       sigm = matrix resistivity
%       sigf = fluid resistivity
%       m = connectivity parameters
%
%
if sigf<=sigm
    sigf = sigm;
    disp('The volume fraction must be associated with the more conductive material. Try (1-phi) instead')
end

p = (log10(1-phi.^m))./(log10(1-phi));
sig = sigf.*(phi).^m+ sigm.*(1-phi).^p;


end