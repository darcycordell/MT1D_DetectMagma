function [sigf] = MAL_solve_sigf(sigm,sigb,m,phi)
% Solve Modified Archie's Law (Glover et al., 2000) for melt resistivity
%
% Usage: sigf = MAL_solve_sigf(sigm, sigb, m, phi)
%
% Inputs:
%       phi = volume fraction of more conductive material
%       sigm = matrix resistivity
%       sigf = fluid resistivity
%       m = connectivity parameters
%
%
if sigm>=sigb
    error('The volume fraction must be associated with the more conductive material. Try (1-phi) instead')
end

p = (log10(1-phi.^m))./(log10(1-phi));
sigf = (sigb - sigm.*(1-phi).^p)./(phi.^m);


end