function [sigplus,sigminus] = HS(sigm,sigf,phi)
% Solve Hashin Shtrickman Upper and Lower Bounds
%
% Usage: sig = HS(sigm, sigf, m, phi)
%
% Inputs:
%       phi = volume fraction of more conductive material
%       sigm = matrix resistivity
%       sigf = fluid resistivity
%
%
if sigf<=sigm
    error('The volume fraction must be associated with the more conductive material. Try (1-phi) instead')
end

phi1 = phi;
phi2 = 1-phi1;

sigplus = sigf + phi2./((1./(sigm-sigf))+(phi1./(3*sigf)));
sigminus = sigm + phi1./((1./(sigm+sigf))+(phi2./(3*sigm)));


%Alternate form of HS+ from Glover 2000
%sigplus = (sigf).*(1-(3*(1-phi).*(sigf-sigm))./(3*sigf-phi*(sigf-sigm)));

end