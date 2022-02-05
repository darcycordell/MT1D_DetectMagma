function [sig] = MAL3(sigc,sigm,sigg,phic,phim,phig,mm,mg)
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
if phic == 0
    phic = 10^-4;
end

Gm = phim.^mm;
Gg = phig.^mg;

mc = log(1-Gg-Gm)./log(phic);

Gc = phic.^mc;

Gchecksum = Gc + Gm + Gg;

sig = sigc.*Gc + sigm.*Gm + sigg.*Gg;

if abs(Gchecksum-1) > 10^-3
   error('G not equal 1. Something wrong.')
end

if imag(1/sig)>10^-3
    error('Imaginary conductivies. Something wrong.')
else
    sig = real(sig);
end
    

%     phichecksum = phic + phim + phig;
%     if abs(phichecksum-1) > 10^-3
%         disp('phi not equal 1')
%     end






end