function [sig] = MAL3(sigc,sigm,sigg,phic,phim,phig)
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
m = 1;

Gc = phic.^m;
Gm = phim.^m;
Gg = phig.^m;

Gchecksum = Gc + Gm + Gg;

phichecksum = phic + phim + phig;

sig = sigc.*Gc + sigm.*Gm + sigg.*Gg;

if abs(Gchecksum-1) > 10^-3
    disp('G not equal 1')
end

if abs(phichecksum-1) > 10^-3
    disp('phi not equal 1')
end





end