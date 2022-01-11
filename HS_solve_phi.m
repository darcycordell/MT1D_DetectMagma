function [phi] = HS_solve_phi(rho_b,rho_h,rho_f)
% Function which solves for melt fraction using HS+
% given a known bulk resistivity, matrix resistivity, and fluid resistivity
%
% Usage: [phi,i] = HS_solve_phi(rho_b,rho_m,rho_f)
%
%   Inputs: rho_b = bulk resistivity (Ohm m)
%           rho_h = host rock resistivity (Ohm m)
%           rho_f = fluid resistivity (Ohm m)
%
%   Outputs: phi = melt fraction (or porosity)
%
% HS+ formula:
%
% sig_b = sig_f + (1-phi)/((1/(sigm-sigf))+phi/(3*sigf))
%
%
% The formula cannot be solved algebraically so the Newton-Raphson method
% is used to solve the equation numerically. Since 0<phi<1, a hard-coded
% starting guess of 0.5 is used.

% Check for input errors first:
flag = 0;
if rho_b > rho_h
    disp('Bulk resistivity cannot be greater than matrix resistivity')
    flag = 1;
end

if rho_b < rho_f
    disp('Bulk resistivity cannot be less than fluid resistivity')
    flag = 1;
end

if rho_f > rho_h
    disp('Fluid resistivity cannot be greater than matrix resistivity')
    flag = 1;
end

sigf = 1./rho_f;
sigm = 1./rho_h;
sigb = 1./rho_b;


phi = (3*sigf.*(sigm-sigb))./(3*sigf.*(sigm-sigf+(1/3)*sigf)+sigm.*sigb-sigf.*sigm-sigb.*sigf);

