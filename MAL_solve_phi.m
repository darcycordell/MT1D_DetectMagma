function [phi,i] = MAL_solve_phi(rho_b,rho_h,rho_f,m)
% Function which solves for melt fraction using Modified Archie's Law
% given a known bulk resistivity, matrix resistivity, fluid resistivity,
% and connectivity exponent.
%
% Usage: [phi,i] = solve_phi_MAL(rho_b,rho_m,rho_f,m)
%
%   Inputs: rho_b = bulk resistivity (Ohm m)
%           rho_h = host rock resistivity (Ohm m)
%           rho_f = fluid resistivity (Ohm m)
%           m = connectivity parameter (m>0, unitless)
%
%   Outputs: phi = melt fraction (or porosity)
%            i = number of iterations to reach solution
%
% Modified Archie's Law (MAL) formula:
%
%   rho_b*rho_h*phi^m + rho_b*rho_f(1-phi)^p = rho_f*rho_h
%   
%       p = log10(1-phi^m)/log10(1-phi)
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


N_max = 1000; %Maximum # of iterations (rarely need >10)
epsilon = 10^-5; %Stopping criterion (rarely need melt fractions with greater precision than this)
x = zeros(N_max,1); %Initialize unknown
x(1) = 0.25; %Starting guess is 0.5 just because 0 < phi < 1

i = 2;
while 1   
    p = log10(1-x(i-1).^m)./log10(1-x(i-1)); %set p
    f = rho_f*rho_h*(rho_h*x(i-1).^m+rho_f*(1-x(i-1)).^p).^-1-rho_b; %MAL
    %Derivate of f is ridiculously complicated because of the "p" (e.g. phi to an
    %exponent containing phi).
    da = m*rho_h*x(i-1).^(m-1);
    db = rho_f*(1-x(i-1)).^p.*((log10(rho_f)*log10(1-x(i-1).^m).^2-m*x(i-1).^(m-1).*log10(rho_f).*log10(1-x(i-1)).^2-m*x(i-1).^(m-1).*log10(1-x(i-1)).^3)./(log10(1-x(i-1).^m).*log10(1-x(i-1)).^3));
    fp = -rho_f*rho_h*(rho_h*x(i-1).^m+rho_f*(1-x(i-1)).^p).^-2.*(da+db);


    x(i) = x(i-1) - f/fp; %Newton-Raphson
    
    if x(i)>=1 %If the algorithm finds a resistivity >1 then this is non-physical so maybe try a different starting guess
        x(i) = 0.99;
    end
    
    if x(i)<=0 %Similarly if the algorithm finds a resistivity <0 then try a different starting guess
        x(i) = 0.01;
    end

    %Stopping criteria
    if abs(x(i)-x(i-1))<epsilon
        break
    end

    if i>1000
        break
    end

    %x(i)
    
    i = i+1;
    
    

end

phi = x(i);

if flag == 1
    phi = NaN; i = NaN;
end