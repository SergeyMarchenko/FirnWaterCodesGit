function [ k ] = keffcalonne( r )
%keffsturm( r ) returns the effective thermal conductivity of snow/firn (k) in [J (s m K)^-1] 
% as a function of density (r) in [kg m^-3] following:
% Calonne N, Milliancourt L, Burr A, Philip A, Martin CL, Flin F and Geindreau C (2019)
% Thermal conductivity of snow, firn, and porous ice from 3-d image-based computations. Geophysical Research Letters, 46(22), 13079-13089
% (doi: 10.1029/2019GL085228)

k = 0.024 - 1.23*10^-4.*r + 2.5*10^-6.*r.^2;
end