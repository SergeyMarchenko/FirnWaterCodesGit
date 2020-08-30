function [ k ] = keffriche( r )
%keffrische( r ) returns the effective thermal conductivity of snow/firn (k) in [J (s m K)^-1] 
% as a function of density (r) in [kg m^-3] following:
% Riche, F. and Schneebeli, M.: Thermal conductivity of snow measured
% by three independent methods and anisotropy considerations,
% The Cryosphere, 7, 217–227, https://doi.org/10.5194/tc-7-
% 217-2013, 2013.
k   = 3   * 10^-6 * r.^2 - 1.06 .* 10^-5 .* r + 0.024;
end