function [ k ] = keffcalonne2019( r )
%keffcalonne2019( r ) returns the effective thermal conductivity of snow/firn (k) in [J (s m K)^-1] 
% as a function of density (r) in [kg m^-3] following:
% Calonne, N., Milliancourt, L., Burr, A.,
% Philip, A., Martin, C. L., Flin, F.,
% & Geindreau, C. (2019). Thermal
% conductivity of snow, firn, and
% porous ice from 3-D image-based
% computations. Geophysical Research
% Letters, 46. https://doi.org/10.1029/
% 2019GL085228
k   = 2.107 + 0.003618.*(r-910);
end

