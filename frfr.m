
function [ Zctt ] = frfr( data, Z, thr )
%FRFR(data, Z, thr) find depths of the freezing front
%   given:
%          temperature matrix (depth changes along columns, time along rows) (data),
%          depth vector                                                      (Z),
%          a definition of the freezing front in terms of temperature        (thr)
%   the function returns the depth of freezing front (Zctt) for every time step:
%         linearly interpolated depth of the uppermost point where temperature reaches the threshold value
    for i = 1:size(data,2)                                          % loop over data columns (moments in time)
        ind = find( data(:,i) < thr );                              % indexes of layers with temperature < threshold
        if     isempty(ind);                      Zctt(i) = NaN;    % return NaN if: everything is temperate
        elseif min(ind) ~= 1                      Zctt(i) = NaN;    %                the uppermost layer is temperate
        elseif length(ind) == 1;                  Zctt(i) = interp1(data(ind:ind+1,i), Z(ind:ind+1), thr);
        else
               if    max( diff(ind) ) == 1;  ind = max( ind);
               else;                         ind = min( find(diff(ind)>1) );
               end

               if   ind == size(Z,1); Zctt(i) = Z(end);
               else                   Zctt(i) = interp1(data(ind:ind+1,i), Z(ind:ind+1), thr);
               end
        end
    end

end

