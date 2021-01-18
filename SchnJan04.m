function [w] = SchnJan04(r, gv)
%SchnJan04
% given density (r) returns the irreducible water content of snow/firn as either
% a gravimetric (gv = "g") or a volumetric (gv = "v")
% fraction following Schneider and Jansson, 2004

i = 1;
n (:,i) = 1 - (r - 0)/910;               % initial guesses of: porosity
w (:,i) = 0.0143.*exp( 3.3 .* n(:,i) );  %                     the gravimetric water content
dw(:,i) = repmat(0.5, length(r), 1);     %                     the differernce in gravimetric water content for neighboring steps

while max(dw) > 0.0001
      i = i + 1;
      n(:,i) = 1 - (r - w(:,i-1).*r )/910;    % update the estimate of porosity
      w(:,i) = 0.0143.*exp(3.3.*n(:,i));      % update the estimate of gravimetric water content 
      dw(:,i) = abs( w(:,i) - w(:,i-1) );
end

w = w(:,end);

if gv == 'v'             % optionally convert to volumetric water content
    w = r.*w./1000;
end

end

