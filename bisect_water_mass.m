function  [m,d] = bisect_water_mass(f, L, H, tol)


% L = -2; % H = 8; % lim = 1;
% function [y] = fun( x, k, b, lim )
% y = k.*x + b; y(y<lim) = lim;
% end
% f = @(x)fun( x, -1, 5, 1 );
% x = -2:8; k = -1; b = 5; lim = 1; y = fun(x,k,b, lim); close all; figure; plot(x,y, 'o')
% [L m H; y_L y_m y_H]

i = 0;
y_L = feval(f, L);
y_H = feval(f, H);
d = nan(1000, 8);
if      y_L       < 0;      m = 0;  y_m = feval(f, m);    i = i+1;    d(i,:) = [i  L  y_L  m  y_m  H  y_H 1];                      return

elseif  y_H       > 0;      while (abs(H - L) >= tol)
                                      m = (H + L)/2;
                                    y_m = feval(f, m);    i = i+1;    d(i,:) = [i  L  y_L  m  y_m  H  y_H 2];

                                    if          y_m == 0 ;      return
                                    elseif      y_m > y_H;      L = m;  y_L = y_m;
                                    else;                       H = m;  y_H = y_m;
                                    end
                            end

elseif  y_L * y_H < 0;      while (abs(H - L) >= tol)
                                      m = (H + L)/2;
                                    y_m = feval(f, m);    i = i+1;    d(i,:) = [i  L  y_L  m  y_m  H  y_H 3];
                                    if          y_m == 0;       return
                                    elseif y_L * y_m < 0;       H = m;  y_H = y_m;
                                    else                        L = m;  y_L = y_m;
                                    end
                            end
end

d( i+1:end,: ) = [];

end
