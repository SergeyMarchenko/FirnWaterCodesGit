function [Ts, Wm] = Tkw(t, z, TIC, R, K, Wm, thr)
% simulate the evolution of subsurface temperature [Ts] based on
% t   - time vector
% z   - depth vector
% TIC - temperature initial conditions
% R   - density profile,
% K   - thermal conductivity profile,
% Wm  - water mass profile,
% thr - temperature of water freezing

Ts = TIC;
up = ones(size(t))+1;

dn = ones(size(t)) .* (length(z)-1);      % lower BC at lowermost depth
Ts( up:dn, 2:end ) = NaN;

Wm( isnan(Wm) ) = 0;
Wm = [Wm nan(size(Ts(2:end,2:end)))];

for it = 2:size(t, 2)

% if rem(it,10) == 0
% disp(['t = ' num2str(it) ' of ' num2str(size( t, 2))] );
% end

dz  =   z( 2  ) - z( 1  );
dt  = ( t( it ) - t(it-1) ) * 24 * 3600;

% new time step to ensure that CFL<0.45
dtn = max(  dt .* K   (    up(it)-1:dn(it)+1     )   ./ ...  % "worst" (=largest) CFL parameter for all fin diff equations
                  R   (    up(it)-1:dn(it)+1     )   ./ ...
                  Cice( Ts(up(it)-1:dn(it)+1,it-1) ) ./ ...
                  dz^2  );
dtn = ceil(  dtn / 0.45  );                                 % number of additional steps in time required
dtn = dt / dtn;                                             % new time step in [seconds]

tmpW = [Wm(:,it-1)                    NaN( dn(it)-up(it)+2, round(dt/dtn  ) )                            ];
tmpT = [Ts(up(it)-1:dn(it)+1,it-1)    NaN( dn(it)-up(it)+3, round(dt/dtn-1) )    Ts(up(it)-1:dn(it)+1,it)];

% tmp = NaN( size(z,1), round(dt/dtn+1) );

tmpT( 1  , : ) = interp1([t(it-1) t(it)], [ tmpT(  1 , 1 )  tmpT(  1 , end )], [t(it-1):dtn/3600/24:t(it)]);
tmpT( end, : ) = interp1([t(it-1) t(it)], [ tmpT( end, 1 )  tmpT( end, end )], [t(it-1):dtn/3600/24:t(it)]);

% tmp ( : , 1 ) = Ts ( : , it-1);

% temperature evolution
for  tn = 2:size(tmpT,2)
for  iz = 2:size(tmpT,1)-1
    T11 = tmpT( iz-1, tn-1 );
    T21 = tmpT( iz  , tn-1 );
    T31 = tmpT( iz+1, tn-1 );
    R2  = R  ( iz         );
    K1  = K  ( iz-1       );
    K2  = K  ( iz         );
    K3  = K  ( iz+1       );
    C   =    Cice(T21 );
tmpT(iz,tn) = T21  + dtn / (2*R2) / C / dz^2 * ( (K3+K2)*(T31  - T21 ) - (K2+K1)*(T21  - T11 ));
end

tmpW(:, tn) = tmpW(:, tn-1);                         % water from the previous time step

% the effect of pore water refreezing on subsurface temperature and water content
ind_cold = find( tmpT(:,tn) <= thr );       % cold points
ind_wet  = find( tmpW(:,tn) >   0  );       % wet layers
ind = intersect( ind_cold, ind_wet );       % indexes of cold points WITH water in the layer just above
if isfinite(ind)
    Ec = abs(tmpT(ind, tn)-thr) .* R(ind) .* Cice( tmpT(ind,tn) ) .* dz; % energy needed to warm up [ind] points to temperature [thr]  (J/m^2)
    Ew =     tmpW(ind, tn)      .* 333500;                             % energy released if all water in [ind] layers refreezes      (J/m^2)

    % update current temperature and water
    tmpT(ind, tn) = min( thr, tmpT(ind, tn) + Ew ./ R(ind) ./ Cice( tmpT(ind,tn) ) ./ dz );
    tmpW(ind, tn) = max( 0  , tmpW(ind, tn) - Ec ./ 333500                               );
end

end
Ts (up(it)-1:dn(it)+1,it) = tmpT(:,end  );        % temperature on the current time step after conductive heat exchange and water refreezing
Wm (        :        ,it) = tmpW(:,end  );        % water content after refreezing
clear tmp* dz dt dtn tn iz T11 T21 T31 R2 K1 K2 K3 C ind* Ec Ew

end

end