%% calculation of the water content using the forward model
%  as the difference in increment of the refreezing capacity on every time step between measurements and simulation
clearvars
clc; close all;

% Assign parameter values
p  = 0;   % to plot or not to plot, that is the question
load('W_data.mat');

%% Temperature profile: simulated from measurements at the previous time steps accounting for conduction only  -  measured
out.thr = thr;
for is = 1:11
    
if p; figure('units','normalized','outerposition',[0.5 0 0.5 1]); hold on; grid on; set(gca, 'YDir', 'reverse'); c = lines(2); end

out.t {is} = T.t{is};
out.z {is} = T.z{is};
out.T {is} = T.T{is};
out.F0{is} = T.F0{is}(1,:); out.F0{is}(2,1) = out.F0{is}(1,1);
out.K {is} = T.or_kR_f{1*(is == 1) + 2*(is ~= 1)}(:,2);
out.R {is} = T.or_kR_f{1*(is == 1) + 2*(is ~= 1)}(:,3);
W  = zeros(size(out.z{is})); W(end) = [];
out.dT{is}  = nan(size(T.T{is}));
out.dT{is}(:,1) = 0;
for it = 2:size(out.t{is}, 2)

disp(['string = ' num2str(is) ', t = ' num2str(it) ' of ' num2str(size( out.t{is}, 2))] )

Tm  = T.T{is}(:,it     );
Ts  = T.T{is}(:,it-1:it);
Ts(2:end-1, 2:end ) = NaN;

[Ts, ~] = Tkw( out.t{is}(it-1:it), out.z{is}, Ts, out.R{is}, out.K{is}, W, thr );

out.F0{is}(2,it) = frfr(Ts(:,2), out.z{is}, thr);    % depth of ctt in simulated data
out.dT{is}(:,it) = abs( Ts(:,2) ) - abs(Tm); % difference of absolute temperatures

if p;    plot( Ts(:,1)         , out.z{is}                                , '.-' , 'Color', 'k'   , 'DisplayName', 'current IC', 'LineWidth', 1)
         plot( Tm              , out.z{is}                                , '.-' , 'Color', c(1,:), 'DisplayName', 'measured profile')
         plot( Ts(:,2)         , out.z{is}                                , '.-' , 'Color', c(2,:), 'DisplayName', 'model');
         plot( out.dT{is}(:,it), out.z{is}                                , '.-' , 'Color', 'g'   , 'DisplayName', 'model - measurements');
         plot([-1 0]           , [ out.F0{is}(1,it-1) out.F0{is}(1,it-1) ]       , 'Color', 'k'   , 'DisplayName', 'ctt IC')
         plot([-1 0]           , [ out.F0{is}(1,it  ) out.F0{is}(1,it  ) ]       , 'Color', c(1,:), 'DisplayName', 'ctt measurements')
         plot([-1 0]           , [ out.F0{is}(2,it  ) out.F0{is}(2,it  ) ]       , 'Color', c(2,:), 'DisplayName', 'ctt model')
         plot([thr thr]        , [ out.z{is}(1) out.z{is}(end) ], '--' , 'Color', 'k'   , 'DisplayName', 'assumed freezing point')
         xlim([-0.1 0.2]); ylim([1 13.4]); legend show
         title( datestr(out.t{is}(it), 'yy-mmm-dd, hh') );
         pause(0.2); cla;
end

end; clear it
end; clear is c Tm Ts W

%% find the depth indexes of upper and lower boundaries of the positive dT anomaly in vicinity or just below the CTT
m1 = 4;     % max half-thickness of                  the positive dT anomaly (given in number of steps along the depth grid)
m2 = 5;     % max distance from the simulated CTT to the positive dT anomaly (given in number of steps along the depth grid)
p = 0;
for is = 1:11
    
    an_up = nan(1, size(out.t{is},2));
    an_dn = nan(1, size(out.t{is},2));
    
for it = 1:size(out.dT{is},2)

           dT = out.dT{is}(:,it);                                   % |Tsim| - |Tmeas| profile                    at time it
           [~, ind] = min( abs( out.z{is} - out.F0{is}(2,it) ) );   % index of the layer closest to CTT simulated at time it

           if   dT(ind) > 0                                % there is a positive dT anomaly at the simulated CTT
                
                for iz = ind  :-1:max(ind-m1, 1         )            % find the uppermost index of the positive dT anomaly
                    if dT(iz) > 0; an_up(it) = iz;
                    else;          break
                    end
                end

                for iz = ind  : 1:min(ind+m1, size(dT,1))            % find the lowermost index of the positive dT anomaly
                    if dT(iz) > 0; an_dn(it) = iz;
                    else;          break
                    end
                end

           else                                            % no         positive dT anomaly at the simulated CTT
               
                for iz = ind+1: 1:min(ind+m2, size(dT,1))            % try to find the uppermost index of the positive dT anomaly below the CTT
                    if dT(iz) > 0; an_up(it) = iz; break;
                    end
                end
                
                if isfinite(an_up(it))                           % find the lowermost index of the positive dT anomaly
                    for iz = an_up(it)  : 1:min(an_up(it)+2*m1, size(dT,1))   
                        if dT(iz) > 0; an_dn(it) = iz;
                        else;          break
                        end
                    end
                end

           end
           
           if p
           hold on
           plot([1 size(dT,1)], [0 0], 'k')
           plot(dT, '.-')
           if dT(ind) > 0;     c = 'r';    else;    c = 'b';    end
           plot(ind      , dT(ind      ), 'o', 'Color',  c, 'LineWidth', 2)
           plot(an_up(it), dT(an_up(it)), 'o', 'Color', 'g'               )
           plot(an_dn(it), dT(an_dn(it)), 'o', 'Color', 'm'               )
           title(['is = ' num2str(is) ', it = ' num2str(it) ', up = ' num2str(an_up(it)) ', dn = ' num2str(an_dn(it))])
           ylim([min(out.dT{is}(:)) max(out.dT{is}(:))])
           pause;
           if or(isnan(an_up(it)), isnan(an_dn(it)));  pause; end
           clear c; cla;
           end
end

ind1 = find( and( an_up == 1                , an_dn == 1                 ) );
ind2 = find( and( an_up == size(out.z{is},1), an_dn == size(out.z{is},1) ) );
an_up( [ind1 ind2] ) = NaN;
an_dn( [ind1 ind2] ) = NaN;

out.an_ind{is} = find(isfinite(an_up));

an_up( isnan(an_up) ) = [];
an_dn( isnan(an_dn) ) = [];

out.an_up{is} = an_up;
out.an_dn{is} = an_dn;

out.anz_up{is} = out.z{is}(an_up)';
out.anz_dn{is} = out.z{is}(an_dn)';

end; clear m1 m2 is it dT ind iz an_*


%%
close all;
% hold on
for is = 1:11
       figure('units','normalized','outerposition',[0 0 1 1]);
       
       ax1 = subplot(2,1,1); hold on;
       imagesc(out.t{is}, out.z{is}, out.dT{is});
       plot(out.t{is}                , out.F0{is}(2,:), 'k')
       plot(out.t{is}(out.an_ind{is}), out.anz_up{is}, 'm.')
       plot(out.t{is}(out.an_ind{is}), out.anz_dn{is}, 'g.')
       
       set(gca, 'YDir', 'reverse',...
           'XLim', datenum([out.t{is}(1)-0.5  out.t{is}(end)+0.5]),...
           'YLim',         [out.z{is}(1)-0.1  out.z{is}(end)+0.1]); axis tight;
       colorbar; colormap(ax1, darkb2r(min(out.dT{is}(:)), ceil(100*max(out.dT{is}(:)))/100));

       ax2 = subplot(2,1,2); hold on;
       imagesc(out.t{is}, out.z{is}, out.T{is});
       plot(out.t{is}                , out.F0{is}(2,:), 'w')
       plot(out.t{is}(out.an_ind{is}), out.anz_up{is}, 'k.')
       plot(out.t{is}(out.an_ind{is}), out.anz_dn{is}, 'g.')
       
       set(gca, 'YDir', 'reverse',...
           'XLim', datenum([out.t{is}(1)-0.5  out.t{is}(end)+0.5]),...
           'YLim',         [out.z{is}(1)-0.1  out.z{is}(end)+0.1]); axis tight;
       colorbar;    colormap(ax2, jet)
       
       pause%(0.5);
end

%% convert temperature differences to water masses
p = 0;
for is = 1:11
    
out.w_z{is} = nan(size(out.an_ind{is},2), 2);
out.w_t{is} = nan(size(out.an_ind{is},2), 2);
out.w_m{is} = nan(size(out.an_ind{is},2), 1);
out.ord{is} = nan(size(out.an_ind{is},2), 1);

for it = 1:size(out.an_ind{is},2) %2:size(out.t{is}, 2)

ind_up = out.an_up {is}(it);
ind_dn = out.an_dn {is}(it);
ind_t  = out.an_ind{is}(it);

% water mass:
w = out.dT{is}(ind_up:ind_dn, ind_t ) .* Cice(out.T{is}(ind_up:ind_dn,ind_t)) .* out.R{is}(ind_up:ind_dn) .* 0.1;
w = max( 0, sum(w) / 333500 );

out.w_t{is}(it,:) =      [ out.t{is}(  ind_t-1)   out.t{is}(  ind_t)];
out.w_z{is}(it,:) = sort([out.F0{is}(1,ind_t-1)  out.F0{is}(1,ind_t)], 'ascend');

if     out.F0{is}(1,ind_t-1) <= out.F0{is}(1,ind_t  ) && out.F0{is}(1,ind_t  ) <= out.F0{is}(2,ind_t  ); out.w_m{is}(it) = w; out.ord{is}(it) = 1;   % IC         -> Field      -> Simulation
elseif out.F0{is}(1,ind_t-1) <= out.F0{is}(2,ind_t  ) && out.F0{is}(2,ind_t  ) <= out.F0{is}(1,ind_t  ); out.w_m{is}(it) = 0; out.ord{is}(it) = 2;   % IC         -> Simulation -> Field
elseif out.F0{is}(1,ind_t  ) <= out.F0{is}(1,ind_t-1) && out.F0{is}(1,ind_t-1) <= out.F0{is}(2,ind_t  ); out.w_m{is}(it) = w; out.ord{is}(it) = 3;   % Field      -> IC         -> Simulation
elseif out.F0{is}(1,ind_t  ) <= out.F0{is}(2,ind_t  ) && out.F0{is}(2,ind_t  ) <= out.F0{is}(1,ind_t-1); out.w_m{is}(it) = 0; out.ord{is}(it) = 4;   % Field      -> Simulation -> IC
elseif out.F0{is}(2,ind_t  ) <= out.F0{is}(1,ind_t-1) && out.F0{is}(1,ind_t-1) <= out.F0{is}(1,ind_t  ); out.w_m{is}(it) = 0; out.ord{is}(it) = 5;   % Simulation -> IC         -> Field
elseif out.F0{is}(2,ind_t  ) <= out.F0{is}(1,ind_t  ) && out.F0{is}(1,ind_t  ) <= out.F0{is}(1,ind_t-1); out.w_m{is}(it) = 0; out.ord{is}(it) = 6;   % Simulation -> Field      -> IC
end

end; clear it ind* w

for i = 1:6
out.order(is,i) = numel( find(out.ord{is} == i ) ) / numel(out.ord{is}) * 100;
end; clear i

if p
    close all; figure; hold on; box on;
    c = jet(size(out.w_m{is},1));
    [w_m_s,ind] = sort( out.w_m{is}(1:end), 'ascend' );
    for i=2:size(out.w_t{is},1)
        plot( [ mean(out.w_t{is}(i,:)) mean(out.w_t{is}(i,:)) ], [out.w_z{is}(i,1) out.w_z{is}(i,2) ], 'Color', c(ind(i),:), 'LineWidth', 3 )
    end; clear i c ind
    plot(out.t{is}, out.F0{is}(1,:))
    caxis([1 size(out.w_m{is},1)])
    colorbar; colormap jet; title(num2str(is)); pause;
end

end; clear is


%% average data
for is = 1:11
w_z = out.w_z{is}; w_z(1,1) = w_z(1,1)+0.000001;
w_m = out.w_m{is};

Za  = out.z{is};   
%floor(w_z(1,1)*10)/10 : 0.1 : ceil(w_z(end,2)*10)/10 ]'; % vertical grid for vertically averaged data
for iz = 1:size(Za,1)-1
    % find w_z depth intervals that 
    ind1 = find( [w_z(:,1) < Za(iz)  ]  &  [Za(iz)   < w_z(:,2)] & [w_z(:,2) < Za(iz+1)] ); % sticking up
    ind2 = find( [Za(iz)   < w_z(:,1)]  &                          [w_z(:,2) < Za(iz+1)] ); % entirely within
    ind3 = find( [Za(iz)   < w_z(:,1)]  &  [w_z(:,1) < Za(iz+1)] & [Za(iz+1) < w_z(:,2)] ); % sticking down
    ind4 = find( [w_z(:,1) < Za(iz)  ]  &                          [Za(iz+1) < w_z(:,2)] ); % sticking both up and down   from the Za(iz)...Za(iz+1) depth interval
    
    Wm(iz) = 0;
    if     isfinite(ind2);       Wm(iz)  = Wm(iz) + nansum(w_m(ind2));
    elseif isfinite(ind1)
        for i = 1:length(ind1);  Wm(iz)  = Wm(iz) + w_m(ind1(i)) * ( w_z(ind1(i),2) - Za(iz)         )  /  ( w_z(ind1(i),2) - w_z(ind1(i),1) );
        end; clear i
    elseif isfinite(ind3)
        for i = 1:length(ind3);  Wm(iz)  = Wm(iz) + w_m(ind3(i)) * ( Za(iz+1)       - w_z(ind3(i),1) )  /  ( w_z(ind3(i),2) - w_z(ind3(i),1) );
        end; clear i
    elseif isfinite(ind4)
        for i = 1:length(ind4);  Wm(iz)  = Wm(iz) + w_m(ind4(i)) * ( Za(iz+1)       - Za(iz)         )  /  ( w_z(ind4(i),2) - w_z(ind4(i),1) );
        end; clear i        
    end
        
end; clear iz w_* *a ind*

out.Wm {is} = Wm';
clear Wm;
end; clear is

%%

save('W_frwd.mat', 'out')