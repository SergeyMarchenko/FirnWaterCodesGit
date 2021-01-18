%% calculation of the water content using the forward model
%  as the difference in increment of the refreezing capacity on every time step between measurements and simulation
clearvars
clc; close all;
% for nit = 1:30
% Assign parameter values
p  = 0;   % to plot or not to plot, that is the question
% water01_input
load('W_data.mat');
out.dTcutoff = 0; %0.006;
is_list = 1:11;
%% Temperature profile: simulated from measurements at the previous time steps accounting for conduction only  -  measured
out.thr = thr;
for is = is_list
    
if p; figure('units','normalized','outerposition',[0.5 0 0.5 1]); hold on; grid on; set(gca, 'YDir', 'reverse'); c = lines(2); end

out.t {is} = T.t{is};
out.z {is} = T.z{is};
out.T {is} = T.T{is};
out.F0{is} = T.F0{is}(1,:); out.F0{is}(2,1) = out.F0{is}(1,1);

out.R {is} = T.or_kR_f{1*(is == 1) + 2*(is ~= 1)}(:,3);     % optimized density
% out.R {is} = T.R{1*(is == 1) + 2*(is ~= 1)};              % measured density
out.K {is} = T.or_kR_f{1*(is == 1) + 2*(is ~= 1)}(:,2);     % optimized thermal conductivity
% out.K{is} = keffcalonne2019(out.R{is});                   % parameterized thermal conductivity


W  = zeros(size(out.z{is})); W(end) = [];
out.dT{is}  = nan(size(T.T{is}));
out.Ts{is}  = nan(size(T.T{is}));
out.dT{is}(:,1) = 0;
out.Ts{is}(:,1) = T.T{is}(:,1);
for it = 2:size(out.t{is}, 2)

% disp(['nit = ' num2str(nit) ', string = ' num2str(is) ', t = ' num2str(it) ' of ' num2str(size( out.t{is}, 2))] )
disp(['string = ' num2str(is) ', t = ' num2str(it) ' of ' num2str(size( out.t{is}, 2))] )

Tm  = T.T{is}(:,it     );
Ts  = T.T{is}(:,it-1:it);

[Ts, ~] = Tkw( out.t{is}(it-1:it), out.z{is}, Ts, out.R{is}, out.K{is}, W, thr );

out.F0{is}(2,it) = frfr(Ts(:,2), out.z{is}, thr);    % depth of ctt in simulated data
out.dT{is}(:,it) = abs( Ts(:,2) ) - abs(Tm); % difference of absolute temperatures
out.Ts{is}(:,it) = Ts(:,2);

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

%%  find the indexes and depths of upper and lower boundaries of the positive dT anomaly
for is = is_list
    for it = 1:length(out.t{is})
        
        [an_up(it), an_dn(it)] = find_dTan(out.dT{is}(:,it), out.z{is}, out.F0{is}(1,it), out.dTcutoff); % 0.015
    
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
    
    clear ind* an_*

end; clear is it


%%
% close all;
% for is = 1:11
%    
% figure; plotbrowser on;
% 
% z  = out.z{is};
% t  = out.t{is};
% dT = out.dT{is};
% 
% hold on; box on;
% set(gca, 'YDir', 'reverse', 'TickDir', 'out',...
%          'XLim', [t(1)-2    t(end)+2  ],...
%          'YLim', [z(1)-0.1  z(end)+0.1]);
% 
% 
% cl = lines(4);
% imagesc(t, z, dT);
% plot   (t, out.F0{is}(1,:)   , 'Color', 'k');
% plot   (t(out.an_ind{is})  , out.anz_up{is}  , '.', 'Color', cl(1,:), 'MarkerSize', 3);
% plot   (t(out.an_ind{is})  , out.anz_dn{is}  , '.', 'Color', cl(2,:), 'MarkerSize', 3);
% 
% legend('measured CTT', 'z_t top of the positive anomaly', 'z_b bottom of the positive anomaly');
% 
% xlabel('Time, Mmm-dd')
% datetickzoom('x', 'mmm-dd', 'keeplimits')
% 
% cb = colorbar('location', 'southoutside', 'TickDirection', 'out');
% cb.Label.String = '|T^s| - |T^m|, ^oC';
% colormap(darkb2r(min(dT(:)), ceil(100*max(dT(:)))/100))
%     
% end

%% convert temperature differences to water masses
p = 0;
for is = is_list
    
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
for is = is_list
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
clearvars -except out nit
% save(['W_frwd_nit' num2str(nit) '.mat'], 'out')
save('W_dir.mat', 'out')
% save('C:\DATA\4\FirnWaterCodes\SensitivityExperiments\W_dir_KCalonne.mat', 'out')
% save('C:\DATA\4\FirnWaterCodes\SensitivityExperiments\W_dir_RhoMeasured.mat', 'out')
% end
