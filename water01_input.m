%% LOAD INPUT DATA
clear;
clc;
thr = -0.05;  % temperature threshold to define freezing front position

[T] = water00_preprocess;   % load and preprocess temperature data

% cut data before august and above the depth of 1 m

    [~, tind] = min( abs( T.t{ 1} - datenum([ 2015 07 30 0 0 0 ]) ) );
    [~, zind] = min( abs( T.z{ 1} - 1                             ) );
    T.t{ 1} = T.t{ 1}( tind:end );    T.z{1 } = T.z{ 1}( zind:end, : );    T.T{ 1} = T.T{ 1}( zind:end, tind:end );
    
    [~, tind] = min( abs( T.t{ 2} - datenum([ 2014 07 30 0 0 0 ]) ) );
    for is = 2:11
    [~, zind] = min( abs( T.z{is} - 1                             ) );
    T.t{is} = T.t{is}( tind:end );    T.z{is} = T.z{is}( zind:end, : );    T.T{is} = T.T{is}( zind:end, tind:end );
    end; clear is *ind

% resample the temperature data to an even in time 6 h-spaced grid and lowpass it using 0.5 day threshold
dt = 6/24; % min(unique( diff(T14.t{1} )) );
for is = 1:11
    
    tmp_t = T.t{is};
    tmp_T = T.T{is};
    
    T.t{is} = T.t{is}(1) : dt : T.t{is}(end);
    T.T{is} = [interp1(tmp_t, tmp_T', T.t{is})]';
    
    T.T{is} = lowpassRP( T.t{is}, T.T{is}', 0.5 )';
%     T.T{is} = anstmp; figure; hold on; plot( T.t{is}, T.T{is}(10,:) ); plot( T.t{is}, anstmp(10,:) );
    
end; clear is dt tmp_*

% reset temperature readings >0degC to 0degC
for is = 1:size(T.t,2)
    T.T{is}(T.T{is}>0) = 0;
end; clear is
    
% find [thr] and [-2] deg C isotherms
for is = 1:size(T.t,2)
    T.F0{is} = frfr(T.T{is}, T.z{is}, thr);
    T.F2{is} = frfr(T.T{is}, T.z{is}, -2 );
end; clear is

% cut data before the start of freezing front propagation downwards (1:tind)
for is = 1:size(T.t,2)
    tind = max( find( isnan(T.F0{is} )) );
    
    T.t {is}(  1:tind-1) = [];
    T.F0{is}(  1:tind-1) = [];
    T.F0{is}(  1       ) = T.z{is}(1);
    T.F2{is}(  1:tind-1) = [];
    T.T {is}(:,1:tind-1) = [];
end; clear *ind is

% cut data below the k-optimization domain (zinf:end)
% for is = 1:size(T.t,2)
    % [~, zind] = min( abs( 0.1*ceil(max(  T.F2{is} )*10) - T14.z{is} ) ) ; zind = zind+1;
    % T.T{is}(zind:end,:) = [];
    % T.z{is}(zind:end)   = [];
% end; clear *ind is

% cut data after the freezing front has first hit the bottom of the domain
for is = 1:size(T.t,2)
    ind = min( find( T.F0{is} >= T.z{is}(end) ) );
    T.t {is}(:, ind:end) = [];
    T.T {is}(:, ind:end) = [];
    T.F0{is}(:, ind:end) = [];
    T.F2{is}(:, ind:end) = [];
end; clear is ind

% cut data below the deepest position of the freezing front
for is = 1:size(T.t,2)
    ind = find( 0.1*ceil(max(10*T.F0{is})) == 0.1*round(10*T.z{is},1) ); % here smth is wrong. 0.1*round(10*T.z{is},1) is there in place of T.z{is} because otherwise there were very small differences and the find function would not find the right indexes for strings 9 and 11
    T.z {is}(ind+1:end, : ) = [];
    T.T {is}(ind+1:end, : ) = [];
end; clear is ind

for is = 1:size(T.t,2)
    T.F0{is}(2,:) = NaN;
    T.F0{is}(2,1) = T.F0{is}(1,1);
end; clear is


%%
is = 3;       % choose the thermistor strings: 1 - average of all in 2015, 2 - average of all in 2014, 3:11 separate strings for 2014

close all; figure;
ax1 = subplot(1, 100,  1: 60); hold on;
      uimagesc(  T.t{is}, T.z {is}, T.T{is} )
      plot( ax1, T.t{is}, T.F0{is}(1,:), 'w.' )
      plot( ax1, T.t{is}, T.F2{is}     , 'k.' )
      set(gca, 'YDir', 'reverse'); shading flat; axis tight; colormap jet
ax2 = subplot(1, 100, 70:100); hold on;
for it = 1:size(T.t{is},2)
    bar = plot(ax1, T.t{is}(   it), T.F0{is}(1,it), 'ok'); title(ax1, ['ind ' num2str(it) ', time' datestr(T.t{is}(it), 'yy-mmm-dd-hh')]);
    tmp = plot(ax2, T.T{is}(:, it), T.z{is}, '.-', 'Color', lines(1));
    ctt = plot([min(T.T{is}(:)) max(T.T{is}(:))], [T.F0{is}(1,it) T.F0{is}(1,it)], 'k');
    set(gca, 'YDir', 'reverse', 'XLim', [min(T.T{is}(:)) max(T.T{is}(:))]); xlim([min(T.T{is}(:)) max(T.T{is}(:))]); ylim([T.z{is}(1) T.z{is}(end)]);
    pause(0.2); delete(bar); delete(tmp); delete(ctt);
end; clear it

% ylim([ 1.9 7 ])
% xlim([ datenum([2014 8 15 0 0 0]) datenum([2014 12 1 0 0 0]) ])

clearvars -except T thr

%% load measured density data and interpolate it to the depth grid of temperature data
% load('C:\DATA\rho.mat', 'LF14', 'LF15')
load(        'rho.mat', 'LF14', 'LF15')
R{1} = LF15.rho_reg;
R{2} = LF14.rho_reg;

T.S{1} = LF15.strat;
T.S{2} = LF14.strat;

T.R{1} = interp1(R{1}(:,1), R{1}(:,2), T.z{1});
T.R{2} = interp1(R{2}(:,1), R{2}(:,2), T.z{2});

clear R* LF14 LF15

%% load optimization results
% load('C:\DATA\4\ChenGong\codes\invK8rho8_gamma10_maskedBC_longP.mat', 'K_opt', 'rho_opt', 't_data_opt');
load(                         'invK8rho8_gamma10_maskedBC_longP.mat', 'K_opt', 'rho_opt', 't_data_opt');
N_k = 0;
N_R = 0;

% max depths of domains covered by temperature measurements
    zlim15 =      max(T.z{1});
    zlim14 = 0;
for is = 2:11
    zlim14 = max( max(T.z{is}), zlim14 );
end; clear is

% time-limits of domains used for optimization
T.o_t_s{1} = datevec( t_data_opt{4,1}/24/3600 );
T.o_t_f{1} = datevec( t_data_opt{4,2}/24/3600 );
T.o_t_s{2} = datevec( t_data_opt{3,1}/24/3600 );
T.o_t_f{2} = datevec( t_data_opt{3,2}/24/3600 );

% depths, k and rho values in optimization nodes + add noise of amplitude N_k and N_R
T.o_kR_s{1} = [ K_opt{4,1}( : ,1)    K_opt{4,1}( : ,2) + N_k*rand(8,1) - 0.5*N_k    rho_opt{4,1}( : ,2) + N_R*rand(8,1) - 0.5*N_R ];
T.o_kR_f{1} = [ K_opt{4,2}(1:4,1)    K_opt{4,2}(1:4,2) + N_k*rand(4,1) - 0.5*N_k    rho_opt{4,2}(1:4,2) + N_R*rand(4,1) - 0.5*N_R ];
T.o_kR_s{2} = [ K_opt{3,1}( : ,1)    K_opt{3,1}( : ,2) + N_k*rand(8,1) - 0.5*N_k    rho_opt{3,1}( : ,2) + N_R*rand(8,1) - 0.5*N_R ];
T.o_kR_f{2} = [ K_opt{3,2}(1:6,1)    K_opt{3,2}(1:6,2) + N_k*rand(6,1) - 0.5*N_k    rho_opt{3,2}(1:6,2) + N_R*rand(6,1) - 0.5*N_R ];

% k and R values at the bottom of temperature data domains from linear interpolation of the optimization results
klim15s = polyval( polyfit(T.o_kR_s{1}(:,1), T.o_kR_s{1}(:,2), 1), zlim15 );
Rlim15s = polyval( polyfit(T.o_kR_s{1}(:,1), T.o_kR_s{1}(:,3), 1), zlim15 );

klim15f = polyval( polyfit(T.o_kR_f{1}(:,1), T.o_kR_f{1}(:,2), 1), zlim15 );
Rlim15f = polyval( polyfit(T.o_kR_f{1}(:,1), T.o_kR_f{1}(:,3), 1), zlim15 );

klim14s = polyval( polyfit(T.o_kR_s{2}(:,1), T.o_kR_s{2}(:,2), 1), zlim14 );
Rlim14s = polyval( polyfit(T.o_kR_s{2}(:,1), T.o_kR_s{2}(:,3), 1), zlim14 );

klim14f = polyval( polyfit(T.o_kR_f{2}(:,1), T.o_kR_f{2}(:,2), 1), zlim14 );
Rlim14f = polyval( polyfit(T.o_kR_f{2}(:,1), T.o_kR_f{2}(:,3), 1), zlim14 );


T.o_kR_s{1} = [T.o_kR_s{1}; zlim15 klim15s Rlim15s];
T.o_kR_f{1} = [T.o_kR_f{1}; zlim15 klim15f Rlim15f];
T.o_kR_s{2} = [T.o_kR_s{2}; zlim14 klim14s Rlim14s];
T.o_kR_f{2} = [T.o_kR_f{2}; zlim14 klim14f Rlim14f];


% k and R values on a regular grid
T.or_kR_s{1} = [ 1:0.1:zlim15 ]';
T.or_kR_f{1} = [ 1:0.1:zlim15 ]';
T.or_kR_s{2} = [ 1:0.1:zlim14 ]';
T.or_kR_f{2} = [ 1:0.1:zlim14 ]';

T.or_kR_s{1}(:,2:3) = interp1( T.o_kR_s{1}(:,1), T.o_kR_s{1}(:,2:3), T.or_kR_s{1}, 'linear' );
T.or_kR_f{1}(:,2:3) = interp1( T.o_kR_f{1}(:,1), T.o_kR_f{1}(:,2:3), T.or_kR_f{1}, 'linear' );
T.or_kR_s{2}(:,2:3) = interp1( T.o_kR_s{2}(:,1), T.o_kR_s{2}(:,2:3), T.or_kR_s{2}, 'linear' );
T.or_kR_f{2}(:,2:3) = interp1( T.o_kR_f{2}(:,1), T.o_kR_f{2}(:,2:3), T.or_kR_f{2}, 'linear' );

T.N_k{1} = N_k;    T.N_R{1} = N_R;
T.N_k{2} = N_k;    T.N_R{2} = N_R;

clearvars -except o T* thr


%%
save('W_data.mat')   % data necessaty for quantification of the firn water content


