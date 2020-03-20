%% F01_Temp
clearvars
clc; close all;

load('W_data.mat');

p.p_f   = [ 17    3    17  12.5];
p.p_ax1 = [  1.5  7.4  15   4  ];
p.p_ax2 = [  1.5  2.6  15   4  ];

p.p_cb  = [  1.5  1.1  15   0.3];
p.p_lg  = [  3.3  8.2     1   1  ];

p.p_an1 = [ p.p_ax1(1:2) + p.p_ax1(3:4) - [1 1.2]  1  1 ];
p.p_an2 = [ p.p_ax2(1:2) + p.p_ax2(3:4) - [1 1.2]  1  1 ];
p.p_an3 = [ p.p_ax1(1:2) +                [7.3 -0.3]  1  1 ];
p.p_an4 = [ p.p_ax1(1:2) +                [8.5 -0.3]  1  1 ];
p.p_an5 = [ p.p_ax2(1:2) +                [7.3 -0.3]  1  1 ];
p.p_an6 = [ p.p_ax2(1:2) +                [8.5 -0.3]  1  1 ];

close all;
p.f = figure( 'units', 'centimeters', 'Position', p.p_f   );  colormap jet;
p.ax1 = axes( 'units', 'centimeters', 'position', p.p_ax1, 'TickLength',[0.005, 0.01], 'FontSize', 9, 'FontWeight', 'bold', 'YDir', 'reverse', 'TickDir', 'out', 'box', 'on' );
p.ax2 = axes( 'units', 'centimeters', 'position', p.p_ax2, 'TickLength',[0.005, 0.01], 'FontSize', 9, 'FontWeight', 'bold', 'YDir', 'reverse', 'TickDir', 'out', 'box', 'on' );

hold(p.ax1,'on');
hold(p.ax2,'on');

imagesc(p.ax1, T.t{2 }, T.z { 2}+0.9, T.T{2});
imagesc(p.ax2, T.t{1 }, T.z { 1}    , T.T{1});


for is = 3:11
    plot(p.ax1, T.t{is}, T.F0{is}+0.9, 'w', 'LineWidth', 1);
end
    plot(p.ax1, T.t{ 2}, T.F0{ 2}+0.9, 'm', 'LineWidth', 3);

    plot(p.ax2, T.t{ 1}, T.F0{ 1}    , 'w', 'LineWidth', 1);
    
    plot(p.ax1, T.t{ 2}, T.F2{ 2}+0.9, 'k', 'LineWidth', 1);
    plot(p.ax2, T.t{ 1}, T.F2{ 1}    , 'k', 'LineWidth', 1);
    
annotation('textbox', 'String', 'A', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.p_an1);
annotation('textbox', 'String', 'B', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.p_an2);
annotation('textbox', 'String', '2014', 'Color', 'w', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.p_an3);
annotation('textbox', 'String', '2015', 'Color', 'w', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.p_an4);
annotation('textbox', 'String', '2015', 'Color', 'w', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.p_an5);
annotation('textbox', 'String', '2016', 'Color', 'w', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.p_an6);
                  
cb = colorbar(p.ax1, 'location', 'southoutside', 'units', 'centimeters', 'Position', p.p_cb, 'TickDirection', 'out', 'TickLength',[0.005]);
cb.Label.String = 'Temperature, ^oC';
         
caxis(p.ax1, [ -10 0])
caxis(p.ax2, [ -10 0])

     xlabel(p.ax2, 'Time, mmm-dd');
yl = ylabel(p.ax1, 'Depth, m');
set(yl, 'units', 'centimeters', 'Position', [-0.75 -0.7 0]);

xlim        (p.ax1, datenum([2014 8 22 18 0 0; 2015 4 12 12 0 0])); % xlim        (p.ax1, [T.t{2}(1)-1 T.t{2}(end)+1])
xlim        (p.ax2, datenum([2015 8 22 18 0 0; 2016 4 12 12 0 0])); % xlim        (p.ax2, [T.t{1}(1)-1 T.t{1}(end)+1])

datetickzoom(p.ax1, 'x', 'mmm-dd', 'keeplimits'); xtickangle(p.ax1, 0)
datetickzoom(p.ax2, 'x', 'mmm-dd', 'keeplimits');
set(p.ax2, 'xticklabel',{[]} )

ylim        (p.ax1, [T.z{2}(1)-0.2 T.z{2}(end)+0.2]+0.9)
ylim        (p.ax2, [T.z{1}(1)-0.2 T.z{1}(end)+0.2]    )

p.lg0_thick = plot(p.ax1, [1 2], [1 2], 'm', 'LineWidth', 3);
p.lg0_thin  = plot(p.ax1, [1 2], [1 2], 'w', 'LineWidth', 1);
p.lg_fusk   = plot(p.ax1, [1 2], [1 2], 'w', 'LineStyle', 'none');
p.lg2       = plot(p.ax1, [1 2], [1 2], 'k', 'LineWidth', 1);
p.lg = legend( [p.lg0_thick p.lg0_thin, p.lg_fusk, p.lg2], {'spatially averaged data', 'individual T-strings', '', '-2 ^oC isotherm'}, 'box', 'off',...
                'units', 'centimeters', 'Position', p.p_lg, 'TextColor', 'w');
p.lg.Title.String = 'Freezing fronts from:';
drawnow
p.lg.Title.NodeChildren.Position = p.lg.Title.NodeChildren.Position + [-0.15 0 0];

p.f.InvertHardcopy = 'off';  % keep the white font of the legend text entires
set(gcf,'color','w');

% print(gcf, 'C:\DATA\4\text\WCont_TeX\fig\F01_Temp.jpg','-djpeg','-r600');

%% F02_RK
clearvars
clc; %close all;
load('W_data.mat')
bh = load('BHrepr.mat', 'z', 'prob14');
bh.z = -bh.z;

p.p_f    = [ 17    3     17   14 ];
p.p_ax1  = [  1.5  3.5    7   10 ];
p.p_ax2  = [  9    3.5    7   10 ];
p.p_lg1  = [  6.8  4      1    1 ];
p.p_lg2  = [  5    0.7    1    1 ];
p.p_lg3  = [ 12    0.7    1    1 ];

p.p_an1 = [ p.p_ax1(1:2) + p.p_ax1(3:4) - [1 1.2]  1  1 ];
p.p_an2 = [ p.p_ax2(1:2) + p.p_ax2(3:4) - [1 1.2]  1  1 ];

close all;
p.f = figure( 'Units', 'centimeters', 'Position', p.p_f   );  colormap jet;
p.ax1f= axes( 'Units', 'centimeters', 'Position', p.p_ax1, 'YDir', 'reverse', 'TickDir', 'out', 'XLim', [250  910 ], 'YLim', [1 13.4],...
                     'FontSize', 9, 'FontWeight', 'bold', 'YTick', [], 'XTick', [277.5  327.5], 'XTickLabel', [2014 2015], 'TickLength',[0, 0.01], 'Color', 'none' );
p.ax1 = axes( 'Units', 'centimeters', 'Position', p.p_ax1, 'YDir', 'reverse', 'TickDir', 'out', 'XLim', [250  910 ], 'YLim', [1 13.4],...
                     'FontSize', 9, 'FontWeight', 'bold', 'XTick', [400 500 600 700 800 900], 'TickLength',[0.005, 0.01] );

xtickangle(p.ax1f, 90)
p.ax2 = axes( 'Units', 'centimeters', 'Position', p.p_ax2, 'FontSize', 9, 'FontWeight', 'bold',...
    'TickLength',[0.005, 0.01], 'YAxisLocation', 'right', 'YDir', 'reverse', 'TickDir', 'out', 'XLim', [0.55 1.65], 'YLim', [1 13.4] );

hold(p.ax1,'on'); box(p.ax1,'on');
hold(p.ax2,'on'); box(p.ax2,'on');
    
annotation('textbox', 'String', 'A', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.p_an1);
annotation('textbox', 'String', 'B', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.p_an2);

xlabel      (p.ax1, 'Density, kg m^-^3');
xlabel      (p.ax2, 'Thermal conductivity, W m^-^1 K^-^1');
ylabel      (p.ax1, 'Depth, m');

c = [lines( 4 ); [176 215 247;...
                   65 166 249;...
                    3  77 137;...
                  237 177  32]/255];
                
A1 = find(T.or_kR_f{2}(:,1) == T.o_kR_f{2}( end-1, 1));
A2 = find(T.or_kR_f{1}(:,1) == T.o_kR_f{1}( end-1, 1));

set(p.f, 'CurrentAxes', p.ax1);
    fill(350+bh.prob14*200        , bh.z                      + 0.9, c(8,:), 'EdgeColor', c(8,:) );
    plot(T.R      {2}             , T.z      {2}              + 0.9, '-', 'Color', c(1,:)                 )   % measured
    plot(T.o_kR_f {2}(1:end-1  ,3), T.o_kR_f {2}(1:end-1  ,1) + 0.9, 'o', 'Color', c(1,:), 'LineWidth', 2 )   % optimized
    plot(T.o_kR_f {2}(  end    ,3), T.o_kR_f {2}(  end    ,1) + 0.9, 'o', 'Color', 'k'   , 'LineWidth', 2 )   % linearly extrapolated to bottom
    plot(T.or_kR_f{2}(    1:A1 ,3), T.or_kR_f{2}(    1:A1 ,1) + 0.9, '-', 'Color', c(1,:), 'LineWidth', 2 )   % regular grid
    plot(T.or_kR_f{2}( A1+1:end,3), T.or_kR_f{2}( A1+1:end,1) + 0.9, '-', 'Color', 'k'   , 'LineWidth', 2 )   % regular grid extrapolated

    plot(T.R      {1}             , T.z      {1}              + 0.0, '-', 'Color', c(2,:)                 )
    plot(T.o_kR_f {1}(1:end-1  ,3), T.o_kR_f {1}(1:end-1  ,1) + 0.0, 'o', 'Color', c(2,:), 'LineWidth', 2 )
    plot(T.o_kR_f {1}(  end    ,3), T.o_kR_f {1}(  end    ,1) + 0.0, 'o', 'Color', 'k'   , 'LineWidth', 2 )
    plot(T.or_kR_f{1}(    1:A2 ,3), T.or_kR_f{1}(    1:A2 ,1) + 0.0, '-', 'Color', c(2,:), 'LineWidth', 2 )   % regular grid
    plot(T.or_kR_f{1}( A2+1:end,3), T.or_kR_f{1}( A2+1:end,1) + 0.0, '-', 'Color', 'k'   , 'LineWidth', 2 )   % regular grid extrapolated
    
    for iz = 1:max( size(T.S{2},1), size(T.S{1},1) )
        
        if iz <= size(T.S{2},1); tmp2 = T.S{2}(iz,:) + [0.9 0.9 0]; end
        if iz <= size(T.S{1},1); tmp1 = T.S{1}(iz,:)              ; end
        
        tmp2(tmp2(1:2) <= 1.05) = 1.05; tmp2(tmp2(1:2) >= 13.35) = 13.35;
        tmp1(tmp1(1:2) <= 1.05) = 1.05; tmp1(tmp1(1:2) >= 13.35) = 13.35;
        
        fill([255 255 300 300], [tmp2(1) tmp2(2) tmp2(2) tmp2(1) ],...
                   c(tmp2(3)+4,:), 'EdgeColor', c(tmp2(3)+4,:), 'LineWidth', 0.5);
        fill([305 305 350 350], [tmp1(1) tmp1(2) tmp1(2) tmp1(1) ],...
                   c(tmp1(3)+4,:), 'EdgeColor', c(tmp1(3)+4,:), 'LineWidth', 0.5);
    end; clear iz tmp*
    plot([350 350], [1 13.4], 'k', 'LineWidth', 1)
    

set(p.f, 'CurrentAxes', p.ax2);
    plot(T.o_kR_f {2}(1:end-1  ,2), T.o_kR_f {2}(1:end-1  ,1) + 0.9, 'o', 'Color', c(1,:), 'LineWidth', 2 )
    plot(T.o_kR_f {2}(  end    ,2), T.o_kR_f {2}(  end    ,1) + 0.9, 'o', 'Color', 'k'   , 'LineWidth', 2 )
    plot(T.or_kR_f{2}(    1:A1 ,2), T.or_kR_f{2}(    1:A1 ,1) + 0.9, '-', 'Color', c(1,:), 'LineWidth', 2 )   % regular grid
    plot(T.or_kR_f{2}( A1+1:end,2), T.or_kR_f{2}( A1+1:end,1) + 0.9, '-', 'Color', 'k'   , 'LineWidth', 2 )   % regular grid extrapolated
    
    plot(T.o_kR_f {1}(1:end-1  ,2), T.o_kR_f {1}(1:end-1  ,1) + 0.0, 'o', 'Color', c(2,:), 'LineWidth', 2 )
    plot(T.o_kR_f {1}(  end    ,2), T.o_kR_f {1}(  end    ,1) + 0.0, 'o', 'Color', 'k'   , 'LineWidth', 2 )
    plot(T.or_kR_f{1}(    1:A2 ,2), T.or_kR_f{1}(    1:A2 ,1) + 0.0, '-', 'Color', c(2,:), 'LineWidth', 2 )   % regular grid
    plot(T.or_kR_f{1}( A2+1:end,2), T.or_kR_f{1}( A2+1:end,1) + 0.0, '-', 'Color', 'k'   , 'LineWidth', 2 )   % regular grid extrapolated
    
    kC   = plot(keffcalonne(T.o_kR_f {2}(:,3)), T.o_kR_f {2}(:,1) + 0.9    , 'k');
    kR   = plot(keffriche  (T.o_kR_f {2}(:,3)), T.o_kR_f {2}(:,1) + 0.9    , 'r');
    kC19 = plot(keffcalonne2019(T.o_kR_f {2}(:,3)), T.o_kR_f {2}(:,1) + 0.9, 'm');

l1 =  plot([1 1], -[ 7  8], '-' , 'Color', c(1,:) );                   % measured 2014
l2 =  plot([1 1], -[ 8  9], '-' , 'Color', c(2,:) );                   % measured 2015
l3 =  plot([1 1], -[ 9 10], '-o', 'Color', c(1,:), 'LineWidth', 2 );   % optimized 2014
l4 =  plot([1 1], -[10 11], '-o', 'Color', c(2,:), 'LineWidth', 2 );   % optimized 2015
l5 =  plot([1 1], -[11 12], '-o', 'Color', 'k'   , 'LineWidth', 2 );   % linearly extrapolated to bottom
lsnow = fill(-[2 2 1 1], -[2 2 1 1], c(5,:), 'EdgeColor', c(5,:), 'LineWidth', 0.5);
lfirn = fill(-[2 2 1 1], -[2 2 1 1], c(6,:), 'EdgeColor', c(6,:), 'LineWidth', 0.5);
lice  = fill(-[2 2 1 1], -[2 2 1 1], c(7,:), 'EdgeColor', c(7,:), 'LineWidth', 0.5);
lprob = fill(-[2 2 1 1], -[2 2 1 1], c(8,:), 'EdgeColor', c(8,:), 'LineWidth', 0.5);

p.lg1 = legend(p.ax1, [ lsnow lfirn lice lprob], {'snow','firn', 'ice', 'P'},...
    'box', 'off', 'units', 'centimeters', 'Position', p.p_lg1, 'NumColumns', 1);
p.lg2 = legend(p.ax1f, [ l3 l4 l5 l1 l2 ], {'optimized 2014', 'optimized 2015', 'linear extrapolation', 'measured 2014', 'measured 2015'},...
    'box', 'off', 'units', 'centimeters', 'Position', p.p_lg2, 'NumColumns', 2);
p.lg3 = legend(p.ax2, [ kC kC19 kR], {'Calonne et al., 2011', 'Calonne at al., 2019', 'Riche and Schneebeli, 2013'},...
    'box', 'off', 'units', 'centimeters', 'Position', p.p_lg3, 'NumColumns', 1);
p.lg1.Title.String = 'Stratigraphy:';
p.lg2.Title.String = 'Density and thermal conductivity:';
p.lg3.Title.String = 'Parameterizations of thermal conductivity:';


% print(gcf, 'C:\DATA\4\text\WCont_TeX\fig\F02_RK.jpg','-djpeg','-r600');


%% F03 Optimization of water mass
clearvars
clc; close all;
load('W_opt.mat')

p.p_f    = [ 17  3     8.6   5 ];
p.p_ax1  = [  1.4  1    7   3.8 ];

close all;
p.f = figure( 'Units', 'centimeters', 'Position', p.p_f   );
p.ax1 = axes( 'Units', 'centimeters', 'Position', p.p_ax1, 'TickDir', 'out', 'TickLength',[0.005, 0.01], 'box', 'on', 'FontSize', 9, 'FontWeight', 'bold' ); hold on

for is = 1:11
p.p(is) = plot(p.ax1, 2:10, sum(diff(out.Wm{is},1,2)), '.-', 'MarkerSize', 15);
end; clear is
set(p.p(2), 'Color', 'k', 'LineWidth', 2)
xlabel(p.ax1, 'Number of optimization iteration, n' )
ylabel(p.ax1, '\Sigma (m^{*}_{n} - m^*_{n-1}), kg m^{-2}' )

print(gcf, 'C:\DATA\4\text\WCont_TeX\fig\F03_opt.jpg','-djpeg','-r600');


%% F04 Forward calculation of firn water mass
clearvars
clc; close all;
load('W_frwd.mat'); thr = -0.05;
is = 8;

[~, indt] = min( abs( out.t{2} - datenum([2014 10 1 0 0 0]) ) );

z = out.z{2};
t = out.t{2}(  indt:indt+50);
K = out.K{2};
R = out.R{2};
W = zeros(size(z)); W(1) = [];

Tm  = out.T{2}(:,indt:indt+50);    %Tm(indz-10:end,:) = 0;
Ts  = Tm;
Ts(2:end-1,2:end-1) = NaN;
Ts = Tkw( t, z, Ts, R, K, W, 0);
Ts = Ts(:,end); clear t K R W

ctt1   = frfr(Tm(:, 1 ), z, thr);
ctt2m  = frfr(Tm(:,end), z, thr);
ctt2s  = frfr(Ts       , z, thr);

c = [0 0 0; lines(1); 0.4660    0.6740    0.1880]; 
close all;
pf   = [ 5    2    17.5   9.5 ];
pax1 = [ 1    2     4     7   ];
pax2 = [ 5.5  2    11     7   ];
pcb  = [ 5.8  2.3   6     0.5 ];
plg  = [12.5  7.5     1     1   ];

pan1 = [ pax1(1:2) + pax1(3:4) - [0.7 1]  1  1 ];
pan2 = [ pax2(1:2) + pax2(3:4) - [0.7 1]  1  1 ];


f  = figure('Units', 'centimeters', 'Position', pf);
ax1 = axes  ('Units', 'centimeters', 'Position', pax1,...
            'YDir', 'reverse', 'box', 'on', 'TickDir', 'out',...
            'FontSize', 9, 'FontWeight', 'bold',...
            'XLim', [-3.5 0.1], 'YTick', [1:1:10], 'YLim', [out.z{is}(1)-0.1 out.z{is}(end)+0.1]); hold on;
        

ax2 = axes  ('Units', 'centimeters', 'Position', pax2,...
            'YDir', 'reverse', 'box', 'on', 'TickDir', 'out','YAxisLocation','right',...
            'FontSize', 9, 'FontWeight', 'bold', 'YTick', [1:1:10],...
            'XLim', datenum([out.t{is}(1)-0.5 out.t{is}(end)+0.5]),...
            'YLim', [out.z{is}(1)-0.1 out.z{is}(end)+0.1]); hold on;

plot(ax1, Tm(:, 1 ), z    , 'Color', c(1,:), 'LineWidth', 2)
plot(ax1, Tm(:,end), z    , 'Color', c(2,:), 'LineWidth', 2)
plot(ax1, Ts       , z    , 'Color', c(3,:), 'LineWidth', 2, 'LineStyle', '--')
plot(ax1,       1  , ctt1 , 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', c(1,:) , 'Marker', '<', 'MarkerSize', 8, 'LineStyle', 'none')
plot(ax1,     0.1  , ctt1 , 'MarkerFaceColor', c(1,:) , 'MarkerEdgeColor', c(1,:) , 'Marker', '<', 'MarkerSize', 8)
plot(ax1,     0.1  , ctt2m, 'MarkerFaceColor', c(2,:) , 'MarkerEdgeColor', c(2,:) , 'Marker', '<', 'MarkerSize', 8)
plot(ax1,     0.1  , ctt2s, 'MarkerFaceColor', c(3,:) , 'MarkerEdgeColor', c(3,:) , 'Marker', '<', 'MarkerSize', 8)

cl = lines(2);
       imagesc(ax2, out.t{is}                , out.z{is}, out.dT{is}); colormap jet;
CTT  = plot   (ax2, out.t{is}                , out.F0    {is}(1,:)     , 'Color', 'k');
anup = plot   (ax2, out.t{is}(out.an_ind{is}), out.anz_up{is}     , '.', 'Color', cl(1,:), 'MarkerSize', 3);
andn = plot   (ax2, out.t{is}(out.an_ind{is}), out.anz_dn{is}     , '.', 'Color', cl(2,:), 'MarkerSize', 3);

lg1 = legend(ax1, 'T_i^m', 'T_{i+1}^m', 'T_{i+1}^s', 'CTT depths');
set(lg1, 'box', 'off', 'Location', 'SouthWest')
xlabel(ax1, 'Temperature, ^oC')
ylabel(ax1, 'Depth, m')

lg2 = legend(ax2, 'measured CTT', 'top of the positive anomaly', 'bottom of the positive anomaly', 'Units', 'centimeters', 'Position', plg);
lg2.ItemTokenSize = [10; 18];


xlabel(ax2, 'Time, Mmm-dd in 2014')
datetickzoom(ax2, 'x', 'mmm-dd', 'keeplimits')

% cb = colorbar(ax2, 'location', 'west', 'units', 'centimeters', 'Position', pcb, 'TickDirection', 'out', 'TickLength', [0.005],...
%     'AxisLocation', 'in');
cb = colorbar(ax2, 'location', 'south', 'units', 'centimeters', 'Position', pcb, 'TickDirection', 'out', 'TickLength', [0.005],...
    'AxisLocation', 'in');
cb.Label.String = '|T^s| - |T^m|, ^oC';
colormap(ax2, darkb2r(min(out.dT{is}(:)), ceil(100*max(out.dT{is}(:)))/100))

annotation('textbox', 'String', 'A', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', pan1);
annotation('textbox', 'String', 'B', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', pan2);

print(gcf, 'C:\DATA\4\text\WCont_TeX\fig\F04_frwd.jpg','-djpeg','-r600');


%% water content results F05
clearvars
clc; close all;

load('W_opt.mat' ); o = out;
load('W_frwd.mat'); f = out; clear out
bh = load('BHrepr.mat', 'z', 'prob14');
bh.z = -bh.z;

for is = 1:11
    
    Wp{is} = [     o.z{is}(1:end-1)+0.05 o.Wm{is}(:,5) f.Wm{is} ];  % point-wise arranged per T-string
    
    Wp{is}(:,2:3) = Wp{is}(:,2:3);
    
    tmp1 = Wp{is}(:,1);
    tmp2 = Wp{is}(:,2);
    tmp3 = Wp{is}(:,3);
    
    if rem(size(tmp1,1), 2) == 0
        tmp11 = reshape(tmp1, 2, 0.5*size(tmp1,1));
        tmp21 = reshape(tmp2, 2, 0.5*size(tmp2,1));
        tmp31 = reshape(tmp3, 2, 0.5*size(tmp3,1));
        Wa{is+1} = [mean(tmp11',2) mean(tmp21',2) mean(tmp31',2)];
    else
        tmp12 = reshape(tmp1(1:end-1), 2, 0.5*size(tmp1(1:end-1),1));
        tmp22 = reshape(tmp2(1:end-1), 2, 0.5*size(tmp2(1:end-1),1));
        tmp32 = reshape(tmp3(1:end-1), 2, 0.5*size(tmp3(1:end-1),1));
        Wa{is+1} = [mean(tmp12',2) mean(tmp22',2) mean(tmp32',2);...
                    tmp1(end)      tmp2(end)      tmp3(end)      ];
    end
    
end
Wa{1}  = [Wa{2}; Wa{3}; Wa{4}; Wa{5}; Wa{6}; Wa{7}; Wa{8}; Wa{9}; Wa{10}; Wa{11}; Wa{12}];
ind = or(Wa{1}(:,2)>0.5, Wa{1}(:,3)>0.5);
WaB = Wa{1}(ind,:);
clear is tmp* ind

[R ,~] = corrcoef(Wa{1}(:,2), Wa{1}(:,3));      R   = round(R (2)^2,2);
[RB,~] = corrcoef(WaB  (:,2), WaB  (:,3));      RB  = round(RB(2)^2,2);
md  = mean( Wa{1}(:,3) - Wa{1}(:,2) );          md  = round(md     ,2);
mdB = mean( WaB  (:,3) - WaB  (:,2) );          mdB = round(mdB    ,2);

close all;
pf   = [  5     2    17.5   9.5 ];
pax1 = [  1     1.3     5     8   ];
pax2 = [  6.3   1.3     5     8   ];
pax3 = [ 12.4   4.3     5     5   ];

pan1 = [ pax1(1:2) + pax1(3:4) - [0.7 1]  1  1 ];
pan2 = [ pax2(1:2) + pax2(3:4) - [0.7 1]  1  1 ];
pan3 = [ pax3(1:2) + pax3(3:4) - [1.1 1]  1  1 ];
pan4 = [ pax3(1:2) + pax3(3:4) - [4.95   1.2]    2.15  1.1 ];
pan5 = [ pax3(1:2) + pax3(3:4) - [4.95   1.8  ]  2.15  1 ];

plg  = [ 14   1.1     1 1   ];

f  = figure('Units', 'centimeters', 'Position', pf);
hax1 = axes  ('Units', 'centimeters', 'Position', pax1,...
            'YDir', 'reverse', 'YGrid', 'on', 'box', 'on', 'TickDir', 'out',...
            'FontSize', 9, 'FontWeight', 'bold',...
            'XLim', [0 2.1], 'XTick', [0:0.5:2], 'YTick', [1:2:11], 'YLim', [1 12]); hold on;
hax2 = axes  ('Units', 'centimeters', 'Position', pax2,...
            'YDir', 'reverse', 'YGrid', 'on', 'box', 'on', 'TickDir', 'out', 'YAxisLocation','right',...
            'YTickLabel', [], 'FontSize', 9, 'FontWeight', 'bold',...
            'XLim', [-0.1 6], 'YTick', [1:2:11], 'YLim', [1 12]); hold on;
hax3 = axes  ('Units', 'centimeters', 'Position', pax3,...
            'box', 'on', 'TickDir', 'out',...
            'FontSize', 9, 'FontWeight', 'bold',...
            'XLim', [-0.1 3.05], 'YLim', [-0.1 3.05], 'XTick', [0:1:3], 'YTick', [0:1:3]); hold on;
        
c = lines(2);
hwpo15 = plot(hax1, Wp{1}(:,2), Wp{1}(:,1)    , 'o-', 'MarkerSize', 4, 'LineWidth', 1, 'Color', c(2,:)  );
hwpo14 = plot(hax1, Wp{2}(:,2), Wp{2}(:,1)+0.9, 'o-', 'MarkerSize', 4, 'LineWidth', 1, 'Color', c(1,:)  );


% fill(hax2, bh.prob14*6, bh.z+0.9, [237 177  32]/255, 'EdgeColor', [237 177  32]/255 );
c = linspecer(9,'qualitative');
for is = 3:11
    ind = Wp{is}(:,2) >= 0;
    plot(hax2, Wp{is}(ind,2), Wp{is}(ind,1)+0.9, '.-', 'MarkerSize', 5, 'LineWidth', 1, 'Color', c(is-2,:)); %lines(1)
end; clear is ind
hwpf14i = plot(hax2, -[8 9], -[8 9], '.-', 'MarkerSize', 5, 'LineWidth', 1, 'Color', 'k');
hwp = plot(hax2, 100*0.0143*exp( 3.3* [1 - o.R{2}(1:size(o.z{2},1))/900] ).*o.R{2}(1:size(o.z{2},1))./900, o.z{2}+0.9, 'r', 'LineWidth', 2);


        plot(hax3, [0 3], [0 3], 'k', 'LineWidth', 1);
hwof  = plot(hax3, Wa{1}(:,2), Wa{1}(:,3), 'ok', 'MarkerSize', 3);
hwofB = plot(hax3, WaB  (:,2), WaB  (:,3), '.r', 'MarkerSize', 4);
       
xlabel(hax1, '\Theta_w^{v, opt}, vol. %'); ylabel(hax1, 'Depth, m');
xlabel(hax2, '\Theta_w^{v, opt}, vol. %')
xlabel(hax3, '\Theta_w^{v, opt}, vol. %')
ylabel(hax3, '\Theta_w^{v, dir}, vol. %')

hlg = legend(hax1, [hwpo15 hwpo14 hwpf14i hwp hwof],...
'fall 2015, individ. T-string', 'fall 2014, horiz.-averaged data', 'fall 2014, individ. T-strings', ['Schneider and Jansson, 2004' newline 'parameterization'], 'optimized vs direct calculation',...
    'Units', 'centimeters', 'Position', plg);
legend(hax1, 'boxoff');


annotation('textbox', 'String', 'A', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', pan1);
annotation('textbox', 'String', 'B', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', pan2);
annotation('textbox', 'String', 'C', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', pan3);
annotation('textbox', 'String', [{['  R^2    MD']};...
                                 {[num2str(R ) '   ' num2str(md )]}], 'Color', 'k', 'LineStyle', 'none', 'FontSize', 8, 'Units', 'centimeters', 'Position', pan4, 'Margin', 2);
annotation('textbox', 'String', [{[num2str(RB) '   ' num2str(mdB)]}], 'Color', 'r', 'LineStyle', 'none', 'FontSize', 8, 'Units', 'centimeters', 'Position', pan5, 'Margin', 2);

% print(gcf, 'C:\DATA\4\text\WCont_TeX\fig\F05_res.jpg','-djpeg','-r600');

%% CTT patterns F06
clearvars; clc;
load('W_data.mat' );
load('W_opt.mat' ); o = out;
load('W_frwd.mat'); f = out; clear out
close all;

is = 2;

wo = o.Wm{is}(:,5);
wf = f.Wm{is};
w0 = zeros(size(o.z{is},1)-1,1);
E = 0.0143*exp(3.3*(1 - o.R{is}./900)); wp = o.R{is} .* 0.1 .* E ./ (1-E); wp = wp(1:size(o.z{is},1)-1);
kp  = keffcalonne2019(o.R{is});

fr0 = frfr(Tkw(o.t{is}, o.z{is}, o.TIC{is}, o.R{is}, o.K{is}, wo, thr), o.z{is}, thr);
fr1 = frfr(Tkw(o.t{is}, o.z{is}, o.TIC{is}, o.R{is}, o.K{is}, wf, thr), o.z{is}, thr);
fr2 = frfr(Tkw(o.t{is}, o.z{is}, o.TIC{is}, o.R{is}, o.K{is}, w0, thr), o.z{is}, thr);
fr3 = frfr(Tkw(o.t{is}, o.z{is}, o.TIC{is}, o.R{is}, o.K{is}, wp, thr), o.z{is}, thr);
fr4 = frfr(Tkw(o.t{is}, o.z{is}, o.TIC{is}, o.R{is}, kp     , wp, thr), o.z{is}, thr);
clear w E kp

p.p_f   = [ 17    3    17    6.5];
p.p_ax  = [  1.1  1.1  15.5  5  ];
p.p_lg  = [  4    1.5   1    1  ];
p.p_an1 = [ p.p_ax(1:2) + [7.3 -0.3]  1  1 ];
p.p_an2 = [ p.p_ax(1:2) + [8.5 -0.3]  1  1 ];

p.f = figure( 'units', 'centimeters', 'Position', p.p_f   );  colormap jet;
p.ax1 = axes( 'units', 'centimeters', 'position', p.p_ax, 'TickLength',[0.005, 0.01], 'YDir', 'reverse', 'TickDir', 'out', 'FontSize', 9, 'FontWeight', 'bold' );
hold on; box on;

plot(o.t{ is}, 0.9+T.F0{is}(1,:), '-k', 'LineWidth', 2); % measured from hor-aver. data
plot(o.t{ is}, 0.9+fr0     , '-m', 'LineWidth', 1); % measured from hor-aver. data
plot(o.t{ is}, 0.9+fr1     , '-' , 'LineWidth', 2); % w    - directly-calculated
plot(o.t{ is}, 0.9+fr2     , '-' , 'LineWidth', 2); % w    - 0
plot(o.t{ is}, 0.9+fr3     , '-' , 'LineWidth', 2); % w    - parameterization
plot([0 1]  , [0 1]       , 'w'                 ); % fake line for legend
plot(o.t{ is}, 0.9+fr4     , '-' , 'LineWidth', 2); % w, k - parameterization

xlabel('Time in 2014 and 2015, mmm-dd');
ylabel('Depth, m');

ylim        ([o.z{is}(1)-0.2 o.z{is}(end)+0.2]+0.9)
xlim        ([o.t{is}(1)-1   o.t{is}(end)+1])
datetickzoom('x', 'mmm-dd', 'keeplimits')

p.lg  = legend( {'\Theta_w^{v, opt.}',...
                '\Theta_w^{v, dir.}',...
                '\Theta_w^v = 0',...
                '\Theta_w^{v, param.}',...
                '',...
               ['\Theta_w^{v, param.}, ' 'k_{eff}^{param.}']},... % newline 
                'box', 'off', 'units', 'centimeters', 'Position', p.p_lg, 'TextColor', 'k','NumColumns', 3);
p.lg.Title.String = 'CTT from T simulated using:';

p.lg.ItemTokenSize = [10; 18];
drawnow
p.lg.Title.NodeChildren.Position = p.lg.Title.NodeChildren.Position + [-0.2 0 0];

title(['is = ' num2str(is) ', Wo = ' num2str(sum(wo)) ', Wf = ' num2str(sum(wf))])

% p.f.InvertHardcopy = 'off';  % keep the white font of the legend text entires
% set(gcf,'color','w');

% print(gcf, 'C:\DATA\4\text\WCont_TeX\fig\F06_frfr.jpg','-djpeg','-r600');
