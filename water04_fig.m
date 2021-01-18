%% F01_Temp
clearvars
clc;
close all;

load('W_data.mat');

p.p_f   = [ 17    3    17  12.5];
p.p_ax1 = [  1.5  8.3  15   4  ];
p.p_ax2 = [  1.5  3.1  15   4  ];

p.p_cb  = [  1.5  1.1  15   0.3]; % [  1.5  1.1  15   0.3];
p.p_lg  = [  3.8  9.4   1   1  ]; % [  3.3  9     1   1  ];

p.p_an1 = [ p.p_ax1(1:2) + p.p_ax1(3:4) - [1 1.2]     1  1 ];
p.p_an2 = [ p.p_ax2(1:2) + p.p_ax2(3:4) - [1 1.2]     1  1 ];
p.p_an3 = [ p.p_ax1(1:2) +                [7.3 -1.5]  1  1 ];
p.p_an4 = [ p.p_ax1(1:2) +                [8.5 -1.5]  1  1 ];
p.p_an5 = [ p.p_ax2(1:2) +                [7.3 -1.5]  1  1 ];
p.p_an6 = [ p.p_ax2(1:2) +                [8.5 -1.5]  1  1 ];

p.p_anl(1,:) = [ p.p_ax1(1:2) + [7 1.8]  1  1 ];

p.f = figure( 'units', 'centimeters', 'Position', p.p_f   );  colormap jet;
p.ax1 = axes( 'units', 'centimeters', 'position', p.p_ax1, 'TickLength',[0.005, 0.01], 'FontSize', 9, 'FontWeight', 'bold', 'YDir', 'reverse', 'TickDir', 'out', 'box', 'on' );
p.ax2 = axes( 'units', 'centimeters', 'position', p.p_ax2, 'TickLength',[0.005, 0.01], 'FontSize', 9, 'FontWeight', 'bold', 'YDir', 'reverse', 'TickDir', 'out', 'box', 'on' );

hold(p.ax1,'on');
hold(p.ax2,'on');

imagesc(p.ax1, T.t{2 }, T.z { 2}+0.9, T.T{2});
imagesc(p.ax2, T.t{1 }, T.z { 1}    , T.T{1});

for is = 3:11
%     plot(p.ax1, repmat(datenum([2014 11 1+is*5 0 0 0]), size(T.zr{is})), T.zr{is}+0.9, 'ok', 'DisplayName', num2str(is))
    plot(p.ax1, T.t{is}, T.F0{is}(1,:)+0.9, 'w', 'LineWidth', 1, 'DisplayName', num2str(is));
end
    plot(p.ax1, T.t{ 2}, T.F0{ 2}(1,:)+0.9, 'm', 'LineWidth', 3);
    

% c = lines(3);
% plot( p.ax1, datenum([2014 8 24 0 0 0]), T.zr{3}(4)+0.9, '.', 'Color', c(1,:))
% plot( p.ax1, datenum([2014 8 24 0 0 0]), T.zr{3}(5)+0.9, '.', 'Color', c(2,:))
% plot( p.ax1, datenum([2014 8 24 0 0 0]), T.zr{3}(6)+0.9, '.', 'Color', c(3,:))
% clear c

%     plot(p.ax1, T.t{ 2}, frfr(T.T{2}, T.z{2}, 0) +0.9, 'g', 'LineWidth', 1);

    plot(p.ax2, T.t{ 1}, T.F0{ 1}(1,:)    , 'w', 'LineWidth', 1);
    
    plot(p.ax1, T.t{ 2}, T.F2{ 2}     +0.9, 'k', 'LineWidth', 1);
    plot(p.ax2, T.t{ 1}, T.F2{ 1}         , 'k', 'LineWidth', 1);
    plot(p.ax2, T.t{ 1}, frfr(T.T{1}, T.z{1}, 0), 'g', 'LineWidth', 1);
    
annotation('textbox', 'String', 'A', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.p_an1);
annotation('textbox', 'String', 'B', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.p_an2);
annotation('textbox', 'String', '2014', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize',  9, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.p_an3);
annotation('textbox', 'String', '2015', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize',  9, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.p_an4);
annotation('textbox', 'String', '2015', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize',  9, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.p_an5);
annotation('textbox', 'String', '2016', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize',  9, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.p_an6);

cb = colorbar(p.ax1, 'location', 'southoutside', 'units', 'centimeters', 'Position', p.p_cb, 'TickDirection', 'out', 'TickLength',[0.005]);
cb.Label.String = 'Temperature, ^oC';
         
caxis(p.ax1, [ -10 0])
caxis(p.ax2, [ -10 0])

xl = xlabel(p.ax2, 'Time, Mmm-dd');
yl = ylabel(p.ax1, 'Depth, m');
set(xl, 'units', 'centimeters', 'Position', [7  -1   0]);
set(yl, 'units', 'centimeters', 'Position', [-0.75 -0.7 0]);

xlim        (p.ax1, datenum([2014 8 22 6 0 0; 2015 4 10 0 0 0])); % xlim        (p.ax1, [T.t{2}(1)-1 T.t{2}(end)+1])
xlim        (p.ax2, datenum([2015 8 22 6 0 0; 2016 4 10 0 0 0])); % xlim        (p.ax2, [T.t{1}(1)-1 T.t{1}(end)+1])

datetickzoom(p.ax1, 'x', 'mmm-dd', 'keeplimits'); xtickangle(p.ax1, 0)
datetickzoom(p.ax2, 'x', 'mmm-dd', 'keeplimits');
% set(p.ax2, 'xticklabel',{[]} )

ylim        (p.ax1, [T.z{2}(1)-0.2 T.z{2}(end)+0.2]+0.9)
ylim        (p.ax2, [T.z{1}(1)-0.2 T.z{1}(end)+0.2]    )


% add arrows labeled with numbers to highlight the drops in CTT
set(p.ax1, 'units', 'normalized'); corner = get(p.ax1, 'Position'); set(p.ax1, 'units', 'centimeters')
xl = xlim(p.ax1);
yl = ylim(p.ax1);

x(1,:) = datenum([2014   9  25 0 0 0; 2014  10   5 0 0 0]);   y(1,:) = [4.8 4.8];
x(2,:) = datenum([2014  10   5 0 0 0; 2014  10  15 0 0 0]);   y(2,:) = [7.3 7.3];
x(3,:) = datenum([2014  11   1 0 0 0; 2014  10  23 0 0 0]);   y(3,:) = [8.3 8.3];

x = corner(1) + abs( x - xl(1) ) .* corner(3) ./ diff(xl);
y = corner(2) + abs( y - yl(2) ) .* corner(4) ./ diff(yl);

annotation( 'textarrow', x(1,:), y(1,:), 'String', '1', 'Color', 'y', 'HeadLength', 8, 'HeadWidth', 8 )
annotation( 'textarrow', x(2,:), y(2,:), 'String', '2', 'Color', 'y', 'HeadLength', 8, 'HeadWidth', 8 )
annotation( 'textarrow', x(3,:), y(3,:), 'String', '3', 'Color', 'y', 'HeadLength', 8, 'HeadWidth', 8 )

p.lg0_thick = plot(p.ax1, [1 2], [1 2], 'm', 'LineWidth', 3);
p.lg0_thin  = plot(p.ax1, [1 2], [1 2], 'w', 'LineWidth', 1);
p.lg_fusk   = plot(p.ax1, [1 2], [1 2], 'w', 'LineStyle', 'none');
p.lg2       = plot(p.ax1, [1 2], [1 2], 'k', 'LineWidth', 1);
p.lg_real0  = plot(p.ax1, [1 2], [1 2], 'g', 'LineWidth', 1);
p.lg = legend( [p.lg2, p.lg_real0, p.lg0_thick p.lg0_thin ],...
    {'T_0=-2^oC', 'T_0=0^oC', 'T_0=-0.03^oC, spatially averaged data', 'T_0=-0.03^oC, individual T-strings' }, 'box', 'off',...
                'units', 'centimeters', 'Position', p.p_lg, 'TextColor', 'w');

% p.lg = legend( [p.lg0_thick p.lg0_thin, p.lg_fusk, p.lg2], {'spatially averaged data', 'individual T-strings', '', '-2 ^oC isotherm'}, 'box', 'off',...
%                 'units', 'centimeters', 'Position', p.p_lg, 'TextColor', 'w');
p.lg.Title.String = 'CTT with:';
drawnow
p.lg.Title.NodeChildren.Position = p.lg.Title.NodeChildren.Position + [-0.35 0 0];
p.lg.ItemTokenSize = [10; 20];

p.f.InvertHardcopy = 'off';  % keep the white font of the legend text entires
set(gcf,'color','w');

% print(gcf, 'C:\DATA\4\text\WCont_TeX\fig\F01_Temp.jpg','-djpeg','-r600');
% plotbrowser on

%% F02_RK
clearvars
clc; %close all;
load('W_data.mat')
bh = load('BHrepr.mat', 'z', 'prob14');
bh.z = -bh.z;

p.p_f    = [ 17    3     17   14 ];
p.p_ax1  = [  1.5  3.5    7   10 ];
p.p_ax2  = [  9    3.5    7   10 ];
p.p_lg1  = [  3.5  4      1    1 ];
p.p_lg2  = [  4.7  0.9    1    1 ];
p.p_lg3  = [ 12.5  0.9    1    1 ];

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
%     fill(350+bh.prob14*200        , bh.z                      + 0.9, c(8,:), 'EdgeColor', c(8,:) );
    plot(T.R      {2}             , T.z      {2}              + 0.9, '-', 'Color', c(1,:)                 )   % measured
    plot(T.o_kR_f {2}(1:end-1  ,3), T.o_kR_f {2}(1:end-1  ,1) + 0.9, 'o', 'Color', c(1,:), 'LineWidth', 2 )   % optimized
    plot(T.o_kR_f {2}(  end    ,3), T.o_kR_f {2}(  end    ,1) + 0.9, '^', 'Color', 'k'   , 'LineWidth', 2 )   % linearly extrapolated to bottom
    plot(T.or_kR_f{2}(    1:A1 ,3), T.or_kR_f{2}(    1:A1 ,1) + 0.9, '-', 'Color', c(1,:), 'LineWidth', 2 )   % regular grid
    plot(T.or_kR_f{2}( A1+1:end,3), T.or_kR_f{2}( A1+1:end,1) + 0.9, '-', 'Color', c(1,:), 'LineWidth', 2 )   % regular grid extrapolated
                                                                                 % 'k'
    plot(T.R      {1}             , T.z      {1}              + 0.0, '-', 'Color', c(2,:)                 )
    plot(T.o_kR_f {1}(1:end-1  ,3), T.o_kR_f {1}(1:end-1  ,1) + 0.0, 'o', 'Color', c(2,:), 'LineWidth', 2 )
    plot(T.o_kR_f {1}(  end    ,3), T.o_kR_f {1}(  end    ,1) + 0.0, '^', 'Color', 'k'   , 'LineWidth', 2 )
    plot(T.or_kR_f{1}(    1:A2 ,3), T.or_kR_f{1}(    1:A2 ,1) + 0.0, '-', 'Color', c(2,:), 'LineWidth', 2 )   % regular grid
    plot(T.or_kR_f{1}( A2+1:end,3), T.or_kR_f{1}( A2+1:end,1) + 0.0, '-', 'Color', c(2,:), 'LineWidth', 2 )   % regular grid extrapolated
                                                                                 % 'k'
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
    plot(T.o_kR_f {2}(  end    ,2), T.o_kR_f {2}(  end    ,1) + 0.9, '^', 'Color', 'k'   , 'LineWidth', 2 )
    plot(T.or_kR_f{2}(    1:A1 ,2), T.or_kR_f{2}(    1:A1 ,1) + 0.9, '-', 'Color', c(1,:), 'LineWidth', 2 )   % regular grid
    plot(T.or_kR_f{2}( A1+1:end,2), T.or_kR_f{2}( A1+1:end,1) + 0.9, '-', 'Color', c(1,:), 'LineWidth', 2 )   % regular grid extrapolated
                                                                                 % 'k'
    plot(T.o_kR_f {1}(1:end-1  ,2), T.o_kR_f {1}(1:end-1  ,1) + 0.0, 'o', 'Color', c(2,:), 'LineWidth', 2 )
    plot(T.o_kR_f {1}(  end    ,2), T.o_kR_f {1}(  end    ,1) + 0.0, '^', 'Color', 'k'   , 'LineWidth', 2 )
    plot(T.or_kR_f{1}(    1:A2 ,2), T.or_kR_f{1}(    1:A2 ,1) + 0.0, '-', 'Color', c(2,:), 'LineWidth', 2 )   % regular grid
    plot(T.or_kR_f{1}( A2+1:end,2), T.or_kR_f{1}( A2+1:end,1) + 0.0, '-', 'Color', c(2,:), 'LineWidth', 2 )   % regular grid extrapolated
                                                                                 % 'k'
    kC   = plot(keffcalonne(T.o_kR_f {2}(:,3)), T.o_kR_f {2}(:,1) + 0.9    , 'k');
    kR   = plot(keffriche  (T.o_kR_f {2}(:,3)), T.o_kR_f {2}(:,1) + 0.9    , 'r');
    kC19 = plot(keffcalonne2019(T.o_kR_f {2}(:,3)), T.o_kR_f {2}(:,1) + 0.9, 'm');

c14 = fill(-[2 2 1 1], -[2 2 1 1], c(1,:), 'EdgeColor', c(1,:), 'LineWidth', 0.5); % color for 2014 data
c15 = fill(-[2 2 1 1], -[2 2 1 1], c(2,:), 'EdgeColor', c(2,:), 'LineWidth', 0.5); % color for 2015 data
l1 =  plot([1 1], -[1 1], '-', 'Color', 'k'                 );                     % measured 2014
l2 =  plot([1 1], -[1 1], 'o', 'Color', 'k', 'LineWidth', 2 );                     % optimized 2014
l3 =  plot([1 1], -[1 1], '^', 'Color', 'k', 'LineWidth', 2 );                     % linearly extrapolated to bottom
l4 =  plot([1 1], -[1 1], '-', 'Color', 'k', 'LineWidth', 2 );                     % used 2014


lsnow = fill(-[2 2 1 1], -[2 2 1 1], c(5,:), 'EdgeColor', c(5,:), 'LineWidth', 0.5);
lfirn = fill(-[2 2 1 1], -[2 2 1 1], c(6,:), 'EdgeColor', c(6,:), 'LineWidth', 0.5);
lice  = fill(-[2 2 1 1], -[2 2 1 1], c(7,:), 'EdgeColor', c(7,:), 'LineWidth', 0.5);
% lprob = fill(-[2 2 1 1], -[2 2 1 1], c(8,:), 'EdgeColor', c(8,:), 'LineWidth', 0.5);

p.lg1 = legend(p.ax1, [ lsnow lfirn lice], {'snow','firn', 'ice'},... %  lprob  , 'P'
    'box', 'off', 'units', 'centimeters', 'Position', p.p_lg1, 'NumColumns', 1);
legend({['blue' char(10) 'line'],'red line'})
p.lg2 = legend(p.ax1f, [ c14 c15 l1 l2 l3 l4 ], {'2014', '2015', 'measurements', 'optimization', ['linear extrapolation' char(10) 'of optimized values'], ['linear interpolation']},...
    'box', 'off', 'units', 'centimeters', 'Position', p.p_lg2, 'NumColumns', 3);
p.lg2.ItemTokenSize = [15; 18];
p.lg3 = legend(p.ax2, [ kC kC19 kR], {'Calonne et al., 2011', 'Calonne at al., 2019', 'Riche and Schneebeli, 2013'},...
    'box', 'off', 'units', 'centimeters', 'Position', p.p_lg3, 'NumColumns', 1);
p.lg3.ItemTokenSize = [15; 18];
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
ylabel(p.ax1, '\Sigma (m^{opt}_{n} - m^{opt}_{n-1}), kg m^{-2}' )

% print(gcf, 'C:\DATA\4\text\WCont_TeX\fig\F03_opt.jpg','-djpeg','-r600');


%% F04 Forward calculation of firn water mass
clearvars
clc; close all;
load('W_dir.mat'); thr = out.thr;
is = 5;

[~, indt] = min( abs( out.t{2} - datenum([2014 9 25 0 0 0]) ) );

z = out.z{2};
t = out.t{2}(  indt:indt+25);
K = out.K{2};
R = out.R{2};
W = zeros(size(z)); W(1) = [];

Tm  = out.T{2}(:,indt:indt+25);    %Tm(indz-10:end,:) = 0;
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
plg  = [13.5  6     1     1   ];

pan1 = [ pax1(1:2) + pax1(3:4) - [0.7 1]  1  1 ];
pan2 = [ pax2(1:2) + pax2(3:4) - [0.7 1]  1  1 ];
pan3 = [ pax2(1:2) + [5    -1.5]  1  1 ];
pan4 = [ pax2(1:2) + [6.3  -1.5]  1  1 ];

f  = figure('Units', 'centimeters', 'Position', pf);
ax1 = axes  ('Units', 'centimeters', 'Position', pax1,...
            'YDir', 'reverse', 'box', 'on', 'TickDir', 'out',...
            'FontSize', 9, 'FontWeight', 'bold',...
            'XLim', [-3.5 0.1], 'YTick', [1:1:11], 'YLim', [out.z{is}(1)-0.1 out.z{is}(end)+0.1]+0.9); hold on;

ax2 = axes  ('Units', 'centimeters', 'Position', pax2,...
            'YDir', 'reverse', 'box', 'on', 'TickDir', 'out','YAxisLocation','right',...
            'FontSize', 9, 'FontWeight', 'bold', 'YTick', [1:1:11],...
            'XLim', datenum([out.t{is}(1)-0.5 out.t{is}(end)+0.5]),...
            'YLim', [out.z{is}(1)-0.1 out.z{is}(end)+0.1]+0.9); hold on;

plot(ax1, Tm(:, 1 ), z    +0.9, 'Color', c(1,:), 'LineWidth', 2)
plot(ax1, Tm(:,end), z    +0.9, 'Color', c(2,:), 'LineWidth', 2)
plot(ax1, Ts       , z    +0.9, 'Color', c(3,:), 'LineWidth', 2, 'LineStyle', '--')
plot(ax1,       1  , ctt1 +0.9, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', c(1,:) , 'Marker', '<', 'MarkerSize', 8, 'LineStyle', 'none')
plot(ax1,     0.1  , ctt1 +0.9, 'MarkerFaceColor', c(1,:) , 'MarkerEdgeColor', c(1,:) , 'Marker', '<', 'MarkerSize', 8)
plot(ax1,     0.1  , ctt2m+0.9, 'MarkerFaceColor', c(2,:) , 'MarkerEdgeColor', c(2,:) , 'Marker', '<', 'MarkerSize', 8)
plot(ax1,     0.1  , ctt2s+0.9, 'MarkerFaceColor', c(3,:) , 'MarkerEdgeColor', c(3,:) , 'Marker', '<', 'MarkerSize', 8)

cl = lines(2);
       imagesc(ax2, out.t{is} , out.z{is}      +0.9, out.dT{is}); colormap jet;
CTT  = plot   (ax2, out.t{is} , out.F0{is}(1,:)+0.9, 'Color', 'k');
ind = out.anz_up{is} == out.anz_dn{is};
t_an = out.t{is}(out.an_ind{is});
anup   = plot (ax2, t_an(~ind), out.anz_up{is}(~ind)+0.9, '.', 'Color', cl(1,:), 'MarkerSize', 3);
andn   = plot (ax2, t_an(~ind), out.anz_dn{is}(~ind)+0.9, '.', 'Color', cl(2,:), 'MarkerSize', 3);
anupdn = plot (ax2, t_an( ind), out.anz_up{is}( ind)+0.9, '.', 'Color', c( 3,:), 'MarkerSize', 3);

lg1 = legend(ax1, 'T_k', 'T_{k+1}', 'T_{k+1}^s', 'CTT depths');
set(lg1, 'box', 'off', 'Location', 'SouthWest')
xlabel(ax1, 'Temperature, ^oC');
ylabel(ax1, 'Depth, m')

lg2 = legend(ax2, [CTT anup andn anupdn],...
    'measured CTT', 'z_t  top of the dT anomaly', 'z_b bottom of the dT anomaly', 'z_t = z_b', 'Units', 'centimeters', 'Position', plg);
lg2.ItemTokenSize = [10; 18];


hxl = xlabel(ax2, 'Time, Mmm-dd');
set(hxl, 'units', 'centimeters', 'Position', [5  -1   0]);
datetickzoom(ax2, 'x', 'mmm-dd', 'keeplimits')

cb = colorbar(ax2, 'location', 'south', 'units', 'centimeters', 'Position', pcb, 'TickDirection', 'out', 'TickLength', [0.005],...
    'AxisLocation', 'in');
cb.Label.String = '|T^s| - |T|, ^oC';
colormap(ax2, darkb2r(min(out.dT{is}(:)), ceil(100*max(out.dT{is}(:)))/100))

annotation('textbox', 'String', 'A', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', pan1);
annotation('textbox', 'String', 'B', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', pan2);
annotation('textbox', 'String', '2014', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize',  9, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', pan3);
annotation('textbox', 'String', '2015', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize',  9, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', pan4);

% print(gcf, 'C:\DATA\4\text\WCont_TeX\fig\F04_frwd.jpg','-djpeg','-r600');


%% F05 readings from several individual sensors at T-strings 1, 4, 5 (2014) resulting in jumps of CTT depths
clear; clc;
load('W_data.mat'); Td = T; clear T;

% data from T-strings: time, depth, temperature
is = [1 4 5];
t = Td.tr{1+2};
z{1} = Td.zr{1+2}(4: 6  ) + 0.9;
z{2} = Td.zr{4+2}(5: 9  ) + 0.9;
z{3} = Td.zr{5+2}(7:12  ) + 0.9;
T{1} = Td.Tr{1+2}(4: 6,:);
T{2} = Td.Tr{4+2}(5: 9,:);
T{3} = Td.Tr{5+2}(7:12,:);

mi = min([z{1}; z{2}; z{3}]);
ma = max([z{1}; z{2}; z{3}]);
cz = jet( round( (ma-mi) / 0.1 + 1 ) );
ic{1} = cz(round( (z{1}-mi)./0.1 )+1,:); % rgb color coding the depth
ic{2} = cz(round( (z{2}-mi)./0.1 )+1,:);
ic{3} = cz(round( (z{3}-mi)./0.1 )+1,:);
    
p.f    = [  7     1.5  17   16.5 ];
p.ax1  = [  1.5  12.3  15    4   ];
p.ax2  = [  1.5   7.5  15    4   ];
p.ax3  = [  1.5   2.5  15    4   ];
p.cb   = [  1.5   1    15    0.3 ];
p.an1 = [ p.ax1(1:2) + p.ax1(3:4) - [0.7 1]     1  1 ];
p.an2 = [ p.ax2(1:2) + p.ax2(3:4) - [0.7 1]     1  1 ];
p.an3 = [ p.ax3(1:2) + p.ax3(3:4) - [0.7 1]     1  1 ];
                  

close all;
h.f     = figure( 'Units', 'centimeters', 'Position', p.f   );
h.ax{1} = axes(   'Units', 'centimeters', 'Position', p.ax1, 'TickLength',[0.005, 0.01], 'FontSize', 9, 'FontWeight', 'bold', 'TickDir', 'out', 'box', 'on', 'Ygrid', 'on' ); hold on;
h.ax{2} = axes(   'Units', 'centimeters', 'Position', p.ax2, 'TickLength',[0.005, 0.01], 'FontSize', 9, 'FontWeight', 'bold', 'TickDir', 'out', 'box', 'on', 'Ygrid', 'on' ); hold on;
h.ax{3} = axes(   'Units', 'centimeters', 'Position', p.ax3, 'TickLength',[0.005, 0.01], 'FontSize', 9, 'FontWeight', 'bold', 'TickDir', 'out', 'box', 'on', 'Ygrid', 'on' ); hold on;

% plot( h.ax1, [t{1}(1) t{1}(end)], [ 0   0 ], 'k')
% plot( h.ax2, [t{1}(1) t{1}(end)], [ 0   0 ], 'k')
% plot( h.ax3, [t{1}(1) t{1}(end)], [ 0   0 ], 'k')

h.thr = plot( h.ax{1}, [t(1) t(end)], [thr thr], 'k--', 'DisplayName', 'ice melt temperature T_0');
        plot( h.ax{2}, [t(1) t(end)], [thr thr], 'k--');
        plot( h.ax{3}, [t(1) t(end)], [thr thr], 'k--');

for i = 1:3
for k = 1:length(z{i})
    h.T{i,k} = plot( h.ax{i}, t, T{i}(k,:), '.-', 'Color', ic{i}(k,:),...
        'DisplayName', ['String ' num2str(is(i)) ', depth ' num2str(z{i}(k)) ' m']);
end
datetickzoom(h.ax{i}, 'x', 'mmm-dd')
end

xlim( h.ax{1}, datenum([ 2014  7 23  0 0 0; 2014 11 10  0 0 0]) )
xlim( h.ax{2}, datenum([ 2014  8 23  0 0 0; 2015  1 20  0 0 0]) )
xlim( h.ax{3}, datenum([ 2014  7 23  0 0 0; 2015  1 20  0 0 0]) )
ylim( h.ax{1}, [-0.4 0.1] )
ylim( h.ax{2}, [-0.4 0.1] )
ylim( h.ax{3}, [-0.4 0.1] )

xlabel(h.ax{3}, 'Time, mmm-dd');
ylabel(h.ax{1}, 'Temperature, ^oC');
ylabel(h.ax{2}, 'Temperature, ^oC');
ylabel(h.ax{3}, 'Temperature, ^oC');

h.cb = colorbar(h.ax{3}, 'Location', 'SouthOutside', 'Units', 'centimeters', 'Position', p.cb, 'TickDirection', 'out', 'TickLength',[0.005] );
caxis(h.ax{3}, [mi ma]); colormap(h.ax{3}, 'jet')
xlabel(h.cb, 'Depth of sensor, m')

annotation('textbox', 'String', 'A', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.an1);
annotation('textbox', 'String', 'B', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.an2);
annotation('textbox', 'String', 'C', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.an3);
                      
% print(gcf, ['C:\DATA\4\text\WCont_TeX\fig\F05_sensors.jpg'], '-djpeg', '-r600');

%% F06 water content results
clearvars
clc; close all;

load('W_opt.mat');  o = out;
load('W_dir.mat');  f = out; 
clear out

bh = load('BHrepr.mat', 'z', 'prob14');
bh.z = -bh.z;

for is = 1:11
    
    Wp{is} = [     o.z{is}(1:end-1)+0.05 o.Wm{is}(:,4) f.Wm{is} ];  % point-wise arranged per T-string
        
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
ind = find(sum([Wa{1}(:,2:3)>0.5]')' > 0);
WaB = Wa{1}(ind,:);
clear is tmp* ind

[R ,~] = corrcoef(Wa{1}(:,2), Wa{1}(:,3));      R   = round(R (2)^2,2);
[RB,~] = corrcoef(WaB  (:,2), WaB  (:,3));      RB  = round(RB(2)^2,2);
md  = mean( Wa{1}(:,3) - Wa{1}(:,2) );          md  = round(md     ,2);
mdB = mean( WaB  (:,3) - WaB  (:,2) );          mdB = round(mdB    ,2);

k1 = polyfit( Wa{1}(:,2), Wa{1}(:,3), 1)
k2 = polyfit( WaB(  :,2), WaB(  :,3), 1)

close all;
pf   = [  5     5    17.5  11 ];
pax1 = [  1     1.5     5     8   ];
pax2 = [  6.3   1.5     5     8   ];
pax3 = [ 12.4   4.5     5     5   ];

pan1 = [ pax1(1:2) + pax1(3:4) - [0.7 1]  1  1 ];
pan2 = [ pax2(1:2) + pax2(3:4) - [0.7 1]  1  1 ];
pan3 = [ pax3(1:2) + pax3(3:4) - [1.1 1]  1  1 ];
pan4 = [ pax3(1:2) + pax3(3:4) - [4.95   1.2]    2.15  1.1 ];
pan5 = [ pax3(1:2) + pax3(3:4) - [4.95   1.8  ]  2.15  1 ];

plg  = [ 14    1.3     1    1 ];

hf  = figure('Units', 'centimeters', 'Position', pf);
hax1 = axes('Units', 'centimeters', 'Position', pax1,...
            'YDir', 'reverse', 'YGrid', 'on', 'box', 'on', 'TickDir', 'out',...
            'FontSize', 9, 'FontWeight', 'bold',...
            'XLim', [0 2.2], 'XTick', [0:0.4:2.4], 'YTick', [1:2:11], 'YLim', [1 12]); hold on;
hax2 = axes('Units', 'centimeters', 'Position', pax2,...
            'YDir', 'reverse', 'YGrid', 'on', 'box', 'on', 'TickDir', 'out', 'YAxisLocation','right',...
            'YTickLabel', [], 'FontSize', 9, 'FontWeight', 'bold',...
            'XLim', [-0.1 3.8], 'YTick', [1:2:11], 'YLim', [1 12]); hold on;
        tmp = get(hax2,'XLim');
% hax2f= axes('Units', 'centimeters', 'Position', pax2,...%[pax2(1)+(2-tmp(1))*pax2(3)/diff(tmp) pax2(2)+pax2(4) (tmp(2)-2)*pax2(3)/diff(tmp) 0.01],...
%             'YDir', 'reverse', 'YTickLabel', [], 'YTick', [1:2:11], 'YLim', [1 12],...
%             'XDir', 'reverse', 'XLim', [0 diff(tmp)], 'XTick', [0 0.9 1.8], 'XTickLabel', [0 0.5 1],...
%             'Color', 'none', 'TickDir', 'out', 'XAxisLocation', 'top', 'FontSize', 9, 'FontWeight', 'bold');
hax3 = axes('Units', 'centimeters', 'Position', pax3,...
            'box', 'on', 'TickDir', 'out',...
            'FontSize', 9, 'FontWeight', 'bold',...
            'XLim', [-0.1 3.05], 'YLim', [-0.1 3.05], 'XTick', [0:1:3], 'YTick', [0:1:3]); hold on;

c = lines(7);
hwpo14 = plot(hax1, Wp{2}(:,2), Wp{2}(:,1)+0.9, '.-', 'MarkerSize', 4, 'LineWidth', 1, 'Color', 'k'  );
hwpf14 = plot(hax1, Wp{2}(:,3), Wp{2}(:,1)+0.9, '.-', 'MarkerSize', 4, 'LineWidth', 1, 'Color', c(6,:)  );

tmp = get(hax2, 'XLim');
% hP = fill(hax2, tmp(2)-1.8*bh.prob14, bh.z+0.9, [237 177  32]/255, 'EdgeColor', [237 177  32]/255 ); clear tmp; %set(hP, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5); 
for is = 3:11
    ind = Wp{is}(:,2) >= 0;
    hwpo14s = plot(hax2, Wp{is}(ind,2), Wp{is}(ind,1)+0.9, '.-', 'MarkerSize', 5, 'LineWidth', 1, 'Color', c(1,:)); %c = linspecer(9,'qualitative'); c(is-2,:)
end; clear is ind
hwpo15  = plot(hax2, Wp{1}(:,2), Wp{1}(:,1), '.-', 'MarkerSize', 4, 'LineWidth', 1, 'Color', c(2,:) );
hwpo141 = plot(hax2, Wp{2}(:,2), Wp{2}(:,1)+0.9, '.-', 'MarkerSize', 4, 'LineWidth', 1, 'Color', 'k'  );
hwp = plot(hax2, 100*SchnJan04(o.R{2}(1:size(o.z{2},1)), 'v'), o.z{2}+0.9, 'Color', c(5,:), 'LineWidth', 2);

        plot(hax3, [0 3], [0 3], 'k', 'LineWidth', 0.5);
hwof  = plot(hax3, Wa{1}(:,2), Wa{1}(:,3), 'ok', 'MarkerSize', 3);
hwofB = plot(hax3, WaB  (:,2), WaB  (:,3), '.r', 'MarkerSize', 4);
% plot(hax3, Wa{1}(:,2), polyval(k1, Wa{1}(:,2)), 'b')
% plot(hax3, WaB  (:,2), polyval(k2, WaB  (:,2)), 'r')
       
xlabel(hax1, '\Theta^{opt, dir}, vol. %'); ylabel(hax1, 'Depth, m');
xlabel(hax2, '\Theta^{opt}, vol. %');
xlabel(hax3, '\Theta^{opt}, vol. %')
ylabel(hax3, '\Theta^{dir}, vol. %')
% hxl = xlabel(hax2f, 'P');
% set(hxl, 'units', 'centimeters', 'Position', [3.8  8.7   0]);

hlg = legend(hax2, [hwpo14 hwpf14 hwpo14s hwpo15 hwp hwof hwofB],... % hP
            'lat.-aver., fall 2014, optimization',...
            'lat.-aver., fall 2014, "direct"',...
            'single T-str., fall 2014',...
            'single T-str., fall 2015',...
            '\Theta_{ir}, Schn., Jan., 2004',...                     % 'P',...
            'all results',...
            '\Theta^{opt} or \Theta^{dir} > 0.5%',...
            'Units', 'centimeters', 'Position', plg, 'NumColumns', 1);
legend(hax2, 'boxoff');
hlg.ItemTokenSize = [20; 10];

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


% print(gcf, 'C:\DATA\4\text\WCont_TeX\fig\F06_Wres.jpg','-djpeg','-r600');

for is = 1:11
m(is,:) = [sum(o.Wm{is}(:,4)) sum(f.Wm{is})];
end
open m
mean(m(3:11,:))
std(m(3:11,:))

%% F07 sensitivity of water content results to the used thermal conductivity (k) and density (r) profiles
clearvars
clc; close all;

load('W_opt.mat'                                                             ); z  = out.z{ 2}(1:end-1)+0.05;
                                                                                o  = out.Wm{2}(:,4);
                                                                                k  = out.K{ 2}(1:length(z));
                                                                                r  = out.R{ 2}(1:length(z));

load('W_dir.mat'                                                            ); f  = out.Wm{2};
load('C:\DATA\4\FirnWaterCodes\SensitivityExperiments\W_opt_KCalonne.mat'    ); ok = out.Wm{2}(:,4);
                                                                                k2 = out.K{ 2}(1:length(z));
load('C:\DATA\4\FirnWaterCodes\SensitivityExperiments\W_opt_RhoMeasured.mat' ); or = out.Wm{2}(:,4);
                                                                                r2 = out.R{ 2}(1:length(z));
load('C:\DATA\4\FirnWaterCodes\SensitivityExperiments\W_dir_KCalonne.mat'   ); fk = out.Wm{2};
load('C:\DATA\4\FirnWaterCodes\SensitivityExperiments\W_dir_RhoMeasured.mat'); fr = out.Wm{2};
clear out

dok = (ok - o)./o*100;
dor = (or - o)./o*100;
dfk = (fk - f)./f*100;
dfr = (fr - f)./f*100;
dk  = (k2 - k)./k*100;
dr  = (r2 - r)./r*100;

close all;
pf   = [  5     2    12   10.5 ];
pax1 = [  1     2.3   5     8   ];
pax2 = [  6.3   2.3   5     8   ];

pan1 = [ pax1(1:2) + pax1(3:4) - [0.7 1]  1  1 ];
pan2 = [ pax2(1:2) + pax2(3:4) - [0.7 1]  1  1 ];

plg1  = [ 1.4  0.1    8.5    1 ];


f  = figure('Units', 'centimeters', 'Position', pf);
hax1 = axes  ('Units', 'centimeters', 'Position', pax1,...
            'YDir', 'reverse', 'XGrid', 'on', 'YGrid', 'on', 'box', 'on', 'TickDir', 'out',...
            'FontSize', 9, 'FontWeight', 'bold',...
            'XLim', [-103 20], 'XTick', [-100:20:100], 'YTick', [2:2:12], 'YLim', [1.8 12.8]); hold on;
hax2 = axes  ('Units', 'centimeters', 'Position', pax2,...
            'YDir', 'reverse', 'XGrid', 'on', 'YGrid', 'on', 'box', 'on', 'TickDir', 'out', 'YAxisLocation','right',...
            'FontSize', 9, 'FontWeight', 'bold',...
            'XLim', [-43 63], 'XTick', [-100:20:100], 'YTick', [2:2:12], 'YLim', [1.8 12.8]); hold on;
% 103 180
c = lines(2);
hdk  = plot(hax1, dk , z+0.9, '.-', 'Color', 'k'   , 'DisplayName', 'dk' );
hdok = plot(hax1, dok, z+0.9, '.-', 'Color', c(1,:), 'DisplayName', 'dok');
hdfk = plot(hax1, dfk, z+0.9, '.-', 'Color', c(2,:), 'DisplayName', 'dfk');

hdr  = plot(hax2, dr , z+0.9, '.-', 'Color', 'k'   , 'DisplayName', 'dr' );
hdor = plot(hax2, dor, z+0.9, '.-', 'Color', c(1,:), 'DisplayName', 'dor');
hdfr = plot(hax2, dfr, z+0.9, '.-', 'Color', c(2,:), 'DisplayName', 'dfr');

linkaxes([hax1 hax2], 'y')

hxl = xlabel(hax1, 'Relative change, %'); ylabel(hax1, 'Depth, m');
set(hxl, 'units', 'centimeters', 'Position', [5.5  -0.8   0]);

hlg1 = legend(hax1, [hdk hdok hdfk],...
['A: thermal conductivity' newline 'B: density'], '\Theta^{opt}', '\Theta^{dir}',...
    'Units', 'centimeters', 'Position', plg1, 'Box', 'off','NumColumns', 3);


annotation('textbox', 'String', 'A', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', pan1);
annotation('textbox', 'String', 'B', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', pan2);

% print(gcf, 'C:\DATA\4\text\WCont_TeX\fig\F07_sensitivity.jpg','-djpeg','-r600');

[mean(    dk (isfinite(dk )))...
 mean(    dok(isfinite(dok)))...
 mean(    dfk(isfinite(dfk)))...
  std(    dok(isfinite(dok)))...
  std(    dfk(isfinite(dfk)))... 
 mean(abs(dr (isfinite(dr ))))...
 mean(abs(dor(isfinite(dor))))...
 mean(abs(dfr(isfinite(dfr))))...
  std(    dor(isfinite(dor)))...
  std(    dfr(isfinite(dfr)))]

%% F08 CTT patterns
clearvars; clc;
load('W_data.mat');
load('W_opt.mat' ); o = out;
load('W_dir.mat' ); f = out; clear out
close all;

is = 2;

wo = o.Wm{is}(:,4);
wf = f.Wm{is};
w0 = zeros(size(o.z{is},1)-1,1);
wp = 100 .* SchnJan04(o.R{is}, 'v'); wp = wp(1:size(o.z{is},1)-1);

kp  = keffcalonne2019(o.R{is});

fr0 = frfr(Tkw(o.t{is}, o.z{is}, o.TIC{is}, o.R{is}, o.K{is}, wo, thr), o.z{is}, thr);
fr1 = frfr(Tkw(o.t{is}, o.z{is}, o.TIC{is}, o.R{is}, o.K{is}, wf, thr), o.z{is}, thr);
fr2 = frfr(Tkw(o.t{is}, o.z{is}, o.TIC{is}, o.R{is}, o.K{is}, w0, thr), o.z{is}, thr);
fr3 = frfr(Tkw(o.t{is}, o.z{is}, o.TIC{is}, o.R{is}, o.K{is}, wp, thr), o.z{is}, thr);
fr4 = frfr(Tkw(o.t{is}, o.z{is}, o.TIC{is}, o.R{is}, kp     , wp, thr), o.z{is}, thr);
clear E kp

p.p_f   = [ 17    3    17    7];
p.p_ax  = [  1.1  1.5  15.5  5  ];
p.p_lg  = [  4    1.8   1    1  ];
p.p_an1 = [ p.p_ax(1:2) + [7.3  -0.3]  1  1 ];
p.p_an2 = [ p.p_ax(1:2) + [8.5  -0.3]  1  1 ];
p.p_an3 = [ p.p_ax(1:2) + [7.8  -1.5]  1  1 ];
p.p_an4 = [ p.p_ax(1:2) + [9    -1.5]  1  1 ];

p.f = figure( 'units', 'centimeters', 'Position', p.p_f   );  colormap jet;
p.ax1 = axes( 'units', 'centimeters', 'position', p.p_ax, 'TickLength',[0.005, 0.01], 'YDir', 'reverse', 'TickDir', 'out', 'FontSize', 9, 'FontWeight', 'bold' );
hold on; box on;

c = lines(6);
p.h1 = plot(o.t{ is}, 0.9+T.F0{is}(1,:), 'Color', c(2,:), 'LineWidth', 3);  % measured from hor-aver. data
p.h2 = plot(o.t{ is}, 0.9+fr0          , 'Color', 'k', 'LineWidth', 2);  % w    - optimized              %c(1,:)
p.h3 = plot(o.t{ is}, 0.9+fr1          , 'Color', c(6,:), 'LineWidth', 2);  % w    - forward calculation    %c(2,:)
p.h4 = plot(o.t{ is}, 0.9+fr2          , 'Color', c(3,:), 'LineWidth', 2);  % w    - 0                      %c(3,:)
p.h5 = plot(o.t{ is}, 0.9+fr3          , 'Color', c(5,:), 'LineWidth', 2);  % w    - parameterization       %c(4,:)
p.h6 = plot([0 1]  , [0 1]             , 'Color', 'w'   , 'LineWidth', 2);  % fake line for the legend
p.h7 = plot(o.t{ is}, 0.9+fr4          , 'Color', 'm', 'LineWidth', 2);  % w, k - parameterization       %c(5,:)

hxl = xlabel('Time, Mmm-dd');
set(hxl, 'units', 'centimeters', 'Position', [8  -1   0]);
ylabel('Depth, m');

ylim        ([o.z{is}(1)-0.2 o.z{is}(end)+0.2]+0.9)
xlim        ([o.t{is}(1)-1   o.t{is}(end)+1])
datetickzoom('x', 'mmm-dd', 'keeplimits')

p.lg  = legend( [p.h2        p.h3        p.h4     p.h5          p.h6 p.h7                          ],...
                {'\Theta^{opt.}', '\Theta^{dir.}', '\Theta = 0', '\Theta^{param.}', '',  '\Theta^{param.}, k_{eff}^{param.}'},... % newline 
                'box', 'off', 'units', 'centimeters', 'Position', p.p_lg, 'TextColor', 'k', 'NumColumns', 3);
p.lg.Title.String = 'CTT from T simulated using:';

p.lg.ItemTokenSize = [10; 18];
drawnow
p.lg.Title.NodeChildren.Position = p.lg.Title.NodeChildren.Position + [-0.2 0 0];

annotation('textbox', 'String', '2014', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize',  9, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.p_an3);
annotation('textbox', 'String', '2015', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize',  9, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.p_an4);

% title(['is = ' num2str(is) ', Wo = ' num2str(sum(wo)) ', Wf = ' num2str(sum(wf))])

% print(gcf, 'C:\DATA\4\text\WCont_TeX\fig\F08_frfr.jpg','-djpeg','-r600');



%%
clear; clc;
T0 = water00_preprocess(0, 0, 0, 0);   % load and preprocess temperature data
T2 = water00_preprocess(2, 0, 0, 0);

% for is = 2:10
% T0.zr{is} = T0.zr{is+1};    T0.Tr{is} = T0.Tr{is+1};    T0.off{is} = T0.off{is+1};
% T2.zr{is} = T2.zr{is+1};    T2.Tr{is} = T2.Tr{is+1};    T2.off{is} = T2.off{is+1};
% end
T0.zr{12} = T0.zr{1};    T0.Tr{12} = T0.Tr{1};    T0.off{12} = T0.off{1};
T2.zr{12} = T2.zr{1};    T2.Tr{12} = T2.Tr{1};    T2.off{12} = T2.off{1};
clear is

%%
p.p_f  = [  3  3    7   13.5 ];
p.p_ax = [  1  3    5   10   ];
p.p_lg = [  3  0.5  1    1   ];

close all;
h.f = figure( 'units', 'centimeters', 'Position', p.p_f );
h.ax = axes( 'units', 'centimeters', 'position', p.p_ax, 'TickLength',[0.005, 0.01], 'FontSize', 9, 'FontWeight', 'bold', 'TickDir', 'out', 'box', 'on',...
    'YDir', 'reverse', 'YTick', [1:12], 'XLim', [0 11], 'XTick', [1:10], 'XTickLabel', [1:9 1]);
hold(h.ax,'on');    c = lines(2);

for is = [3:12]

    plot( h.ax, repmat(is-2, size(T0.zr{is})), T0.zr{is}+0.9*[is<12], 'o-', 'MarkerSize', 3, 'Color', c( 1*[is<12] + 2*[is==12],:) );
    plot( h.ax, repmat(is-2, size(T2.zr{is})), T2.zr{is}+0.9*[is<12], 'o-', 'MarkerSize', 3, 'Color', c( 1*[is<12] + 2*[is==12],:), 'MarkerFaceColor', c( 1*[is<12] + 2*[is==12],:) );

end
h.l1 = plot(-[1 2],[1 2], 'Color', c(1,:), 'LineWidth', 3);
h.l2 = plot(-[1 2],[1 2], 'Color', c(2,:), 'LineWidth', 3);
h.l3 = plot(-[1 2],[1 2], 'o', 'MarkerSize', 5, 'Color', 'k' );
h.l4 = plot(-[1 2],[1 2], 'o', 'MarkerSize', 5, 'Color', 'k', 'MarkerFaceColor', 'k' );

xlabel(h.ax, 'Number of T-string');
ylabel(h.ax, 'Depth referenced to the glacier surface in April 2015, m');

p.lg = legend( h.ax, [h.l1 h.l2 h.l4 h.l3],...
               {'installed in April 2014',...
                'installed in April 2015',...
                'data used in the study',...
                'data NOT used in the study'},...
                'box', 'off', 'Units', 'centimeters', 'Position', p.p_lg);
p.lg.ItemTokenSize = [10; 20];

print(gcf, 'C:\DATA\4\text\WCont_TeX\fig\S01_Tstr_z.jpg','-djpeg','-r600');


%% statistics on calibration offsets and positive measured temperature values (load data)
clear; clc;
load('W_data.mat');    t = T.t;    clear T;
T0 = water00_preprocess(0, 0, 0, 0);   % load and preprocess temperature data
T2 = water00_preprocess(2, 0, 0, 0);   % load and preprocess temperature data

% calibration offsets
     z = T2.zr{1};
     off15 =         [ repmat([2015 1 ], size(z))  z  T2.off{1 } ];
     off14 = [];
for is = 3:11
    
     z = T2.zr{is};
     off14 = [off14; [ repmat([2014 is], size(z)) 0.9+z T2.off{is} ]];
           
end; clear z is

u14 = unique(off14(:,end));
u15 = unique(off15(:,end));
for i = 1:length(u14)
   
    u14n(i) = sum(off14(:,end) == u14(i));
    
end; clear i
for i = 1:length(u15)
   
    u15n(i) = sum(off15(:,end) == u15(i));
    
end; clear i


% all used temperature values
for is = [1 3:11] % cut data as in calculation domains
    
     [~, ind1] = min(abs( T0.t{is} - t{is}( 1 ) ));
     [~, ind2] = min(abs( T0.t{is} - t{is}(end) ));
     
     T0_d{is} = T0.Tr{is}(:, ind1:ind2);
     T2_d{is} = T2.Tr{is}(:, ind1:ind2);
    
end; clear is ind*

T0_15 = T0_d{1}(:);
T2_15 = T2_d{1}(:);

T0_14 = [ T0_d{3}(:); T0_d{4}(:); T0_d{5}(:); T0_d{6}(:); T0_d{7}(:); T0_d{8}(:); T0_d{9}(:); T0_d{10}(:); T0_d{11}(:) ];
T2_14 = [ T2_d{3}(:); T2_d{4}(:); T2_d{5}(:); T2_d{6}(:); T2_d{7}(:); T2_d{8}(:); T2_d{9}(:); T2_d{10}(:); T2_d{11}(:) ];
clear *_d

e014 = round(min(T0_14),2) : 0.1 : round(max(T0_14),2);    N014 = histcounts(T0_14, e014);
e214 = round(min(T2_14),2) : 0.1 : round(max(T2_14),2);    N214 = histcounts(T2_14, e214);
e015 = round(min(T0_15),2) : 0.1 : round(max(T0_15),2);    N015 = histcounts(T0_15, e015);
e215 = round(min(T2_15),2) : 0.1 : round(max(T2_15),2);    N215 = histcounts(T2_15, e215);

%% statistics on calibration offsets and positive measured temperature values (actual plotting)
p.p_f   = [ 17    3    17  15 ];
p.p_ax1 = [  2  9    14   5 ];
p.p_ax2 = [  2  2    14   5 ];

p.p_lg1 = [  3    12.5   1   1 ];
p.p_lg2 = [  3.8  5   1   1 ];

p.p_an1 = [ p.p_ax1(1:2) + p.p_ax1(3:4) - [1 1.2]     1  1 ];
p.p_an2 = [ p.p_ax2(1:2) + p.p_ax2(3:4) - [1 1.2]     1  1 ];

close all;
h.f = figure( 'units', 'centimeters', 'Position', p.p_f );
h.ax1 = axes( 'units', 'centimeters', 'position', p.p_ax1, 'TickLength',[0.005, 0.01], 'FontSize', 9, 'FontWeight', 'bold', 'TickDir', 'out', 'box', 'on');
h.ax2 = axes( 'units', 'centimeters', 'position', p.p_ax2, 'TickLength',[0.005, 0.01], 'FontSize', 9, 'FontWeight', 'bold', 'TickDir', 'out', 'box', 'on', 'YScale', 'log',...
    'YTick', [10^1 10^2 10^3 10^4 10^5], 'XLim', [-11.3 1], 'XTick', [-11:1]);
hold(h.ax1,'on');
hold(h.ax2,'on');

annotation('textbox', 'String', 'A', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.p_an1);
annotation('textbox', 'String', 'B', 'Color', 'k', 'LineStyle', 'none',...
                      'FontSize', 10, 'FontWeight', 'bold',...
                      'Units', 'centimeters', 'Position', p.p_an2);

c = [lines(2); [0 127/255 0]];
h.off14 = plot( h.ax1, u14, u14n, '.', 'MarkerSize', 10                       );
h.off15 = plot( h.ax1, u15, u15n, '.', 'MarkerSize', 10                       );
h.T014  = plot( h.ax2, [e014(1:end-1) + e014(2:end)]/2, N014, 'k'             );
h.T214  = plot( h.ax2, [e214(1:end-1) + e214(2:end)]/2, N214, 'Color', c(1,:) );
h.T015  = plot( h.ax2, [e015(1:end-1) + e015(2:end)]/2, N015, 'Color', c(3,:) );
h.T215  = plot( h.ax2, [e215(1:end-1) + e215(2:end)]/2, N215, 'Color', c(2,:) );

xlabel(h.ax1, 'Calibration offset, ^oC');
ylabel(h.ax1, 'Number of sensors');
xlabel(h.ax2, 'Measured temperature, ^oC')
ylabel(h.ax2, ['Number of readings' newline '(log scale)'])

p.lg1 = legend( h.ax1, [h.off14  h.off15], {'2014', '2015'},...
                'box', 'on', 'Units', 'centimeters', 'Position', p.p_lg1);
p.lg2 = legend( h.ax2, [h.T014 h.T214 h.T015 h.T215], {'2014 raw', '2014 processed', '2015 raw', '2015 processed'},...
                'box', 'on', 'Units', 'centimeters', 'Position', p.p_lg2);

print(gcf, 'C:\DATA\4\text\WCont_TeX\fig\S02_histograms.jpg','-djpeg','-r600');

