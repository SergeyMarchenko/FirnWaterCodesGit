function p = plot_background(p, t, z, TIC, F0)
%plot_1 background for the figure to be updated during optimization

close all;    p.fh = figure;
p.ax1 = subplot(4, 100, [       1: 65  100+[ 1: 65] ]);  set(p.ax1, 'YDir', 'reverse', 'YGrid', 'on', 'YMinorGrid', 'on', 'TickDir', 'out');
p.ax2 = subplot(4, 100, [      70:100  100+[70:100] ]);  set(p.ax2, 'YDir', 'reverse', 'YGrid', 'on', 'YMinorGrid', 'on', 'YAxisLocation', 'right');
p.ax3 = subplot(4, 100, [ 200+[ 1: 65] 300+[ 1: 65] ]);  set(p.ax3, 'YDir', 'reverse', 'YGrid', 'on', 'YMinorGrid', 'on', 'TickDir', 'out');
p.ax4 = subplot(4, 100, [ 200+[70:100] 300+[70:100] ]);  set(p.ax4,                    'YGrid', 'on', 'YMinorGrid', 'on', 'YAxisLocation', 'right'); % , 'XScale' , 'log'

hold(p.ax1,'on');
hold(p.ax2,'on');
hold(p.ax3,'on');
hold(p.ax4,'on');

set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1], 'CurrentAxes', p.ax1);

im1 = imagesc(p.ax1, t, z, TIC); alpha(im1, 0.5); colormap jet; colorbar(p.ax1, 'westoutside');
      plot   (p.ax1, t, F0, 'w', 'LineWidth', 3);

im2 = imagesc(p.ax3, t, z, TIC); alpha(im2, 0.5); colormap jet;
      plot   (p.ax3, t, F0, 'w', 'LineWidth', 3);

datetickzoom(p.ax1, 'x',  'yy-mmm-dd');
datetickzoom(p.ax3, 'x', 'yy-mmm-dd');
xlabel      (p.ax1, 'Time, yy-mmm-dd');
xlabel      (p.ax2, 'Water content, vol. %')
xlabel      (p.ax3, 'Time, yy-mmm-dd');
xlabel      (p.ax4, 'Guesses of the water mass in current layer, kg m^-^2');
ylabel      (p.ax1, 'Depth, m');
ylabel      (p.ax2, 'Depth, m');
ylabel      (p.ax3, 'Depth, m');
ylabel      (p.ax4, 'Cost function \Sigma (Z_{frfr sim} - Z_{frfr obs}), m');
xlim        (p.ax1, [t(1) t(end)])
ylim        (p.ax1, [z(1) z(end)])
xlim        (p.ax2, [0 3])
ylim        (p.ax2, [z(1) z(end)])

linkaxes([p.ax1,p.ax2],'y')

end