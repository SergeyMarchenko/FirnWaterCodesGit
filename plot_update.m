function [p] = plot_update(p, t, z, TIC, R, K, F0, to, Wm, a, i, to_f, t_i, TIC_i, d, thr)
%plot_2 updated figure during optimization
if i > 1
      delete(p.fi);    delete(p.oi1);    delete(p.oi2);    delete(p.oi3);
      delete(p.fs1);  delete(p.fs3);
      delete(p.wo); delete(p.frs); delete(p.Q);
end

title(p.ax1, [ {['index of optimized water profile value: ' num2str( i ) ]};...
               {['frfr compared at: ' datestr( t(to(i,1)), 'yy-mmm-dd-hh') ' to ' datestr( t(to(i,2)), 'yy-mmm-dd-hh') ' and '...
                                      num2str(F0(to(i,1)))                 ' to ' num2str(F0(to(i,2))) ' m' ]} ]);

% compared freezing front points, depths of layer where water mass is optimized
p.fi = plot( p.ax3,  t(to(i,1):   to(i,2)) , F0(to(i,1):to(i,2)), 'k.-', 'MarkerSize', 8, 'LineWidth', 1 );
p.oi1 = plot( p.ax3, [t(to(i,2)) t(to(i,2))], [z(i); z(i+1)]    , 'k-' , 'MarkerSize', 8 );
p.oi2 = plot( p.ax3, [t(to(i,2))           ], [z(i)        ]    , 'k^' , 'MarkerSize', 8 );
p.oi3 = plot( p.ax3, [t(to(i,2))           ], [      z(i+1)]    , 'kv' , 'MarkerSize', 8 );
xlim(p.ax3, [ t(max(1,to (i,1) - 10))  t(min(length(t), to (i,2) + 10)) ])
ylim(p.ax3, [ z(max(1,i        -  2))  z(min(length(z), i        +  5)) ])

% freezing front simulated using the optimized water masses    
frfr_i = frfr( TIC_i, z, thr); frfr_i(isnan(frfr_i)) = z(1);
p.fs1 = plot( p.ax1, t(1:to(i,2)), frfr_i, 'LineWidth', 2, 'Color', lines(1));
p.fs3 = plot( p.ax3, t(1:to(i,2)), frfr_i, 'LineWidth', 2, 'Color', lines(1));

% optimized water masses
wv(1:2:2*i,1) = Wm(1:i  )/1000/(z(2) - z(1))*100;
wv(2:2:2*i,1) = Wm(1:i  )/1000/(z(2) - z(1))*100;
wv(1:2:2*i,2) = z (1:i  );
wv(2:2:2*i,2) = z (2:i+1);
p.wo = plot(p.ax2, wv(:,1), wv(:,2), 'LineWidth', 2, 'Color', lines(1));

% "central" guesses for the optimized water mass applied in the bisection routine VS corresponding cost function values AND
% corresponding freezing fronts simulated
c = lines( size(d{i,a},1) );
for it = 1:size(d{i,a},1)
    Wm_it = Wm(:,a); Wm_it( to(i,3) ) = d{i,a}(it,4);
    frfr_it = frfr( Tkw(t(1:to(i,2)), z, TIC(:,1:to(i,2)), R, K, Wm_it, thr), z, thr); frfr_it(isnan(frfr_it)) = z(1);
    
    if   it <  size(d{i,a},1); p.Q(it)   = plot( p.ax4, d{i,a}( it,4), d{i,a}(it,5)     , 'o', 'Color', c(it,:));
                               p.frs(it) = plot( p.ax3, t(1:to(i,2)), frfr_it, 'LineWidth', 1, 'Color', c(it,:));
    else;                      p.Q(it)   = plot( p.ax4, d{i,a}( it,4), d{i,a}(it,5)     , 'o', 'Color', 'k', 'MarkerFaceColor', 'k');
                               p.frs(it) = plot( p.ax3, t(1:to(i,2)), frfr_it, 'LineWidth', 1, 'Color', 'k');
    end
    
    
    
end
% print('-r300', p.fh, ['D:\_PHD\presentations\2019_12_16_ICE\opt\opt_' num2str(i) ], '-djpeg');
% 
clear frfr_i wv it Wm_it frfr_i*
end