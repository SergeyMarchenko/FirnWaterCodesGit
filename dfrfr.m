function out = dfrfr(t, z, TIC, R, K, Wm, thr, to_s, to_f, g, i)
% difference between depths of simulated and measured freezing fronts

Wm( i ) = g;        % apply the water mass guess [g] to where it belongs [g_loc] and assign 0 to whatever nans are left

[Ts, ~] = Tkw( t(1:to_f), z, TIC(:,1:to_f), R, K, Wm, thr );

F0s = frfr(Ts , z, thr); F0s( isnan(F0s) ) = z(1);
F0m = frfr(TIC, z, thr); F0m( isnan(F0m) ) = z(1);
out = F0s(to_s:to_f) - F0m(to_s:to_f);
out = nansum(out);

end