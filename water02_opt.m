% for nit = 1:2
tic
clearvars -except nit
clc; close all;

% Assign parameter values
A = 1 ;        % number of iterations of the whole optimization routine
N = 3;          % max number of layers included in the forward model above the optimized layer (the less the more is the reduction of the domain length at the start)
L = 0;          % min possible water mass in a 0.1 m thick layer
H = 10;         % max possible water mass in a 0.1 m thick layer (tested with 10 first)
tol  = 0.001;   % stopping ctireria of the bisection routine. Bisection stops  if it attempts to take 
p.f  = 0;       % to plot (1) or not to plot (0), that is the question
is_list = 1:11; % string number: 1 - single string from 2015, 2 - averaged data from 2014, 3...11 - strings 1...9 from 2014

load('W_data.mat');

out.Wm = cell( 1, 11 );
for is = is_list

t   = T.t{is};
z   = T.z{is};
TIC = T.T{is};

% reset subfreezing temperature in IC for all depths but the upper BC
tmp = TIC(2:end,1); tmp( tmp<thr ) = thr;
      TIC(2:end,1) = tmp; clear tmp;

F0 = frfr(TIC, z, thr)';
F0( isnan(F0) ) = z(1);

R = T.or_kR_f{1*(is == 1) + 2*(is ~= 1)}(:,3);              % optimized density
% R = T.R{1*(is == 1) + 2*(is ~= 1)};                       % measured density
K = T.or_kR_f{1*(is == 1) + 2*(is ~= 1)}(:,2);              % optimized thermal conductivity
% K = keffcalonne2019(R);                                   % parameterized thermal conductivity

% prescribe initial firn water content and the "DO or DO NOT optimize" flag
Wm   = nan  ( numel( z )-1, 1 );
Wm_f = false( numel( z )-1, 1 );            % 0 -> optimize, 1 -> do not optimize the corresponding Wmi value
Wm  ( z(1:end-1)>F0(end) ) = 0;             % below the deepest point of the freezing front: water mass - 0, do not optimize
Wm_f( z(1:end-1)>F0(end) ) = true;
d = cell(size(Wm,1),A);

% indexes of time vector elements when the freezing front first enters 0.1 m depth intervals
to = nan(size(Wm,1),2);
for i = 1:sum(Wm_f == 0)-1
    to(i,2) = find( [F0 - z(1+i)] < 0, 1, 'last');
end
to(  i+1,2) = length( F0 );
to(1    ,1) = 2;
to(2:i+1,1) = to(1:i,2)+1;
to  (:  ,3) = 1:size(to,1);
clear i

ind = find(diff(to(:,1:2),1,2)<0);           % time indexes that "do not go forward in time"
Wm  (ind)   = 0;
Wm_f(ind)   = true;
to  (ind,:) = nan;
to(isnan(to(:,1)),:) = [];
Wm  (:,2:A) = NaN;

%%
% tic
for a = 1:A;    if a > 1;  Wm(:,a) = Wm(:,a-1);  end
    
    if  p.f; close all; p = plot_background(p, t, z, TIC, F0); end

for i = 1:length(to)
    
%     i = i+1;
    if  i<N;    Wm_i =   Wm(:,a);
    else;       Wm_i =   Wm_i(:,to(i-N+1,2));                                   end
   
    if  i<N;   TIC_i =                        TIC(:,1            : to(i,2)) ;
    else;      TIC_i = [TIC_i(:,to(i-N+1,2))  TIC(:,to(i-N+2,1)  : to(i,2))];     end

    if  i<N;     t_i =                          t(  1            : to(i,2)) ;
    else;        t_i =                          t(  to(i-N+2,1)-1: to(i,2)) ;     end

    if  i<N;   to_s  = to(i,1);
    else;      to_s  = to(i,1) - (to(i-N+2,1)-2);                                 end
               to_f  = to_s + diff(to(i,1:2));
               
    f = @(g)dfrfr(t_i, z, TIC_i, R, K, Wm_i, thr, to_s, to_f, g, to(i,3));  % handle to anonymous function quantifying the cost function
    [ Wm(to(i,3),a), d{i,a}] = bisect_water_mass(f, L, H, tol);
    d{i,a}( sum( isfinite(d{i,a}(:,1)) )+1:end,: ) = [];

    [TIC_i, Wm_i] = Tkw(t(1:to(i  ,2)), z, TIC(:,1:to(i  ,2)), R, K, Wm(:,a), thr);
    
%     disp(['nit = ' num2str(nit) ', string ' num2str(is) ' of 11, a = ' num2str(a) ', iteration ' num2str(i) ' of ' num2str(length(to))])
    disp(['string ' num2str(is) ' of 11, a = ' num2str(a) ', iteration ' num2str(i) ' of ' num2str(length(to))])
    
    
    if p.f;  [p] = plot_update(p, t, z, TIC, R, K, F0, to, Wm, Wm_f, a, i, to_f, t_i, TIC_i, d, thr);
       pause
%         print(  gcf, ['C:\DATA\4\FirnWaterCodes\SensitivityExperiments\CTTtooFast\S' num2str(is) '_z' num2str(i) '.jpg'], '-djpeg','-r600');
%         savefig(gcf, ['C:\DATA\4\FirnWaterCodes\SensitivityExperiments\CTTtooFast\S' num2str(is) '_z' num2str(i)       ], 'compact')

    end

end

end
% toc

out.t   {is} = t;
out.z   {is} = z;
out.TIC {is} = TIC;
out.K   {is} = K;
out.R   {is} = R;
out.F0  {is} = F0;
out.to  {is} = to;
out.Wm  {is} = Wm;
out.Wm_f{is} = Wm_f;
out.d   {is} = d;

clear t z TIC K R F0 to d a i *_i to_* f Wm* g

end
% 
out.thr      = thr;
out.tol      = tol;
out.N        = N;

clearvars -except out nit
% save(['W_opt_nit' num2str(nit) '.mat'], 'out')
% save('W_opt.mat', 'out')
% save('C:\DATA\4\FirnWaterCodes\SensitivityExperiments\W_opt_KCalonne.mat', 'out')
% save('C:\DATA\4\FirnWaterCodes\SensitivityExperiments\W_opt_RhoMeasured.mat', 'out')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% overview of the simulation domains and freezing front points to be compared at each stage of optimization
% close all; figure('units','normalized','outerposition',[0 0 1 1])
% hold on; set(gca, 'YDir', 'reverse', 'TickDir', 'out'); colormap jet; colorbar;
% uimagesc(t, z, TIC);
% plot    (t, F0, 'w', 'LineWidth', 5);
% datetickzoom('x', 'yy-mmm-dd'); axis tight;
% for i = 1:length(to)
% 
%     if Wm_f(i); continue; end
%     title([{['run ' num2str(i) ' of '          num2str( sum(Wm_f==0) ) ',']};...
%            {['time domain: 1 to '              datestr( t(to(i,2)), 'yy-mmm-dd-hh') ',']};...
%            {['frfr compared at: '              datestr( t(to(i,1)), 'yy-mmm-dd-hh') ' to ' datestr( t(to(i,2)), 'yy-mmm-dd-hh') ' and '...
%                                                num2str(F0(to(i,1)))                 ' to ' num2str(F0(to(i,2))) ' m' ]};...
%            {['depths of water mass layers: '   num2str( z(   i  ))                  ' to ' num2str( z(i+1)) ' m']}]);
%             
%             p11 = plot( [t(to(i  ,1)-1)  t(to(i  ,2))], repmat(z(i  ),1,2)      , 'k' , 'LineWidth', 3);    % time domain
%     if i>1; if Wm_f(i-1)==0
%             p12 = plot( [t(to(i-1,1)-1)  t(to(i-1,2))], repmat(z(i-1),1,2)      , 'k' , 'LineWidth', 1);    end; end
%     
%             p21 = plot(  t(to(i  ,1)   :   to(i  ,2)) , F0(to(i  ,1):to(i  ,2) ), 'bo', 'LineWidth', 3);        % compared part of the freezing front propagation pattern
%     if i>1; if Wm_f(i-1)==0
%             p22 = plot(  t(to(i-1,1)   :   to(i-1,2)) , F0(to(i-1,1):to(i-1,2) ), 'bo', 'LineWidth', 1);    end; end
% 
%             p31 = plot(  t(                to(i  ,2)) , z(i  :i+1)              , 'g.', 'LineWidth', 3);        % optimized nodes of the water mass profile
%     if i>1; if Wm_f(i-1)==0
%             p32 = plot(  t(                to(i-1,2)) , z(i-1:i  )              , 'g.', 'LineWidth', 1);    end; end
% 
%             p4  = plot(  t(                to(i  ,2)) , [z(Wm_f); z(Wm_f)+0.1]  , 'go', 'LineWidth', 1);    % NOT optimized nodes of the water mass profile
%     
%     drawnow;    pause%(0.5)
%              delete(p11); delete(p21); delete(p31); delete(p4);
%     if  i>1; delete(p12); delete(p22); delete(p32);  end
%     
% end; clear i p*

% end
toc