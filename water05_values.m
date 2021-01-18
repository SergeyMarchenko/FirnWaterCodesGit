%% irreducible volumetric water content values from literature
clear; clc;

r = [400 700]';

100*SchnJan04(r, 'v')

Si = 0.02;                                  % set the irreducible water saturation (Jordan et al., 2008: 0.01; Colbeck, Anderson, 1982: 0.07)
i = 1;
n (:,i) = 1 - (r - 0)/910;                  % initial guesses of: porosity
w (:,i) = repmat(0.5, length(r), 1);        %                     the volumetric water content
dw(:,i) = repmat(0.5, length(r), 1);        %                     the differernce in volumetric water content for neighboring steps

while max(dw) > 0.0001
      i = i + 1;
      n(:,i) = 1 - (r - w(:,i-1).*r )/910;  % update the estimate of porosity
      w(:,i) = Si .* n(:,i);                % update the estimate of volumetric water content 
      dw(:,i) = abs( w(:,i) - w(:,i-1) );
end

w = 100*w(:,end);

%% statistics on calibration offsets
clear; clc;
T = water00_preprocess(2, 0, 1, 1);   % load and preprocess temperature data
off = [];
for is = [1 3:11]
    
    off = [off; repmat(is, size(T.zr{is})) T.zr{is}+0.9 T.off{is} T.off_fraction{is}];
           
end; clear is

disp(['Number of sensors: '      num2str(length(    off        ))]);
disp(['Number of zero offsets: ' num2str(   sum(    off(:,3)==0))]);
disp(['Mean absolute offset: '   num2str(  mean(abs(off(:,3)   ))) '^oC']);
disp(['Mean offset: '            num2str(  mean(    off(:,3)    )) '^oC']);
disp([ 'Max offset: '            num2str(   max(    off(:,3)    )) '^oC']);



%% rates of CTT propagation
clearvars
clc; close all;
load('W_data.mat', 'T');

% figure; hold on
for is = 1:11
    F0{is}  = T.F0{is}(1,1:4:end);
    dF0{is} = 100*diff(F0{is});
%     dF0{is}(dF0{is}>30) = NaN;
    s(is,:) = [min(dF0{is}) nanmean(dF0{is}) max(dF0{is}) std(dF0{is},'omitnan')];
%     pause
%     plot(dF0{is}, '.-')
end
s(is+1,:) = [min(s(1:is,1)) mean(s(1:is,2)) max(s(1:is,3)) mean(s(1:is,4))];

disp(['The average rate of daily CTT subsidence for all temperature datasets is $'...
    num2str(s(end,2)) '$~\unit{cm~day^{-1}}'...
    'with a prominent variability across time and different T-strings: '...
    'from $'   num2str(s(end,1))...
    '$ to $ '  num2str(s(end,3))...
    '$~\unit{cm~day^{-1}}. While the standard deviation in daily CTT propagation rates averaged for individual T-strings is not large, $ '...
    num2str(std(s(3:11,2))) '$~\unit{cm~day^{-1}}, '])

mi = min( [T.t{3}( 1 ) T.t{4}( 1 ) T.t{5}( 1 ) T.t{6}( 1 ) T.t{7}( 1 ) T.t{8}( 1 ) T.t{9}( 1 ) T.t{10}( 1 ) T.t{11}( 1 )] );
ma = max( [T.t{3}(end) T.t{4}(end) T.t{5}(end) T.t{6}(end) T.t{7}(end) T.t{8}(end) T.t{9}(end) T.t{10}(end) T.t{11}(end)] );
CTTA(1,:) = mi:0.25:ma;
CTTA(2:10,:) = nan;
for is = 3:11

    [~,ind] = min(abs( CTTA(1,:) - T.t{is}(1) ));
    CTTA(is-1,ind:ind+length(T.F0{is})-1) = T.F0{is}(1,:);

end
CTTA(11,:) = std(CTTA(2:10,:),'omitnan');


figure; hold on
plot(CTTA(1,:), -CTTA(2:end,:)', 'k'); datetickzoom('x','mmm-dd')
disp(['the mean standard deviation of simultaneous CTT depths for different T-strings installed in 2014 ($ '...
    num2str(100*mean(CTTA(end,:))) '$~\unit{cm}).'])


% average CTT propagation rates 
[~, n] = min(abs( T.t{1} - datenum([2015 11 12 0 0 0]) ));
tmp = [mean(diff(T.F0{1}(1,  1:4: n )))...
       mean(diff(T.F0{1}(1,n+1:4:end)))];
disp(['Before that date the daily CTT propagation rate is \textit{ca} $'...
    num2str(tmp(1)*100) '$~\unit{~cm~day^{-1}}, while '...
    'after it is $' num2str(tmp(2)*100) '$~\unit{cm~day^{-1}}.'])

%% optimization: % of cases when the cost function is <0 even for water mass = 0
clear; clc;
load('W_opt.mat')
a = nan(1,10);
count = 0;
for is = 1:11
    d_loc = out.d{is}(:,5);
    for iz = 1:length(d_loc)
        if isempty(d_loc{iz}); break; end
        
        a = [a; [is iz d_loc{iz}(1,:)] ];    % outcomes of first iterations of all optimization steps
        
    end
end; clear count is iz d_loc
a(1,:) = [];

outcome_n = [sum(a(:,end) == 1) sum(a(:,end) == 2) sum(a(:,end) == 3)]; % number of outomes of different type

outcome_percent = outcome_n/length(a)*100;                              % percentage of outomes of different type

%%
d = a(find(a(:,end) == 1), : );                                         % outcomes when the cost function is Qi is < 0 even for water mass = 0
is_list = unique(d(:,1));
close all;
for is = 1:length(is_list)
    
    % draw measured temperature evolution and freezing front propagation
    figure; hold on; box on; plotbrowser on; set(gca, 'YDir', 'reverse', 'TickDir', 'out'); colormap jet; 
    datetickzoom('x', 'yyyy-mm-dd'); axis tight;
    imagesc(out.t{is_list(is)}, out.z{is_list(is)}, out.TIC{is_list(is)});
    plot(out.t{is_list(is)}, out.F0{is_list(is)}, 'w', 'LineWidth', 2, 'DisplayName', 'Measured frfr');
    
    iz_list = d(d(:,1) == is_list(is),2);                               % list of numbers of optimization steps
    
    for iz = 1:length(iz_list)
        t_ind_start = out.to{is_list(is)}(iz_list(iz),1);
        t_ind_end   = out.to{is_list(is)}(iz_list(iz),2);
                
        plot(out.t{is_list(is)}(t_ind_start : t_ind_end),...
            out.F0{is_list(is)}(t_ind_start : t_ind_end), '.-k', 'DisplayName', 'Parts of the measured frfr where simulation is too shallow even if water = 0')
        
        Wm = out.Wm{ is_list(is)}(:,5);
%         Wm(iz_list(iz)) = 0.5;
        Tsim = Tkw(out.t{  is_list(is)}(  1:t_ind_end),...              % run forward model up to the optimization step in question
                   out.z{  is_list(is)}               ,...
                   out.TIC{is_list(is)}(:,1:t_ind_end),...
                   out.R{  is_list(is)}               ,...
                   out.K{  is_list(is)}               ,...
                   Wm, out.thr);
              plot(out.t{  is_list(is)}(  1:t_ind_end),...
                   frfr(Tsim,...
                   out.z{  is_list(is)}, out.thr), 'y')
               
           imagesc(out.t{  is_list(is)}(  1:t_ind_end),...
                   out.z{  is_list(is)}               ,...
                   Tsim)
        
        
    end
    
end

%% depth difference between frfr for different ice melt T

clear; clc;
load('W_data.mat');
close all;
for is = 1:11
    CTT_0{  is} = frfr( T.T{ is }, T.z { is }, 0   );
    CTT_thr{is} = frfr( T.T{ is }, T.z { is }, thr );
    dCTT{is} = CTT_0{  is} - CTT_thr{is};
    ind = find(dCTT{is}<1.5);
    dCTTa{is} = dCTT{is}(ind);
    
    figure; hold on; set(gca, 'YDir', 'reverse'); title(is)
    plot(T.t{ is }, CTT_0{  is})
    plot(T.t{ is }, CTT_thr{is})
    plot(T.t{ is }, dCTT{is})
%     plot(T.t{ is }(ind), dCTTa{is})
    
    dCTT_m( is) = nanmean(dCTT{is});
    dCTTa_m(is) = nanmean(dCTTa{is});
end; clear is
mean(dCTT_m)
mean(dCTTa_m)


%% +T calculation
R1 = 332.62;   % resistance expected from the sensor at 1°C in [kOhm]
R0 = 351.02;   % resistance expected from the sensor at 0°C in [kOhm]

% 1/R1 = 1/R0 + 1/Rw
% 1/Rw = 1/R1 - 1/R0
Rw = (R1*R0)/(R0-R1); % resistance of the water film that goes around the sensor in [kOhm]

L = 0.1;   % length    of the water film in [m] (roughly the length of the outer heat shrink tube piece)
H = 10^-4; % thickness of the water film in [m] (0.1 mm)
W = 10^-3; % width     of the water film in [m] (1   mm)

rho =  Rw * H * W / L * 1000; % resistivity of the water film [Ohm m]
C   = 1 / rho;                % electrical conductivity in [S/m]

C    = C*1000/100;  % convert electrical conductivity from [S/m] to [mS/cm]
t    = 0;
p    = 0;
long = 0;
lat  = 0;
S    = gsw_SA_from_SP(gsw_SP_from_C(C,t,p),p,long,lat);    % salinity in [ppt]



%% fractions of different optimization outcomes

clear; clc;
load('W_opt.mat');
d = out.d;

a = nan(10^6,3);
c = 0;
for is = 1:11
for it = 1:size(d{is},1)
    if isempty(d{is}{it,10}); continue; end
    c = c+1;
    a(c,:) = [is it d{is}{it,10}(end, end)];
    
end
end
a(sum(isnan(a)')==3,:) = [];

[ sum(a(:,3) == 1)/length(a) sum(a(:,3) == 2)/length(a) sum(a(:,3) == 3)/length(a) ]

%% errors in temperature and CTT depth simulated using water masses from "direct" calculation and optimization
clear; clc;
load W_data.mat;
load W_dir.mat; f = out;
load W_opt.mat; o = out; clear out;

thr = o.thr;
for is = 1:11
    
    
%     dt{is} = f.t{is} - o.t{is};
%     dz{is} = f.z{is} - o.z{is};
%     dR{is} = f.R{is} - o.R{is};
%     dK{is} = f.K{is} - o.K{is};
%     dT{is} = f.T{is} - o.TIC{is};
%     imagesc( dT{is}~=0 ); colorbar; title([num2str(is) ', error: ' num2str(sum(sum(abs(dT{is}))))])
%     pause; clc
    TIC = o.TIC{is};
    wf  = f.Wm{is};
    wo  = o.Wm{is}(:,5);
    
    Tf{is} = Tkw(f.t{is}, f.z{is}, TIC, f.R{is}, f.K{is}, wf, thr);
    To{is} = Tkw(o.t{is}, o.z{is}, TIC, o.R{is}, o.K{is}, wo, thr);
    cttm{is} = frfr(TIC   , f.z{is}, thr);
    cttf{is} = frfr(Tf{is}, f.z{is}, thr);
    ctto{is} = frfr(To{is}, f.z{is}, thr);

      dTf(is) = sqrt(   mean(mean((Tf{is}   - TIC     ).^2)));
      dTo(is) = sqrt(   mean(mean((To{is}   - TIC     ).^2)));
    dCTTf(is) = sqrt(nanmean(     (cttf{is} - cttm{is}).^2));
    dCTTo(is) = sqrt(nanmean(     (ctto{is} - cttm{is}).^2));
    dCTTf_d(is) =    nanmean(      cttf{is} - cttm{is}    ) ;
    dCTTo_d(is) =    nanmean(      ctto{is} - cttm{is}    ) ;
    
end

figure; subplot(2,1,1); hold on; plot(dTf  ); plot(dTo  ); ylabel('RMSD T_{sim} - T_{meas}, ^oC'  ); legend('dir', 'opt');
        subplot(2,1,2); hold on; plot(dCTTf); plot(dCTTo); ylabel('RMSD CTT_{sim} - CTT_{meas}, m'); legend('dir', 'opt');

[mean(  dTf )  mean(  dTo )]
[mean(dCTTf )  mean(dCTTo )]
[mean(dCTTf_d )  mean(dCTTo_d )]

[ min( dTf )  max( dTf ) min( dTo )  max( dTo ) ]


%%

clear
clc
load('C:\DATA\4\FirnWaterCodes\SensitivityExperiments\W_opt_005.mat'); out005 = out; clear out;
load('W_opt.mat');


close all;
for is = 1:11

    a(is,:) = [sum(out.Wm{is}(:,5)) sum(out005.Wm{is}(:,5))];
    
    figure; plotbrowser on;
    ax1 = subplot(2,1,1); hold on; title(ax1, ['[0.05] - [0.02] = ' num2str(a(is,2) - a(is,1)) ])
    ax2 = subplot(2,1,2); hold on;
    
    plot(ax1,    out.z{is}(2:end), cumsum(   out.Wm{is}(:,5)));
    plot(ax1, out005.z{is}(2:end), cumsum(out005.Wm{is}(:,5)));

    
    plot(ax2,    out.F0{is},    out.t{is});
    plot(ax2, out005.F0{is}, out005.t{is});
 
end




