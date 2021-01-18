function [T] = water00_preprocess(cal, n, r, f)
% arguments:
% cal = 0 - raw data, no processing
% cal = 1 - calibation as in TC paper
% cal = 2 - reviewed calibation

% n - the standard deviation of the zero-centered random normally-distributed noise added to the temperature readings
% r and f are the logical variable switching on (=1) or off (=0)
%   r - the interpolation of data to a regular 6 h spaced time grid and 
% 	f - low pass filtering in time using the 12 h threshold
% clear; clc; cal = 0; n = 0; r = 1; f = 1;

%% load temperature data from T-strings from:
% unzip(websave('data_zip', 'http://store.pangaea.de/Publications/Marchenko-etal_2019b/LF_Temp.zip')) % the Pangaea server
load('LF_Temp.mat')                                                                                   % or locally

%  save data to local variables
    T.Tr{ 1} = LF{4}.T.T{1};
    T.off{1} = LF{4}.T.system.off{1   };
    T.zr{ 1} = LF{4}.T.z{1};
    T.tr{ 1} = LF{4}.T.t';
for is = 3:11
    T.Tr{is} = LF{3}.T.T{is-2};
    T.off{is}= LF{3}.T.system.off{is-2   };
    T.zr{is} = LF{3}.T.z{is-2};
    T.tr{is} = LF{3}.T.t';    
end; clear is LF
T.Tr{11}(end,25) = mean(T.Tr{11}(end,[24 26]));

%%
if cal == 0             %  go back to raw data without calibation offsets, reset offsets to zeros
    
    for is = [1 3:11]

        T.Tr{is} = T.Tr{is} + repmat(T.off{is}, 1, size(T.Tr{is},2));
        T.off{is}(:) = 0;
        T.off_fraction{is} = sum(T.Tr{is} == T.off{is}, 2)./length(T.Tr{is})*100; % fraction of all values that are equal to the calibration offset
        
    end
    
end
    
if cal == 1             % stay at the first variant of the calibration offsets

    % introduce offset correction to the sensor n.9 (original depth 5.95 m) at
    % T-string 5 installed in april 2014
    [~, ind1] = min( abs( T.tr{5+2} - datenum( [2014 10 18  0 0 0] ) ) );
    [~, ind2] = min( abs( T.tr{5+2} - datenum( [2014 10 28  0 0 0] ) ) );
    T.off{5+2}(9) = mode( T.Tr{5+2}(9, ind1:ind2) );
    T.Tr{ 5+2}(9, :) = T.Tr{5+2}(9, :) - T.off{5+2}(9); clear ind*
    % fraction of all values that are equal to the calibration offset
    for is = [1 3:11]
    T.off_fraction{is} = sum(T.Tr{is} + T.off{is} == T.off{is}, 2)./length(T.Tr{is})*100;
    end
end

if cal == 2             % go further to the updated calibration offsets and:
    
    % remove the calibration offsets from 2014 data
    for is = [3:11]
            T.Tr{is} = T.Tr{is} + repmat(T.off{is}, 1, size(T.Tr{is},2));
    end
    
    % EXCLUDE data from sensors at:
    is = [   1     2    2     2     2    4     6     6     7     8    9    9     9    9    9    9    ] + 2  ;    % T-str number and
    z  = [  12.75  7.9  9.9  11.9  12.9  8.65  7.9  11.9  11.45  7.9  3.1  9.1  10.1 11.1 12.1 13.1  ] - 0.9;    % depth
    for i = 1:length(is)
            [~, ind] = min(abs( T.zr{ is(i) } - z(i) ));
            T.Tr{  is(i) }(ind,:) = [];
            T.zr{  is(i) }(ind  ) = [];
            T.off{ is(i) }(ind  ) = [];
    end; clear is z i ind
    
    % RESET calibration offset to 0 for sensors at:
    is = [  1      1     1     4    4     4     5     5     5    5     ] + 2  ;    % T-str number and
    z  = [  9.75  10.75 11.75  5.65 6.65  7.65  5.85  6.85  7.85 8.85  ] - 0.9;    % depth
    for i = 1:length(is)
            [~, ind] = min(abs( T.zr{is(i) } - z(i)));
            T.off{ is(i) }(ind  ) = 0;
    end; clear is z i ind
    
    % sensors for which the calibration is to be redefined:
    is = [  4    7 ] + 2;    % T-str number and
    z  = [  9.65 8.45] - 0.9;  % depth
    t1 = [2014 11 26 0 0 0; 2014  9 1 0 0 0;  ]; % start of the calibration period
    t2 = [2014 12 23 0 0 0; 2015  1 1 0 0 0;  ]; % end    ...
    for i = 1:length(is)
            [~, ind_t1] = min(abs( T.tr{ is(i) } - datenum(t1(i,:)) ));
            [~, ind_t2] = min(abs( T.tr{ is(i) } - datenum(t2(i,:)) ));
            [~, ind   ] = min(abs( T.zr{ is(i) } - z(i)        ));
            T.off{ is(i) }(ind  ) = mode( T.Tr{ is(i) }(ind, ind_t1:ind_t2 ) );
    end; clear is z *t* i ind
    
    % go forward and apply the revised calibration offsets
    for is = [3:11]
            T.Tr{is} = T.Tr{is} - repmat(T.off{is}, 1, size(T.Tr{is},2));
    end
    
    % fraction of all values that are equal to the calibration offset
    for is = [1 3:11]
        T.off_fraction{is} = sum(T.Tr{is} + T.off{is} == T.off{is}, 2)./length(T.Tr{is})*100;
    end
    
end



%% apply noise
rng('shuffle');
for is = [1 3:11]
    
    T.Tr{is} = T.Tr{is} + n*randn(size(T.Tr{is},1), size(T.Tr{is},2));
    
end

%%
for is = [1 3:11]
    
    T.t{is} = T.tr{is};
    T.T{is} = T.Tr{is};
    if any(isnan(T.T{is}(:)))
            T.T{is} = fillmissing(T.T{is}, 'linear');
    end
        
    if r                                                              % resample the temperature data to an even in time 6 h-spaced grid
        dt = 6/24; % min(unique( diff(T14.t{1} )) );
        tmp_t = T.t{is};
        tmp_T = T.T{is};
        
        T.t{is} = T.t{is}(1) : dt : T.t{is}(end);
        T.T{is} = [interp1(tmp_t, tmp_T', T.t{is})]';
        
%         figure; hold on;
%         plot(tmp_t  , tmp_T(5,:)  , '.-k')
%         plot(T.t{is}, T.T{is}(5,:), 'o-r')
        
        clear dt tmp_*
    end
    
%     c = lines(2); figure; hold on; box on;
%     h1 = plot(T.t{is},                     T.T{is}         , 'Color', c(1,:));
%     h2 = plot(T.t{is}, lowpassRP( T.t{is}, T.T{is}', 0.5 )', 'Color', c(2,:));
%     datetickzoom('x', 'yyyy-mmm-dd'); xlabel('Time, yyyy-mmm-dd'); ylabel('Temperature, ^oC');
%     legend([h1(1) h2(1)], {'no filter', 'lowpass filtered'})
    
    if f                                                              % lowpass the temperature data using 0.5 day threshold
        T.T{is} = lowpassRP( T.t{is}, T.T{is}', 0.5 )';
%         plot(T.t{is}, T.T{is}(5,:), '*-b')
    end

    T.z{is} = [ ceil(T.zr{is}( 1 )*10)/10 : 0.1 :...
               floor(T.zr{is}(end)*10)/10]';
    tmp =  T.T{is};
    T.T{is} = interp1(T.zr{is}, tmp, T.z{is}, 'pchip', 5);            % intepolate data to a regular 0.1 m spaced grid

%     for it = 1:size(T.t{is},2)
%         plot(T.Tr{is}(:,it), T.zr{is}, 'o'); hold on
%         plot(T.T{ is}(:,it), T.z{ is}, 'LineWidth', 2, 'Color', 'r')
%         title(['String ' num2str(is) ', time ' num2str(it)])
%         set(gca, 'YDir', 'reverse')
%         xlim([-25 2]); ylim([0 12])
%         pause(0.1);  cla;
%     end; clear it

end; clear is

% on string n.9 installed in April 2014 cut interpolated data from 8.3 m and below.
% This excludes data from sensors at the depths of:
%  9.2, 10.2, 11.2,    12.2 m (referenced to glacier surface in April 2014)
% or
% 10.1, 11.1, 12.1 and 13.1 m (referenced to glacier surface in April 2015)
% and helps to avoid a drop in CTT depth
% T.Tr{11}(11:14,:) = [];
% T.zr{11}(11:14,:) = [];
[~, ind] = min(abs( T.z{11} - 8.3 ));
T.z{11}(ind:end  ) = [];
T.T{11}(ind:end,:) = [];

%% horizontal averaging

% vertical and temporal grid
for is = 3:11
    up(is-2) = T.z{is}( 1 );
    dn(is-2) = T.z{is}(end);
end
T.z{2} = [min(up):0.1:max(dn)]';
T.t{2} = T.t{3};
clear is up dn

tmp = nan(size( T.z{2},1), 9, size(T.t{3},2 )); % matrix to be filled by data from 9 individual t-strings
for is = 3:11
    [~, b] = min( abs( T.z{2} - T.z{is}( 1 ) ) );   % index of the upper element
    [~, f] = min( abs( T.z{2} - T.z{is}(end) ) );   %              lower
    tmp(b:f,is-2,:) = T.T{is};
end; clear is b f
T.T{2}  = squeeze( nanmean(tmp,2) );
T.T_sd  = squeeze( std(tmp, 0, 2, 'omitnan') );   % include    omit

% T.tmp = tmp;

end