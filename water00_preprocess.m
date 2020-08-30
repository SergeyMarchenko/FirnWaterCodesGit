function [T] = water00_preprocess
%water00_prepocess
% - load temperature data from Pangea server
% - interpolate data to a regular 0.1 m spaced vertical grid
% - average 2014 T data horizonatally

%% load temperature data from T-strings and save to local variables
% unzip(websave('data_zip', 'http://store.pangaea.de/Publications/Marchenko-etal_2019b/LF_Temp.zip'))
load('LF_Temp.mat')

    T.Tr{1} = LF{4}.T.T{1};
    T.zr{1} = LF{4}.T.z{1};
    T.tr{1} = LF{4}.T.t';
for is = 3:11
    T.Tr{is} = LF{3}.T.T{is-2};
    T.zr{is} = LF{3}.T.z{is-2};
    T.tr{is} = LF{3}.T.t';
end; clear is LF
T.Tr{11}(end,25) = mean(T.Tr{11}(end,[24 26]));


% introduce offset correction to the sensor n.9 (original depth 5.95 m) at
% T-string 5 installed in april 2014
[~, ind1] = min( abs( T.tr{5+2} - datenum( [2014 10 18  0 0 0] ) ) );
[~, ind2] = min( abs( T.tr{5+2} - datenum( [2014 10 28  0 0 0] ) ) );
T.Tr{5+2}(9, :) = T.Tr{5+2}(9, :) - mode( T.Tr{5+2}(9, ind1:ind2) ); clear ind*

% intepolate data to a regular 0.1 m spaced grid
for is = [1 3:11]
    
    tmp = T.Tr{is};
    if any(isnan(tmp(:)))
            tmp = fillmissing(tmp, 'linear');
    end
    
    T.t{is} = T.tr{is};
    T.z{is} = [ ceil(T.zr{is}( 1 )*10)/10 : 0.1 :...
               floor(T.zr{is}(end)*10)/10]';
    
    T.T{is} = interp1(T.zr{is}, tmp, T.z{is}, 'pchip', 5);

%     for it = 1:size(T.tr{is},1)
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

tmp = nan(size( T.z{2},1), 9, size(T.tr{3},2 )); % matrix to be filled by data from 9 individual t-strings
for is = 3:11
    [~, b] = min( abs( T.z{2} - T.z{is}( 1 ) ) );   % index of the upper element
    [~, f] = min( abs( T.z{2} - T.z{is}(end) ) );   %              lower
    tmp(b:f,is-2,:) = T.T{is};
end; clear is b f
T.T{2}  = squeeze( nanmean(tmp,2) );
T.T_sd  = squeeze( std(tmp, 0, 2, 'omitnan') );   % include    omit

T.tmp = tmp;

end