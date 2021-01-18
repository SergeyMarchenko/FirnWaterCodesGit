function [an_up, an_dn] = find_dTan(dT, z, F0, cutoff)
% dT = out.dT{is}(:,it); z = out.z{is}; F0 = out.F0{is}(1,it);
% IN:
% dT   - |Tsim| - |Tmeas| profile
% z    - depth vector
% F0o  - observed  CTT depth
% F0s  - simulated CTT depth

% OUT:
% an_up  - upper index of the dT anomaly
% an_dn  - lower index ...

[~, ind] = min( abs( z - F0 ) );                    % index of the layer closest to the CTT

m_up = max( ind-5, 1         );                     % indexes of layers within 0.5 m above and below simulated CTT
m_dn = min( ind+5, length(z) );

[~, ind_max] = max(dT(m_up:m_dn));                  % index of the layer with max dT anomaly within 0.5 from the simulated CTT
ind_max = m_up + ind_max - 1;

an_up = ind_max;
an_dn = ind_max;

for iz = ind_max-1:-1:1                                % find the upper index of the positive dT anomaly
    if    dT(iz) > cutoff && dT(iz) < dT(an_up)
          an_up = iz;
    else;          break
    end
end

for iz = ind_max+1:length(z)                        % find the lowermost index of the positive dT anomaly
    if    dT(iz) > cutoff && dT(iz) < dT(an_dn)
          an_dn = iz;
    else;          break
    end
end

end