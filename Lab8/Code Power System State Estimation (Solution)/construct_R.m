function R = construct_R(sd_PMU_V,sd_PMU_I,idx_PMU_V,idx_PMU_I)
% R = construct_R(sd_PMU_V,sd_PMU_I,idx_PMU_V,idx_PMU_I)
%
% INPUT
% - sd_PMU_V    standard deviation of voltage measurement noise
% - sd_PMU_I    standard deviation of current measurement noise
% - idx_PMU_V   polyphase indices of voltage measurement locations
% - idx_PMU_I   polyphase indices of current measurement locations
%
% OUTPUT
% - R           measurement noise covariance matrix

%% Block for voltage measurements

R_V = sd_PMU_V^2 * eye(2*length(idx_PMU_V));

%% Block for current measurements

R_I = sd_PMU_I^2 * eye(2*length(idx_PMU_I));

%% Combine blocks

R = blkdiag(R_V,R_I);

end