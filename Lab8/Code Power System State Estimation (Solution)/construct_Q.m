function Q = construct_Q(Q_var,Q_cov,idx_SE)
% Q = construct_Q(Q_var,Q_cov,idx_SE)
%
% INPUT
% - Q_var   process noise variance cov(w_i,w_i)
% - Q_cov   process noise covariance cov(w_i,w_j), where i~=j
% - idx_SE  polyphase indices of estimated buses
%
% OUTPUT
% - Q       process noise covariance matrix

n_states = 2*length(idx_SE);

K_var = eye(n_states); % ones on diagonal
K_cov = ones(n_states) - eye(n_states); % ones off diagonal

Q = Q_var*K_var + Q_cov*K_cov;

end