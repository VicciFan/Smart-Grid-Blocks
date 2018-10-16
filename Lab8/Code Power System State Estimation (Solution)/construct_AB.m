function [A,B] = construct_AB(idx_SE)
% [A,B] = construct_AB(idx_SE)
%
% INPUT
% - idx_SE  polyphase indices of the estimated buses
% 
% OUTPUT
% - A       state-to-state matrix of the process model
% - B       control-to-state matrix of the process model

n_states = 2 * length(idx_SE);

A = eye(n_states);
B = zeros(n_states,n_states);

end