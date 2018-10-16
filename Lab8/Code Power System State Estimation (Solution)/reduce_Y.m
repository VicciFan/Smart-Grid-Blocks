function Y_SE = reduce_Y(Y_LF,Y_sc,idx_LF_slack,idx_LF_neighbor)
% Y_SE = reduce_Y(Y_LF,Y_sc,idx_LF_slack,idx_LF_neighbor)
% 
% INPUT
% - Y_LF            nodal admittance matrix used for the load flow
%                   (i.e., indcluding the slack node)
% - Y_sc            short-circuit admittance of the slack node
% - idx_LF_slack    polyphase index of the slack node
% - idx_LF_neighbor polyphase index of the neighbour of the slack node
% 
% OUTPUT
% - Y_SE            nodal admittance matrix used for state estimation
%                   (i.e., excluding the slack node)

Y_SE = Y_LF;
n_phases = length(idx_LF_slack);

% remove Y_sc from the neighbor node
Y_SE(idx_LF_neighbor,idx_LF_neighbor) = ...
    Y_SE(idx_LF_neighbor,idx_LF_neighbor) - Y_sc*eye(n_phases);

% remove rows/columns related to the slack node
Y_SE(idx_LF_slack,:) = [];
Y_SE(:,idx_LF_slack) = [];

end
