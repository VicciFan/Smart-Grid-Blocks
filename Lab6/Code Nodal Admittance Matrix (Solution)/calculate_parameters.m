function [A,Y_L,Y_T] = calculate_parameters(pi_circuits)
% [A,Y_T,Y_L] = calculate_parameters(pi_circuits)
% 
% INPUT
% - pi_circuits (array of struct)
%       pi-section equivalent parameters, stored in fields
%       i       start node
%       j       end node
%       Y_ij    branch admittance
%       Y_i_ij  shunt admittance on start-node side
%       Y_j_ij  shunt admittance on end-node side
%
% OUTPUT
% - A (matrix)
%       branch incidence matrix
% - Y_L (matrix)
%       primitive branch admittance matrix
% - Y_T (matrix)
%       primitive shunt admittance matrix

i = [pi_circuits.i]; % start node
j = [pi_circuits.j]; % end node
Y_ij   = [pi_circuits.Y_ij];    % branch admittance
Y_i_ij = [pi_circuits.Y_i_ij];  % shunt admittance on start-node side
Y_j_ij = [pi_circuits.Y_j_ij];  % shunt admittance on end-node side

n_branches = length(pi_circuits);
n_nodes = max(max(i),max(j));

%% A

A = zeros(n_branches,n_nodes);

for k=1:n_branches
    A(k,i(k)) = +1;
    A(k,j(k)) = -1;
end

%% Y_L

Y_L = zeros(n_branches,1);

for k=1:n_branches
    Y_L(k) = Y_ij(k);
end

Y_L = diag(Y_L);

%% Y_T

Y_T = zeros(n_nodes,1);

for k=1:n_branches
    Y_T(i(k)) = Y_T(i(k)) + Y_i_ij(k);
    Y_T(j(k)) = Y_T(j(k)) + Y_j_ij(k);
end

Y_T = diag(Y_T);

end