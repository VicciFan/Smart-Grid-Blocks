function [V,J,n_iter] = do_LF(Y,S_star,V_0,idx_slack,idx_pq,idx_pv,tol,n_max)
% [V,J,n_iter] = do_LF(Y,S_star,V_0,idx_slack,idx_pq,idx_pv,tol,n_max)
%
% INPUT
% - Y           nodal admittance matrix
% - S_star      given complex power (active and/or reactive power)
% - V_0         initial voltages (phasors)
% - idx_slack   index of the slack bus
% - idx_pq      indices of the PQ buses
% - idx_pv      indices of the PV buses
% - tol         tolerance for convergence criterion
% - n_max       maximum number of iterations
%
% OUTPUT
% - V           solution voltages (phasors)
% - J           Jacobian at the solution
% - n_iter      number of iterations

Y_abs = abs(Y);
Y_arg = angle(Y);
n_variables = length(V_0);

% Initialization
V_abs = abs(V_0);
V_arg = angle(V_0);
J = [];

for k=1:n_max
    n_iter = k;
    
    % Compute nodal voltages/currents/power
    V = complex(V_abs.*cos(V_arg),V_abs.*sin(V_arg));
    I = Y*V;
    S = V.*conj(I);
    
    %% Mismatch calculation
    
    % Compute the mismatches for the entire network.
    dS = S_star-S;
    dP = real(dS);
    dQ = imag(dS);
    
    % Keep only the relevant mismatches.
    dP(idx_slack) = [];
    dQ(sort([idx_pv;idx_slack])) = [];
    
    dF = [dP;dQ]; % mismatch of the power flow equations
    
    %% Convergence check
    
    if(max(abs(dF))<tol)
        % disp('NR algorithm has converged to a solution!');
        break;
    elseif(k==n_max)
        disp('NR algorithm reached the maximum number of iterations!');
    end
    
    %% Jacobian construction
    
    % For the sake of simplicity, the blocks of J are constructed
    % for the whole network (i.e., with size n_nodes x n_nodes).
    % The unnecessary rows/columns are removed subsequently
    
    % Extract magnitude/angle
    V_abs = abs(V);
    V_arg = angle(V);
    
    % Initialization
    J_PE = zeros(n_variables,n_variables); % derivative: P versus E_abs
    J_PT = zeros(n_variables,n_variables); % derivative: P versus E_arg (theta)
    J_QE = zeros(n_variables,n_variables); % derivative: Q versus E_abs
    J_QT = zeros(n_variables,n_variables); % derivative: Q versus E_arg (theta)
    
    % Construction
    for i=1:n_variables
        
        % Diagonal elements (terms outside the sum)
        J_PE(i,i) =  2*Y_abs(i,i)*V_abs(i)*cos(Y_arg(i,i));
        J_QE(i,i) = -2*Y_abs(i,i)*V_abs(i)*sin(Y_arg(i,i));
        
        for j=1:n_variables
            if i ~= j
                % Diagonal elements (terms inside the sum)
                J_PE(i,i) = J_PE(i,i) + Y_abs(i,j)*V_abs(j)*cos(V_arg(i)-V_arg(j)-Y_arg(i,j));
                J_QE(i,i) = J_QE(i,i) + Y_abs(i,j)*V_abs(j)*sin(V_arg(i)-V_arg(j)-Y_arg(i,j));
                J_PT(i,i) = J_PT(i,i) - V_abs(i)*Y_abs(i,j)*V_abs(j)*sin(V_arg(i)-V_arg(j)-Y_arg(i,j));
                J_QT(i,i) = J_QT(i,i) + V_abs(i)*Y_abs(i,j)*V_abs(j)*cos(V_arg(i)-V_arg(j)-Y_arg(i,j));

                % Offdiagonal elements
                J_PE(i,j) =  Y_abs(i,j)*V_abs(i)*cos(V_arg(i)-V_arg(j)-Y_arg(i,j));
                J_QE(i,j) =  Y_abs(i,j)*V_abs(i)*sin(V_arg(i)-V_arg(j)-Y_arg(i,j));
                J_PT(i,j) =  Y_abs(i,j)*V_abs(i)*V_abs(j)*sin(V_arg(i)-V_arg(j)-Y_arg(i,j));
                J_QT(i,j) = -Y_abs(i,j)*V_abs(i)*V_abs(j)*cos(V_arg(i)-V_arg(j)-Y_arg(i,j));
            end
        end
    end
    
    % Remove extra rows (i.e., unnecessary equations)
    % slack bus: P & Q, PV buses: Q
    
    J_PE(idx_slack,:) = [];
    J_PT(idx_slack,:) = [];
    
    J_QE([idx_pv;idx_slack],:) = [];
    J_QT([idx_pv;idx_slack],:) = [];
    
    % Remove extra columns (i.e., variables)
    % slack bus: E_abs & E_arg, PV nodes: E_abs
    
    J_PE(:,[idx_pv,idx_slack]) = [];
    J_QE(:,[idx_pv,idx_slack]) = [];
    
    J_PT(:,idx_slack) = [];
    J_QT(:,idx_slack) = [];
    
    % Combination
    J = [J_PE,J_PT;J_QE,J_QT];
    
    %% Solution update
    
    % Solve
    dx = J \ dF;
    
    % Reconstruct the solution
    
    dV_abs = zeros(length(V_abs),1);
    dV_abs(idx_pq,1) = dx(1:length(idx_pq));
    
    dV_arg = zeros(length(V_arg),1);
    dV_arg(sort([idx_pq;idx_pv]),1) = dx((length(idx_pq)+1):end);
    
    % Update
    V_abs = V_abs + dV_abs;
    V_arg = V_arg + dV_arg;
end

V = V_abs .* exp(1i*V_arg);

end