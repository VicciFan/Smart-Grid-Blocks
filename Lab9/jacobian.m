function [J_invc] = jacobian(Y,V0)
n_nodes = length(V0)+1;
% Extract magnitude/angle
    E_abs = [1;abs(V0)];
    E_arg = [0;angle(V0)];
    Y_abs = abs(Y);
    Y_arg = angle(Y);
    % Initialization
    J_PE = zeros(n_nodes,n_nodes); % derivative: P versus E_abs
    J_PT = zeros(n_nodes,n_nodes); % derivative: P versus E_arg (theta)
    J_QE = zeros(n_nodes,n_nodes); % derivative: Q versus E_abs
    J_QT = zeros(n_nodes,n_nodes); % derivative: Q versus E_arg (theta)
    
    % Construction
for i=1:n_nodes 
        % Diagonal elements (terms outside the sum)
        J_PE(i,i) = 2*Y_abs(i,i)*E_abs(i)*cos(Y_arg(i,i)); % Complete
        J_QE(i,i) = -2*Y_abs(i,i)*E_abs(i)*sin(Y_arg(i,i)); % Complete
        
        for j=1:n_nodes
            if i ~= j
             % Diagonal elements (terms inside the sum)
                J_PE(i,i) = J_PE(i,i)+Y_abs(i,j)*E_abs(j)*cos(E_arg(i)-E_arg(j)-Y_arg(i,j)); % Complete
                J_QE(i,i) = J_QE(i,i)+Y_abs(i,j)*E_abs(j)*sin(E_arg(i)-E_arg(j)-Y_arg(i,j)); % Complete
                J_PT(i,i) = J_PT(i,i)-E_abs(i)*Y_abs(i,j)*E_abs(j)*sin(E_arg(i)-E_arg(j)-Y_arg(i,j)); % Complete
                J_QT(i,i) = J_QT(i,i)+E_abs(i)*Y_abs(i,j)*E_abs(j)*cos(E_arg(i)-E_arg(j)-Y_arg(i,j)); % Complete

                % Offdiagonal elements
                J_PE(i,j) = Y_abs(i,j)*E_abs(i)*cos(E_arg(i)-E_arg(j)-Y_arg(i,j)); % Complete
                J_QE(i,j) = Y_abs(i,j)*E_abs(i)*sin(E_arg(i)-E_arg(j)-Y_arg(i,j)); % Complete
                J_PT(i,j) = Y_abs(i,j)*E_abs(i)*E_abs(j)*sin(E_arg(i)-E_arg(j)-Y_arg(i,j)); % Complete
                J_QT(i,j) = -Y_abs(i,j)*E_abs(i)*E_abs(j)*cos(E_arg(i)-E_arg(j)-Y_arg(i,j)); % Complete
            end
        end
end
    
J = [J_PE(2:end,2:end),J_PT(2:end,2:end);J_QE(2:end,2:end),J_QT(2:end,2:end)];
J_inv = inv(J);
J_invc = J_inv(:,[4,8,16]);
    