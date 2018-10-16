function H = construct_H(Y_SE,idx_PMU_V,idx_PMU_I,idx_SE)
% H = construct_H(Y_SE,idx_PMU_V,idx_PMU_I,idx_SE)
%
% INPUT
% - Y_SE        nodal admittance matrix (excl. slack node)
% - idx_PMU_V   polyphase indices of voltage measurement locations
% - idx_PMU_I   polyphase indices of current measurement locations
% - idx_SE      polyphase indices of estimated buses
%
% OUTPUT
% - H           measurement model matrix

%% Block for voltage measurements

U = eye(size(Y_SE));
H_V = blkdiag(U(idx_PMU_V,idx_SE),U(idx_PMU_V,idx_SE));

%% Block for current measurements

G = real(Y_SE);
B = imag(Y_SE);

H_I = cell(2,2);
H_I{1,1} =  G(idx_PMU_I,idx_SE);
H_I{1,2} = -B(idx_PMU_I,idx_SE);
H_I{2,1} =  B(idx_PMU_I,idx_SE);
H_I{2,2} =  G(idx_PMU_I,idx_SE);
H_I = cell2mat(H_I);

%% Combine the blocks

H = [H_V;H_I];

end