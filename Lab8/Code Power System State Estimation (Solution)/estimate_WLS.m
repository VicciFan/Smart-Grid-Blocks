function x = estimate_WLS(z,H,R_inv)
% x = estimate_WLS(z,H,R_inv)
%
% INPUT
% - z       measurement vector
% - H       state-to-measurement matrix of the measurement model
% - R_inv   inverse of the measurement noise covariance matrix
%
% OUTPUT
% - x       estimated state vector

G = H.' * R_inv * H;
x = G \ (H.' * R_inv * z);

end