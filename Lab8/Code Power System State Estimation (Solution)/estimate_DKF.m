function [x_new,P_new] = estimate_DKF(x_old,P_old,u,A,B,Q,z,H,R)
% [x_new,P_new] = estimate_DKF(x_old,P_old,A,B,Q,z,H,R)
%
% INPUT
% - x_old   a-posteriori estimated state vector at t-1
% - P_old   a-posteriori estimation error covariance matrix at t-1
% - u       control vector at t-1
% - A       state-to-state matrix of the process model
% - B       control-to-state matrix of the process model
% - Q       process noise covariance matrix
% - z       measurement vector at t-1
% - H       state-to-measurement matrix of the measurement model
% - R       measurement noise covariance matrix
% 
% OUTPUT
% - x_new   a-posteriori estimated state vector at t
% - P_new   a-posteriori estimation error covariance matrix at t

%% A-Priori Estimation (Prediction)

x_priori = A*x_old + B*u;
P_priori = A*P_old*A.' + Q;

%% A-Posteriori Estimation (Estimation)

U = eye(size(P_old));

K = P_priori * H.' / (H*P_priori*H.' + R);
x_posteriori = x_priori + K*(z-H*x_priori);
P_posteriori = (U-K*H) * P_priori;

%% Return

x_new = x_posteriori;
P_new = P_posteriori;

end