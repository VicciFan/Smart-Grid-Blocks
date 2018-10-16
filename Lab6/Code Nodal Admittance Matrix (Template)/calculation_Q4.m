clear all; close all; clc;
% Base value and matrix define
S_B = 5e6; t = 1.2;
V_B1 = 4.16e3;
V_B2 = V_B1*t;
Y_T = zeros(6);
Y_L = zeros(7);

% A matrix
A = [1 -1 0 0 0 0
 0 1 -1 0 0 0
 0 0 1 -1 0 0
 0 0 1 0 -1 0
 0 0 1 0 0 -1
 0 0 0 1 -1 0
 0 0 0 1 0 -1];

% Line 1
S_sc = 300e6;
r_sc = 0.2;
Z_sc = V_B1^2/S_sc;
R1 = Z_sc*sin(atan(r_sc));
X1 = Z_sc*cos(atan(r_sc));
Z_pu(1) = R1*S_B/V_B1^2 + X1*S_B/V_B1^2*i;
Y_L(1,1) = 1/Z_pu(1);

% Transformer
Rt = 0.1; Xt = 0.8; Vt = 5e3; At = 10e6; t = 1.2;
Rt_pu = Rt * Vt^2/At * S_B/V_B2^2; Xt_pu = Xt * Vt^2/At * S_B/V_B2^2;
Yt_pu = 1/(Rt_pu+Xt_pu*1i);
Y_L(2,2) = Yt_pu;

% Other line parameters
Ru(3) = 0.069; Xu(3) = 0.208; Bu(3) = 1.115e-6; l(3) = 1.1;
Ru(4) = 0.173; Xu(4) = 0.312; Bu(4) = 1.205e-6; l(4) = 1.1;
Ru(5) = 0.173; Xu(5) = 0.312; Bu(5) = 1.205e-6; l(5) = 1.1;
Ru(6) = 0.115; Xu(6) = 0.331; Bu(6) = 1.224e-6; l(6) = 0.6;
Ru(7) = 0.115; Xu(7) = 0.331; Bu(7) = 1.224e-6; l(7) = 0.6;

for i = 3:1:7
    R(i) = Ru(i)*l(i); X(i) = Xu(i)*l(i); B(i) = Bu(i)*l(i);
    B_pu(i) = B(i)*V_B2^2/S_B*1i;
    R_pu(i) = R(i)*S_B/V_B2^2; X_pu(i) = X(i)*S_B/V_B2^2;
    Z_pu(i) = R_pu(i)+X_pu(i)*1i;
    Y_L(i,i) = 1/Z_pu(i);
end

% Shunt values
Y_T(3,3) = (B_pu(3)+B_pu(4)+B_pu(5))/2;
Y_T(4,4) = (B_pu(3)+B_pu(6)+B_pu(7))/2;
Y_T(5,5) = (B_pu(4)+B_pu(6))/2;
Y_T(6,6) = (B_pu(5)+B_pu(7))/2;

% Nodal admittance matrix
Y = A'*Y_L*A + Y_T;

Y_L = diag(Y_L);
Y_T = diag(Y_T);
