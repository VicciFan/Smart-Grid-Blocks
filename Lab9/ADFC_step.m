clear all; close all; clc;
%% load data
AFC = importdata('AFC.mat');
Tmax = length(AFC.t); t = 1:1:Tmax; 
P = zeros(13,Tmax); ind = setdiff(1:13,[1;4;5;9;11]);
P(2,:) = AFC.P2t; P(3,:) = AFC.P3t; P(6,:) = AFC.P6t; P(7,:) = AFC.P7t;
P(8,:) = AFC.P8t; P(10,:) = AFC.P10t; P(12,:) = AFC.P12t; P(13,:) = AFC.P13t; 
Psum = sum(P(ind,:));
P4_mean = (AFC.P4t_min+AFC.P4t_max)/2; P11_mean = (AFC.P11t_min+AFC.P11t_max)/2;
P4_actual = AFC.P4t_actual; P11_actual = AFC.P11t_actual;
E9_last = AFC.E9_initial;
%% Construction of constant matrices
% Construction of A
A1 = [-1,-1,0]; A2 = -A1; A3 = [0,-1,0]; A4 = -A3;
A5 = [0,-0.1,0]; A6 = -A5; A7 = [0,-1,0]; A8 = -A7;
A9 = [0,-1,0]; A10 = -A7;
A = [A1;A2;A3;A4;A5;A6;A7;A8;A9;A10];
% Construction of b
b1 = [-1,-1]; b3 = [0,0]; b5 = [0,0]; 
b7 = [0,0]; b9 = [0,-1];
b2 = -b1; b4 = -b3;  b6 = -b5; b8 = -b7; b10 = -b9;
b = [b1';b2';b3';b4';b5';b6';b7';b8';b9';b10'];
% Construction of B
B1 = [0,0,-1;0,0,-1]; B2 = -B1; B3 = B1; B4 = -B1;
B5 = [0,0,-0.1;0,0,-0.1]; B6 = -B5; B7 = B1; B8 = -B1;
B9 = B1; B10 = -B1;
B = [B1;B2;B3;B4;B5;B6;B7;B8;B9;B10];
% Construction of D
D = [[-1;1;-1;1],zeros(4,1);zeros(4,1),[-1;1;-1;1]];
D_cell = repmat({D'}, 1, 10); D = blkdiag(D_cell{:});
% Construction of E
E = [-1, 0, 0; 1, 0, 0];
% Construction of e
e = [0;2];

%% Iterations
cost = zeros(Tmax,1); E9 = zeros(Tmax,1); x = zeros(Tmax,83);
P1 = zeros(Tmax,1); P5 = zeros(Tmax,1); P9 = zeros(Tmax,1);
for i=1:1:Tmax
% Construction of gamma
gamma1 = Psum(i)+1;        %1Ϊt
gamma2 = 1-Psum(i);        %1Ϊt
gamma3 = 1; gamma4 = 1;
gamma5 = 1-E9_last; gamma6 = E9_last;
gamma7 = P(8,i)+1; gamma8 = 1-P(8,i);  %1Ϊt
gamma9 = P(7,i)+P(12,i)+P(8,i)+P(10,i)+P(13,i)+1; gamma10 = 1-P(7,i)-P(12,i)-P(8,i)-P(10,i)-P(13,i);  %1Ϊt
gamma = [gamma1;gamma2;gamma3;gamma4;gamma5;gamma6;gamma7;gamma8;gamma9;gamma10];
% Construction of c
c = zeros(8,1);
c(1) = P(3,i)+1; c(2) = 1-P(3,i);        %1Ϊt
c(3) = -AFC.P4t_min(i); c(4) = AFC.P4t_max(i);    %1Ϊt
c(5) = P(10,i)+P(13,i)+1; c(6) = 1-P(10,i)-P(13,i);     %1Ϊt
c(7) = -AFC.P11t_min(i); c(8) = AFC.P11t_max(i);   %1Ϊt
C_cell = repmat({c'}, 1, 10); C = blkdiag(C_cell{:});
% Construction of f & h
f1 = -1; f2 = -2; f3 = -2*(P4_mean(i)+P11_mean(i));    %1Ϊt
f = [f1;f2;f3];
h = -3*Psum(i)-3*P4_mean(i)-3*P11_mean(i);      %1Ϊt

%% Opitimization
f = [f;zeros(80,1)];
A_linprog = [A,C;E,zeros(2,80)]; b_linprog = [gamma;e];
Aeq = [-B, D]; beq = b;
lb = [-inf*ones(3,1);zeros(80,1)]; ub = [];

[x(i,:),fval] = linprog(f,A_linprog,b_linprog,Aeq,beq,lb,ub);
cost(i) = fval + h; 
P1(i) = -Psum(i)-P4_actual(i)-P11_actual(i)-x(i,1)-(x(i,2)+x(i,3)*(P4_actual(i)+P11_actual(i)));
P5(i) = x(i,1);
P9(i) = x(i,2)+x(i,3)*(P4_actual(1)+P11_actual(1));
E9(i) = E9_last-0.1*P9(i);
E9_last = E9(i);
end

%% Plot
% Battery energy
plot(AFC.t,E9);
axis([10, 100, 0, 1])
xlabel('Time');
ylabel('SOE'); grid on;
% % Power
% plot(AFC.t,P5);hold on; plot(AFC.t,P9);hold on; 
% axis([10, 100, 0, 1])
% legend('P5','P9');
% xlabel('Time'); ylabel('Magnitude'); grid on;