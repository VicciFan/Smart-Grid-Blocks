clear all; close all; clc;
%% Read data
Y = importdata('Y13.mat');
VC = importdata('VC.mat');
V0 = VC.V; V_mag = abs(V0); V_arg = angle(V0); P = real(VC.S); Q = imag(VC.S);
Tmax = VC.NumTimeSlot; idx = [5,9]; E9_last = VC.E9_initial;
%% generate constant matrix (inequality constraint matrix)
A1 = [-eye(12),eye(12),zeros(12,15)];
A2 = [-eye(12),-eye(12),zeros(12,15)];
A3 = [zeros(12),-eye(12),zeros(12,15)]; A4 = -A3;
A5 = [zeros(1,36),-1,0,0]; A6 = -A5;
A7 = [zeros(1,36),0,0.1,0]; A8 = -A7;
A9 = [zeros(1,36),0,-1,0]; A10 = -A9;
A11 = [zeros(1,36),-0.2,0,1];
A12 = [zeros(1,36),-0.2,0,-1];
A = [A1;A2;A3;A4;A5;A6;A7;A8;A9;A10;A11;A12];
%% Iteration
P5 = zeros(Tmax,1); P9 = zeros(Tmax,1);
Q5 = zeros(Tmax,1); Q9 = zeros(Tmax,1);
V = zeros(Tmax,12); E9 = zeros(Tmax,1);
for t = 1:1:Tmax
    % Inequality constraint update
    b = [ones(12,1);-ones(12,1);-0.95*ones(12,1);1.05*ones(12,1);0;0.5;E9_last;1-E9_last;1;1;0;0];   %ÎªÉ¶update£¿
    % Equality constraint
    [Jinv_c] = jacobian(Y,V0(:,t));   %Jc is for P5,P9,Q5
    Aeq = zeros(24,39);
    Aeq(1:24,13:36) = -eye(24);
    Aeq(1:24,37:39) = Jinv_c;
    beq = Jinv_c*[P([4,8],t);Q(4,t)]-[V_mag(:,t);V_arg(:,t)];
    % Cost function
    f = [ones(12,1);zeros(27,1)];
    % Solve linear program
    [x,fval] = linprog(f,A,b,Aeq,beq);
    V(t,:) = x(13:24)';
    P5(t) = x(37); P9(t) = x(38); Q5(t) = x(39);
    E9_last = E9_last-0.1*P9(t);
    E9(t) = E9_last;
end
% update
V0 = V'; V_mag = abs(V0); V_arg = angle(V0);
P([4,8],:) = [P5';P9']; Q(4,:) = Q5'; E9_last = 0.5;                 %ÕâÊÇÉ¶£¿

% for t = 1:1:Tmax
%     % Inequality constraint update
%     b = [ones(12,1);-ones(12,1);-0.95*ones(12,1);1.05*ones(12,1);0;0.5;E9_last;1-E9_last;1;1;0;0];
%     % Equality constraint
%     [Jinv_c] = jacobian(Y,V0(:,t));   %Jc is for P5,P9,Q5
%     Aeq = zeros(24,39);
%     Aeq(1:24,13:36) = -eye(24);
%     Aeq(1:24,37:39) = Jinv_c;
%     beq = Jinv_c*[P([4,8],t);Q(4,t)]-[V_mag(:,t);V_arg(:,t)];
%     % Cost function
%     f = [ones(12,1);zeros(27,1)];
%     % Solve linear program
%     [x,fval] = linprog(f,A,b,Aeq,beq);
%     V(t,:) = x(13:24)';
%     P5(t) = x(37); P9(t) = x(38); Q5(t) = x(39);
%     E9_last = E9_last-0.1*P9(t);
%     E9(t) = E9_last;
% end
%% Plot
 plot(1:Tmax,P5); hold on; 
 plot(1:Tmax,P9); hold on;
% plot(1:Tmax,P(3,:),'g'); hold on;
% plot(1:Tmax,P(10,:),'y'); hold on;
% plot(1:Tmax,-sum(P)'+P(4,:)'+P(8,:)'-P5-P9,'m'); hold on;
% grid on;
 legend('P5','P9');
 xlabel('Time'); ylabel('Magnitude');

% bar(2:13,V(5,:)-V(4,:));

%bar(2:13,V(5,:)-abs(V0(:,5)'));
%xlabel('node index');
%ylabel('voltage difference');

