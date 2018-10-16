clc;
clear all;
close all;

mkdir('YALMIP-master');
addpath(genpath('YALMIP-master'));
%% Data setup
% data time resolution: 30 min

%Model of the building
load bM; % building model data
A = bM.A;
Bu = bM.Bu;
Bd = bM.Bd;

%Model of the storage
load storage; % storage model data
Ast = storage.A;
Bust = storage.Bu;
Bdst = storage.Bd;

%Disturbance
load disturbance.mat;  % Including temperature data 
% Three vaiables are loaded:
% dp - disturbance prediction
% dpup - disturbance upper bound
% dplow - disturbance lower bound

%Prices
load prices.mat % price data 
% Includes the price data variable price_Jan

%Constraints
stmax = 150;             %Maximum power provided to the storage inkW

max_flow = 0.7*[1;0.3;1];  %Instantaneous time-flow in the three zones

Tstmax = 50;            %Maximum allowed temperature of the storage


%% MPC controller Design

N = 48;     % Prediction Horizon in timesteps (24 hours)

% YALMIP symbolic variables definition
x0 = sdpvar(nS,1,'full'); % Initial state 
Tst0 = sdpvar(1,1,'full'); % Initial storage tate 
d = sdpvar(nd,N,'full');    %Disturbance
dup = sdpvar(nd,N,'full');    %Disturbance upper bound
dlow = sdpvar(nd,N,'full');    %Disturbance lower bound
u = sdpvar(nu,N,'full');    %Input to the TABS
ep = sdpvar(nu,N,'full');   %Amount of constraint violations
st = sdpvar(1,N,'full');    %Input to the thermal storage
tf = sdpvar(1,N,'full');    %Tariff scheme (price)
Tlow = sdpvar(1,N,'full');    %Temperature lower bound
Thigh = sdpvar(1,N,'full');    %Temperature upper bound
Sep = sdpvar(1,1,'full'); % violation objective weight
Sinp = sdpvar(1,1,'full'); % price objective weight

normInp = diag(1./stmax);

%% non robust optimizer
x=x0;
Tst = Tst0;
constraints = [];
objective = 0;
for k=1:N
    % define here, iteratively, constraints, state evolution, and objective
    % function, using symbolic variables
    constraints = [constraints, 0 <= u(:,k) <= max_flow.*(Tst*ones(3,1)-x(1:nu))];
    
    x = A*x + Bu*u(:,k) + Bd*d(:,k);
    Tst = Ast*Tst + Bust*[u(:,k);st(k)] + Bdst*d(:,k);
    constraints = [constraints, Tst <= Tstmax; ep(:,k) >= Tlow(k)*ones(3,1)-x(1:nu); ep(:,k)>= x(1:nu)-Thigh(k)*ones(3,1)];

    objective = objective  +  Sep*norm(ep(:,k),1) + Sinp*normInp*tf(:,k)*st(:,k);
end
constraints = [constraints, ep(:)>=0, 0<=st(:)<=stmax];

 options = sdpsettings('verbose',1,'solver','linprog');
cont = optimizer(constraints, objective, options, [x0; Tst0; d(:); tf(:); Tlow(:); Thigh(:); Sep; Sinp], [u; st]);        %Definition of the controller for the building

%% robust optimizer
xup = x0;
xlow = x0;
Tstup = Tst0;
constraints = [];
objective = 0;

for k=1:N
    % define here, iteratively, constraints, state evolution, and objective
    % function, using symbolic variables
    constraints = [constraints, 0 <= u(:,k) <= max_flow.*(Tstup*ones(3,1)-xup(1:nu))];
    
    xup = A*xup + Bu*u(:,k) + Bd*dup(:,k);
    xlow = A*xlow + Bu*u(:,k) + Bd*dlow(:,k);
    Tstup = Ast*Tstup + Bust*[u(:,k);st(k)] + Bdst*dup(:,k);
    
    constraints = [constraints, Tstup <= Tstmax; ep(:,k) >= Tlow(k)*ones(3,1)-xlow(1:nu); ep(:,k)>= xup(1:nu)-Thigh(k)*ones(3,1)];

    objective = objective  +  Sep*norm(ep(:,k),1) + Sinp*normInp*tf(:,k)*st(:,k);
    
end
constraints = [constraints, ep(:)>=0, 0<=st(:)<=stmax];

contRob = optimizer(constraints, objective, options, [x0; Tst0; dup(:); dlow(:); tf(:); Tlow(:); Thigh(:); Sep; Sinp], [u; st]);        %Definition of the controller for the building


%% Simulation setup

% temperature comfort bounds
Tlweek = [repmat([15*ones(16,1);20*ones(20,1);15*ones(12,1)],5,1);15*ones(96+16,1);20*ones(20,1)];
Thweek = [repmat([30*ones(16,1);25*ones(20,1);30*ones(12,1)],5,1);30*ones(96+16,1);25*ones(20,1)];

% tariff scheme
tariffweek = price_Jan; 

% use this for Q1-Q10
% Wep = 10; % violation objective weight
% Winp = 1; % price objective weight

% use this for Q11
Wep = 100; % violation objective weight
Winp = 1; % price objective weight

%% Non-robust simulation

x_init = x_ss;
T_init = 40;

StateHistory = [];
InputHistory = [];
PredHistory = [];
StorageHistory = [];
costNonRob = 0; % total accumulated cost
for i=1:length(dp)-N    
    % ***************************************************
    % simulate here the non-robust MPC policy as required
    % ***************************************************
    pred = dp(i:i+N-1,:)';
    tariff =  tariffweek(i+1:i+N);
    Tl = Tlweek(i+1:i+N);
    Th = Thweek(i+1:i+N);
    [output,infeasible] = cont{[x_init;T_init;pred(:);tariff;Tl;Th;Wep;Winp]};
    u_cont = output(1:nu,:);
    u_st = output(nu+1,:);
    costNonRob = costNonRob + u_st(1)*tariff(1);

    % actual disturbance
    dp_true = dplow(i, :)';

    % Complete the code here...
    xplus = A*x_init+Bu*u_cont(:,1)+Bd*dp_true;
    T_plus = Ast*T_init + Bust*[u_cont(:,1);u_st(1)] + Bdst*dp_true;
    
    x_init = xplus; 
    T_init = T_plus;
     
    % Update the history
    StateHistory = [StateHistory, xplus];
    InputHistory = [InputHistory, [u_cont(:,1); u_st(:,1)]];
    PredHistory = [PredHistory, pred(:,1)];
    StorageHistory = [StorageHistory,T_plus];
end


%% Plotting
figure
plot(StateHistory(1,:));
grid on
hold on
plot(StateHistory(2,:),'r');
plot(StateHistory(3,:),'g');
plot(StorageHistory(1,:),'k');
plot(Tlweek(1:size(StorageHistory,2)),'--m');
plot(Thweek(1:size(StorageHistory,2)),'--m');
legend('Temperature - ZN1','Temperature - ZN2','Temperature - ZN3','Temperature - Storage', 'comfort bounds')

figure
plot(InputHistory(1,:));
grid on
hold on
plot(InputHistory(2,:),'r');
plot(InputHistory(3,:),'g');
plot(InputHistory(4,:),'k');
plot(InputHistory(1,:)+InputHistory(2,:)+InputHistory(3,:),'m');
plot(tariffweek(1:size(InputHistory,2)),'--m');
legend('Power Input - ZN1','Power Input - ZN2','Power Input - ZN3','Power Input - Storage', 'Total Power Input - ZN1,ZN2,ZN3')

%% Thermostat controller simulation


x_init = x_ss;
T_init = 40;

StateHistory = [];
InputHistory = [];
StorageHistory = [];
costThermo = 0; % total accumulated cost

for i=1:length(dp)-N    
    % ***************************************************
    % simulate here the thermostat controller as required
    % ***************************************************
    
     Tmean = 0.5*(Tlweek(i+1)+Thweek(i+1));
     
     e = [0,0,0]; 
     est = 0;
     
     u = [0,0,0]; 
     ust = 0;
  
    
    e(1) = (x_init(1)<Tmean)||((x_init(1)<(Tmean+0.5))&&(deltax(1)>0));
    u(1) = e(1)*max_flow(1)*(T_init-x_init(1));
    
    e(2) = (x_init(2)<Tmean)||((x_init(2)<(Tmean+0.5))&&(deltax(2)>0));
    u(2) = e(2)*max_flow(2)*(T_init-x_init(2));
    
    e(3) = (x_init(3)<Tmean)||((x_init(3)<(Tmean+0.5))&&(deltax(3)>0));
    u(3) = e(3)*max_flow(3)*(T_init-x_init(3));

    est = (T_init<(Tstmax-20))||((T_init<Tstmax-10)&&(deltaT>0));
    ust = est*stmax;

    costThermo = costThermo + ust*tariffweek(i+1);

    % actual disturbance
    dp_true = dplow(i, :)';
    
    % Complete the code here...
    xplus = A*x_init+Bu*u(:)+Bd*dp_true;
    T_plus = Ast*T_init+Bust*[u(:);ust]+Bdst*dp_true;
    deltax = xplus - x_init;
    deltaT = T_plus - T_init;
    x_init = xplus;
    T_init = T_plus;
      
    % Update the history
    StateHistory = [StateHistory, xplus];
    InputHistory = [InputHistory, [u(:); ust(:,1)]];
    StorageHistory = [StorageHistory,T_plus];
end



%% Plotting
figure
plot(StateHistory(1,:));
grid on
hold on
plot(StateHistory(2,:),'r');
plot(StateHistory(3,:),'g');
plot(StorageHistory(1,:),'k');
plot(Tlweek(1:size(StorageHistory,2)),'--m');
plot(Thweek(1:size(StorageHistory,2)),'--m');
legend('Temperature - ZN1','Temperature - ZN2','Temperature - ZN3','Temperature - Storage', 'comfort bounds')

figure
plot(InputHistory(1,:));
grid on
hold on
plot(InputHistory(2,:),'r');
plot(InputHistory(3,:),'g');
plot(InputHistory(4,:),'k');
plot(InputHistory(1,:)+InputHistory(2,:)+InputHistory(3,:),'m');
plot(tariffweek(1:size(InputHistory,2)),'--m');
legend('Power Input - ZN1','Power Input - ZN2','Power Input - ZN3','Power Input - Storage', 'Total Power Input - ZN1,ZN2,ZN3')


%% Robust simulation
x_init = x_ss;
T_init = 40;

StateHistory = [];
InputHistory = [];
PredHistory = [];
StorageHistory = [];
costRob = 0;
for i=1:length(dp)-N    
    % ***************************************************
    % simulate here the robust MPC policy as required
    % ***************************************************
    pred = dp(i:i+N-1,:)';
    predUp = dpup(i:i+N-1,:)';
    predLow = dplow(i:i+N-1,:)';

    tariff =  tariffweek(i+1:i+N);
    Tl = Tlweek(i+1:i+N);
    Th = Thweek(i+1:i+N);
    [output,infeasible] = contRob{[x_init;T_init;predUp(:);predLow(:);tariff;Tl;Th;Wep;Winp]};  
    
    u_contRob = output(1:nu,:);
    u_st = output(nu+1,:);
    costRob = costRob + u_st(1)*tariff(1);
    
    % actual disturbance
    dp_true = dplow(i, :)';

    % complete...
    xplus = A*x_init+Bu*u_contRob(:,1)+Bd*dp_true;
    T_plus = Ast*T_init+Bust*[u_contRob(:,1);u_st(1)]+Bdst*dp_true;
    x_init = xplus; T_init = T_plus;

    
    StateHistory = [StateHistory, xplus];
    InputHistory = [InputHistory, [u_contRob(:,1); u_st(:,1)]];
    PredHistory = [PredHistory, pred(:,1)];
    StorageHistory = [StorageHistory,T_plus];
end

%% Plotting
figure
plot(StateHistory(1,:));
grid on
hold on
plot(StateHistory(2,:),'r');
plot(StateHistory(3,:),'g');
plot(StorageHistory(1,:),'k');
plot(Tlweek(1:size(StorageHistory,2)),'--m');
plot(Thweek(1:size(StorageHistory,2)),'--m');
legend('Temperature - ZN1','Temperature - ZN2','Temperature - ZN3','Temperature - Storage', 'comfort bounds')

figure
plot(InputHistory(1,:));
grid on
hold on
plot(InputHistory(2,:),'r');
plot(InputHistory(3,:),'g');
plot(InputHistory(4,:),'k');
plot(InputHistory(1,:)+InputHistory(2,:)+InputHistory(3,:),'m');
plot(tariffweek(1:size(InputHistory,2)),'--m');
legend('Power Input - ZN1','Power Input - ZN2','Power Input - ZN3','Power Input - Storage', 'Total Power Input - ZN1,ZN2,ZN3')

%% Putting more weight on violation...
Wep = 100; % violation objective weight
Winp = 1; % price objective weight

%% Non-robust simulation

% complete...

%% Robust simulation

% complete...
