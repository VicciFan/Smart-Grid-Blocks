clear all; close all; clc;
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
N = 48;

%% Simulation Setup
% Define initial state
x_init = x_ss;
T_init = 40;
% temperature comfort bounds
Tlweek = [repmat([15*ones(16,1);20*ones(20,1);15*ones(12,1)],5,1);15*ones(96+16,1);20*ones(20,1)];
Thweek = [repmat([30*ones(16,1);25*ones(20,1);30*ones(12,1)],5,1);30*ones(96+16,1);25*ones(20,1)];
Tav = 0.5*(Tlweek+Thweek);
% tariff scheme
tariffweek = price_Jan; 
%% Simulation
x = zeros(length(x_init),length(dp)-N); Tst = zeros(1,length(dp)-N);
x(:,1) = x_init; Tst(1) = T_init; 
e = zeros(3,length(dp)-N); est = zeros(1,length(dp)-N);
u = zeros(3,length(dp)-N); ust = zeros(1,length(dp)-N);
for i=1:length(dp)-N
    if i >= 2
        e(1,i) = (x(1,i)<(Tav(i)))||((x(1,i)<(Tav(i)+0.5))&&(x(1,i)>=x(1,i-1)));
        e(2,i) = (x(2,i)<(Tav(i)))||((x(2,i)<(Tav(i)+0.5))&&(x(2,i)>=x(2,i-1)));
        e(3,i) = (x(3,i)<(Tav(i)))||((x(3,i)<(Tav(i)+0.5))&&(x(3,i)>=x(3,i-1)));
        est(i) = (Tst(i)<(Tstmax-20))||(Tst(i)<Tstmax-10)&&(Tst(i)>Tst(i-1));
    end
    u(:,i) = e(:,i).*max_flow.*(Tst(i)*ones(3,1)-x(1:3,i));
    ust(i) = est(i)*stmax;
    x(:,i+1) = A*x(:,i)+Bu*u(:,i)+Bd*dplow(i,:)';
    Tst(i+1) = Ast*Tst(i)+Bust*[u(:,i);ust(i)]+Bdst*dplow(i,:)';
end    
%% Plot
%% Plotting
figure
plot(x(1,:));
grid on
hold on
plot(x(2,:),'r');
plot(x(3,:),'g');
plot(Tst(1,:),'k');
plot(Tlweek(1:size(Tst,2)),'--m');
plot(Thweek(1:size(Tst,2)),'--m');
legend('Temperature - ZN1','Temperature - ZN2','Temperature - ZN3','Temperature - Storage', 'comfort bounds')

figure
plot(u(1,:));
grid on
hold on
plot(u(2,:),'r');
plot(u(3,:),'g');
plot(ust,'k');
plot(u(1,:)+u(2,:)+u(3,:),'m');
plot(tariffweek(1:size(u,2)),'--m');
legend('Power Input - ZN1','Power Input - ZN2','Power Input - ZN3','Power Input - Storage', 'Total Power Input - ZN1,ZN2,ZN3')