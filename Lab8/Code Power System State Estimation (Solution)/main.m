% main script for linear power system state estimation

clear all;
close all;
clc;

disp(' ');
disp('************************************************************');
disp('* This script solves the state estimation problem *');
disp('************************************************************');
disp(' ');

%%  Step 1: Load nodal admittance matrix and base values
%   The data need to be saved in the folder 'Data_SE'
%   in the same directory as this script.

disp('STEP 1: loading admittance matrix Y_LF and base values A_b + V_b');
disp(' ');

% ! Indicate the number of phases !
n_phases = 1; % 1 (single-phase equivalent)

% ! Put the name of the file containing Y_LF !
Y_LF = importdata('Data_SE/Y_14bus.mat');

% ! Put the names of the files containing A_b and V_b !
A_b = importdata('Data_SE/Base_Power.mat');
V_b = importdata('Data_SE/Base_Voltage.mat');

n_system = size(Y_LF,1)/n_phases; % number of nodes (incl. slack)

disp(['The network has ' int2str(n_system) ' nodes (incl. slack).']);
disp(['The base power is ' num2str(A_b/1e6) ' MVA.']);
disp(['The base voltage is ' num2str(V_b/1e3) ' kV.']);
disp(' ');

Y_b = A_b/V_b^2; % base admittance

%%  Step 2: Construct the reduced nodal admittance matrix
%   Since the voltage of the slack node is fixed,
%   it does not need be estimated.

disp('STEP 2: constructing Y_SE ...');
disp(' ');

% Part A: equivalent circuit of the slack node

% ! Provide the short-circuit parameters !
S_sc = 300e6; % short-circuit power
r_sc = 0.1; % R/X ratio of the short-circuit impedance

X_sc =  V_b^2/(S_sc*sqrt(1+r_sc^2));
R_sc = r_sc*X_sc;
Z_sc = (R_sc + 1i*X_sc) * Y_b; % short-circuit impedance (absolute units)
Y_sc = 1/Z_sc;

% Part B: reduce the nodal admittance matrix

% ! Select the slack node (index) !
sel_LF_slack = 1;

% ! Select the node which is connected to the slack (index) !
sel_LF_neighbor = 2;

% Compute polyphase indices
idx_LF_slack = polyphase_indices(sel_LF_slack,n_phases);
idx_LF_neighbor = polyphase_indices(sel_LF_neighbor,n_phases);

% reduce the nodal admittance matrix
Y_SE = reduce_Y(Y_LF,Y_sc,idx_LF_slack,idx_LF_neighbor);

n_nodes = size(Y_SE,1)/n_phases; % number of nodes (excl. slack)

disp(['The network has ' num2str(n_nodes) ' nodes (excl. slack)']);
disp(' ');

%%  Step 3: Load nodal power profiles
%   The data need to be saved in the folder 'Data_SE'
%   in the same directory as this script.
%   The data need to be given as matrices of size n_nodes x n_timesteps.

disp('STEP 3: loading nodal power profiles ...')
disp(' ');

% ! Enter the names of the files with the absorption profiles !
P_abs = importdata('Data_SE/P_load.mat');
Q_abs = importdata('Data_SE/Q_load.mat');

% ! Enter the names of the files with the injection profiles !
P_inj = importdata('Data_SE/P_generation.mat');
Q_inj = importdata('Data_SE/Q_generation.mat');

n_timesteps = size(P_abs,2);

% Check number of nodes
if(~all([size(P_abs,1),size(Q_abs,1),size(P_inj,1),size(Q_inj,1)]==n_system))
    error('The nodal power profiles must have n_system rows.');
end

% Check number of timesteps
if(~all([size(P_abs,2),size(Q_abs,2),size(P_inj,2),size(Q_inj,2)]==n_timesteps))
    error('The nodal power profiles must have n_timesteps columns.');
end

% Manipulate the power profiles to create a suitable scenario
S_net = 0.5*complex(12*P_inj-P_abs,12*Q_inj-Q_abs)/A_b;

disp(['The profiles correspond to ' int2str(n_timesteps) ' timesteps.']);
disp(' ');

%%  Step 4: Configure the state estimator

disp('Step 4: Configuring SE ...');
disp(' ');

% ! Select which nodes to estimate (node indices) !
sel_SE = 1:n_nodes; % estimate entire network

% Create polyphase indices
idx_SE = polyphase_indices(sel_SE,n_phases);

n_states = 2 * length(idx_SE);

disp(['The SE has ' int2str(n_states) ' states.']);
disp(' ');

%%  Step 5: Construct the measurement model

disp('Step 5a: Configuring PMUs ...');
disp(' ');

% ! Select voltage measurement locations (node incices) !
sel_PMU_V = 1:n_nodes;

% ! Select current measurement locations (node indices) !
sel_PMU_I = 1:n_nodes;

% ! Indicate standard deviation of voltage measurement noise !
sd_PMU_V = 1e-3/3;

% ! Indicate standard deviation of current measurement noise !
sd_PMU_I = 1e-3/3;

% Create polyphase indices
idx_PMU_V = polyphase_indices(sel_PMU_V,n_phases);
idx_PMU_I = polyphase_indices(sel_PMU_I,n_phases);

n_measurements = 2*(length(idx_PMU_V)+length(idx_PMU_I));

disp(['The SE has ' int2str(n_measurements) ' measurements.']);
disp(' ');

% **********

disp('Step 5.2: Constructing H and R ...');
disp(' ');

H = construct_H(Y_SE,idx_PMU_V,idx_PMU_I,idx_SE);
R = construct_R(sd_PMU_V,sd_PMU_I,idx_PMU_V,idx_PMU_I);

R_inv = inv(R); % precompute the inverse for the WLS-SE

% **********

disp('Step 5.3: Verifying observability ...');
disp(' ');

if(n_states>n_measurements)
    error('The necessary observability criterion is violated.');
else
    disp('The necessary observability criterion is satisfied.');
end

if(rank(H)<min(size(H)))
    error('The sufficient observability criterion is violated.');
else
    disp('The sufficient observability criterion is satisifed.');
end

disp(' ');

%%  Step 6: Construct the process model

disp('Step 6: Constructing A, B, and Q ...');
disp(' ');

% ! Configure the process noise (co)variance values !
Q_var = 1e-8;
Q_cov = 0;

Q = construct_Q(Q_var,Q_cov,idx_SE);

[A,B] = construct_AB(idx_SE);

%% Step 7.1: Perform load flow

disp('Step 7.1: Solving load flow ...');
disp(' ');

% ! Select PQ-type nodes for load flow (node indices) !
sel_LF_pq = 2:n_system; % everything except the slack node

% ! Select PV-type nodes for load flow (node indices) !
sel_LF_pv = []; % node

% ! Enter the maximum number of iterations for the NR algorithm !
n_max = 100;

% ! Enter the tolerance for the NR algorithm !
tol = 1e-10;

% Create polyphase indices
idx_LF_slack = polyphase_indices(sel_LF_slack,n_phases);
idx_LF_pq = polyphase_indices(sel_LF_pq,n_phases);
idx_LF_pv = polyphase_indices(sel_LF_pv,n_phases);

% Solve load flow for true nodal voltage phasors

V_LF = zeros(n_system*n_phases,n_timesteps); % true nodal voltage phasors

for t=1:n_timesteps
    S_star = S_net(:,t);
    
    if(t<=1)
        V_0 = flat_start(n_system,n_phases);
    else
        V_0 = V_LF(:,t-1);
    end
    
    [V_LF(:,t),J,n] = do_LF(Y_LF,S_star,V_0,idx_LF_slack,idx_LF_pq,idx_LF_pv,tol,n_max);
end

V_LF(idx_LF_slack,:) = []; % remove slack node
I_LF = Y_SE * V_LF;

%% Step 7.2: Simulate phasor measurement units

disp('Step 7.2: Generating synchrophasor measurements ...');
disp(' ');

V_PMU = add_noise(V_LF(idx_PMU_V,:),0,sd_PMU_V);
I_PMU = add_noise(I_LF(idx_PMU_I,:),0,sd_PMU_I);

z = [real(V_PMU);imag(V_PMU);real(I_PMU);imag(I_PMU)];

%% Step 7.3: Perform state estimation

disp('Step 7.3: Performing state estimation ...');
disp(' ');

u = zeros(n_states,n_timesteps);

x_WLS = zeros(n_states,n_timesteps);
x_DKF = zeros(n_states,n_timesteps);

P_DKF = cell(1,n_timesteps);

t_DKF = zeros(1,n_timesteps);
t_WLS = zeros(1,n_timesteps);

for t=1:n_timesteps %hk is the time-step index
    % WLS
    
    t_start = tic;
    x_WLS(:,t) = estimate_WLS(z(:,t),H,R_inv);
    t_WLS(t) = toc(t_start);
    
    % DKF
    
    if(t<=1)
        % Initialization
        V_init = add_noise(V_LF(:,t),0,sd_PMU_V);
        
        x_DKF(:,t) = [real(V_init);imag(V_init)];
        P_DKF{t} = max(max(R))/max(max(Q)) * Q;
    else
        t_start = tic;
        [x_DKF(:,t),P_DKF{t}] = estimate_DKF(x_DKF(:,t-1),P_DKF{t-1},u(:,t-1),A,B,Q,z(:,t),H,R);
        t_DKF(t) = toc(t_start);
    end
end

%%  Step 8: Post-processing and visualization

disp('Step 8: Analysing and plotting ...');
disp(' ');

folder = pwd();

n_ignore = 100; % ignore startup phenomena (MATLAB, DKF, etc.)
idx_t = (n_ignore+1):n_timesteps;

V_LF = V_LF(:,idx_t);
V_PMU = V_PMU(:,idx_t);

z = z(:,idx_t);
u = u(:,idx_t);
x_WLS = x_WLS(:,idx_t);
x_DKF = x_DKF(:,idx_t);

t_WLS = t_WLS(idx_t);
t_DKF = t_DKF(idx_t);

%% Estimator Performance

% ! Select the PMU to be compared with the SE !
sel_e_PMU = 4;

idx_e_PMU = polyphase_indices(sel_e_PMU,n_phases);
idx_e_SE = idx_PMU_V(idx_e_PMU);

idx_x_re = 1:(n_states/2);
idx_x_im = (n_states/2+1):n_states;
V_WLS = complex(x_WLS(idx_x_re,:),x_WLS(idx_x_im,:));
V_DKF = complex(x_DKF(idx_x_re,:),x_DKF(idx_x_im,:));

h_V = figure(1);
clf;

subplot(2,1,1);

hold on;
plot(real(V_PMU(idx_e_PMU,:)),'-r');
plot(real(V_WLS(idx_e_SE,:)),'-b');
plot(real(V_DKF(idx_e_SE,:)),'-g');
plot(real(V_LF(idx_e_SE,:)),'-m');
hold off,

grid on;

title(['Voltage Phasor at Bus ' int2str(sel_PMU_V(sel_e_PMU)) ' (Real Part)']);
xlabel('Timestep');
ylabel('Value (pu)');
legend({'PMU','WLS','DKF','LF'},'Location','SouthEast');

subplot(2,1,2);

hold on;
plot(imag(V_PMU(idx_e_PMU,:)),'-r');
plot(imag(V_WLS(idx_e_SE,:)),'-b');
plot(imag(V_DKF(idx_e_SE,:)),'-g');
plot(imag(V_LF(idx_e_SE,:)),'-m');
hold off;

grid on;

title(['Voltage Phasor at Bus ' int2str(sel_PMU_V(sel_e_PMU)) ' (Imaginary Part)']);
xlabel('Timestep');
ylabel('Value (pu)');
legend({'PMU','WLS','DKF','LF'},'Location','SouthEast');

saveas(h_V,[folder filesep() 'Timeseries'],'epsc');

%% Estimation Error

x_LF = [real(V_LF);imag(V_LF)];

% WLS

e_WLS = x_WLS - x_LF;
e_DKF = x_DKF - x_LF;

h_e_WLS = plot_e(e_WLS,'WLS',2);
h_e_DKF = plot_e(e_DKF,'DKF',3);

saveas(h_e_WLS,[folder filesep() 'WLS_Estimation_Error'],'epsc');
saveas(h_e_DKF,[folder filesep() 'DKF_Estimation_Error'],'epsc');

% DKF

%% Inferred Measurement Noise

v_WLS = z - H*x_WLS;
v_DKF = z - H*x_DKF;

h_v_WLS = plot_v(v_WLS,'WLS',0,sd_PMU_V,4);
h_v_DKF = plot_v(v_DKF,'DKF',0,sd_PMU_V,5);

saveas(h_v_WLS,[folder filesep() 'WLS_Measurement_Noise'],'epsc');
saveas(h_v_DKF,[folder filesep() 'DKF_Measurement_Noise'],'epsc');

%% Inferred Process Noise

w_DKF = x_DKF(:,2:end) - (A*x_DKF(:,1:(end-1)) + B*u(:,1:(end-1)));

[w_center,f_actual,f_fitted] = get_distribution(w_DKF(:),100);
f_expected = pdf('Normal',w_center,0,sqrt(Q_var));

h_w_DKF = figure(6);
clf;

hold on;
bar(w_center,f_actual,'FaceColor','c','EdgeColor','none');
plot(w_center,f_fitted,'--b');
plot(w_center,f_expected,'--r');
hold off;

title('DKF: Inferred Process Noise Characteristics');
xlabel('Process Noise Variable');
ylabel('Probability Density Function');
labels = {};
labels{1} = 'Actual Distribution';
labels{2} = 'Normal Distribution (Fitted)';
labels{3} = 'Normal Distribution (Expected)';
legend(labels);

grid on;

saveas(h_w_DKF,[folder filesep() 'DKF_Process_Noise'],'epsc');

%% Execution Time

disp(['WLS: mean(t_exec)=' num2str(mean(t_WLS)*1000,'%4.3f') 'ms.']);
disp(['DKF: mean(t_exec)=' num2str(mean(t_DKF)*1000,'%4.3f') 'ms.']);
disp(' ');