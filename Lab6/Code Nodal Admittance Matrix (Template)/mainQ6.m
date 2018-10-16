% admittance matrix calculus

clear all;
close all;

disp(' ');
disp('*************************************************');
disp('* This script computes the Y matrix *');
disp('*************************************************');
disp(' ');

%%  Step 1: line parameters

% Load the line parameters from a text file.
% The lines of the text file need to have the following format:
% Bus Bus R(Ohm/km) X(Ohm/km) B(Siemens/km) length(k).

disp('STEP 1: load line parameters');
disp(' ');

% !!!  put the correct file name here !!!
linedata = load('linedata.txt');

from = linedata(:,1);
to = linedata(:,2);
R_prime = linedata(:,3);
X_prime = linedata(:,4);
B_prime = linedata(:,5);
l = linedata(:,6);

n_lines = length(from);           % number of lines
n_nodes = max(max(from),max(to)); % number of nodes

disp(['The network consists of ', num2str(n_nodes) , ' nodes and ',num2str(n_lines), ' lines.']);
disp(' ');

%%  Step 2: base values & pi-equivalent circuits

% Define the base power and base voltage here.
% The base current and base admittance are computed automatically.

disp('STEP 2: transform line parameters to p.u.');
disp(' ');

% !!! put the correct base values here !!!
A_b = 5e6;      % base value for the power in VA
V_b = 4.16e3;   % base value for the line-to-line voltage in V

% I_b = A_b / (V_b.*sqrt(3));
Y_b = A_b / V_b^2;

% Construct the pi-section equivalent circuits, i.e. structs with fields
% i         start node
% j         end node
% Y_ij      branch admittance
% Y_i_ij    shunt admittance on start-node side
% Y_j_ij    shunt admittacne on end-node side

pi_circuits = cell(n_lines,1);

for k=1:n_lines    
    pi_k = struct();
    
    pi_k.i = from(k);
    pi_k.j = to(k);
    pi_k.Y_ij = 1 / ((R_prime(k) + 1i*X_prime(k)) * l(k)) / Y_b;
    pi_k.Y_i_ij = 1i/2 * B_prime(k) * l(k) / Y_b;
    pi_k.Y_j_ij = 1i/2 * B_prime(k) * l(k) / Y_b;
    
    pi_circuits{k} = pi_k;
end

pi_circuits = [pi_circuits{:}];
%%  Step 3: build primitive admittance matrices and incidence matrix

disp('STEP 3: build primitive admittance matrices and incidence matrix');
disp(' ');

[A,Y_L,Y_T] = calculate_parameters(pi_circuits);

%%  Step 4: compute the nodal admittance matrix

disp('STEP 4: compute the nodal admittance matrix')
disp(' ')

Y = A.' * Y_L * A + Y_T;

disp('admittance matrices (without regulator)');
disp(' ');

print_matrix(Y_L,'Y_L','%4.3e');
disp(' ');

print_matrix(Y_T,'Y_T','%4.3e');
disp(' ');

print_matrix(Y,'Y','%4.3e');
disp(' ');

%%  Step 5: include a regulating transformer (one-base approach)

disp('STEP 5: incorporate regulating transformer')
disp(' ')

A_reg = 1e7; % nominal power
V_reg = 5e3; % nominal voltage
t_reg = 1.2; % tap ratio
R_reg = 0.1; % per-unit short-circuit resistance
X_reg = 0.8; % per-unit short-circuit reactance
k_reg = 2;   % line index

Z_reg = (R_reg + 1i*X_reg) * (V_reg/V_b)^2 * (A_b/A_reg);
Y_reg = 1/Z_reg;

% one-base approach

% **********

    pi_circuits(2).Y_ij = t_reg*Y_reg;
    pi_circuits(2).Y_i_ij =t_reg*(t_reg-1)*Y_reg;
    pi_circuits(2).Y_j_ij = (1-t_reg)*Y_reg;

% **********

[A,Y_L,Y_T] = calculate_parameters(pi_circuits);

Y = A.' * Y_L * A + Y_T;

disp('admittance matrices (with regulator, one base)');
disp(' ');

print_matrix(Y_L,'Y_L','%4.3e');
disp(' ');

print_matrix(Y_T,'Y_T','%4.3e');
disp(' ');

print_matrix(Y,'Y','%4.3e');
disp(' ');

% two-base approach

% **********
    Y_b2 = Y_b/1.44;

k=2;
    pi_circuits(k).Y_ij = Y_reg*1.44;
    pi_circuits(k).Y_i_ij = 0;
    pi_circuits(k).Y_j_ij = 0;

  
    
for k=3:n_lines    

    pi_circuits(k).Y_ij = 1 / ((R_prime(k) + 1i*X_prime(k)) * l(k)) / Y_b2;
    pi_circuits(k).Y_i_ij = 1i/2 * B_prime(k) * l(k) / Y_b2;
    pi_circuits(k).Y_j_ij = 1i/2 * B_prime(k) * l(k) / Y_b2;
    
end


% **********

[A,Y_L,Y_T] = calculate_parameters(pi_circuits);

Y = A.' * Y_L * A + Y_T;

disp('admittance matrices (with regulator, two bases)');
disp(' ');

print_matrix(Y_L,'YL','%4.3e');
disp(' ');

print_matrix(Y_T,'YT','%4.3e');
disp(' ');

print_matrix(Y,'Y','%4.3e');
disp(' ');