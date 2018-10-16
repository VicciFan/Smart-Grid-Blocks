function X_noisy = add_noise(X,mu_N,sigma_N)
% X_noisy = add_noise(X,mu_N,sigma_N)
%
% INPUT
% - X           clean complex values
% - mu_N        mean of normally distributed noise
% - sigma_N     standard deviation of normally distributed noise
%
% OUTPUT
% - X_noisy     complex values with noise in real and imaginary part

[n_rows,n_cols] = size(X);

N_re = mu_N*ones(n_rows,n_cols) + sigma_N*randn(n_rows,n_cols);
N_im = mu_N*ones(n_rows,n_cols) + sigma_N*randn(n_rows,n_cols);

N = N_re + 1i*N_im;

X_noisy = X + N;

end