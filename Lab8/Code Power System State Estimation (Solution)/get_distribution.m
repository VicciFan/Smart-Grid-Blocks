function [x_center,f_actual,f_fitted] = get_distribution(x,n_bin)
% [x_center,f_actual,f_fitted] = get_distribution(x,n_bin)
%
% INPUT
% - x           (supposedly) randomly distributed variable
% - n           number of bins for assessing the distribution
%
% OUTPUT
% - x_center    center values of the bins
% - f_actual    actual probability density function
% - f_fitted    fitted probability density function (normal distribution)

x = x(:);

mu = mean(x);
sigma = std(x);

x_center = 4 * sigma * linspace(-1,1,n_bin);
dx = x_center(2)-x_center(1);
x_count = hist(x,x_center);

f_actual = hist(x,x_center) / (sum(x_count)*dx);
f_fitted = 1/sqrt(2*pi*sigma^2)*exp(-(x_center-mu).^2/(2*sigma^2));

end