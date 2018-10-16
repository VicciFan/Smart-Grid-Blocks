function h_fig = plot_v(v,name,mu,sigma,id_fig)
% h_fig = plot_v(v,z_PMU,name,id_fig)
%
% INPUT
% - v       inferred measurement noise
% - name    state estimator name
% - mu      expected mean value
% - sigma   expected standard deviation
% - id_fig  figure number
% 
% OUTPUT
% - h_fig   figure handle

%%

h_fig = figure (id_fig);
clf;

[v_center,f_actual,f_fitted] = get_distribution(v(:),100);
f_expected = pdf('Normal',v_center,mu,sigma);

hold on;
bar(v_center,f_actual,'FaceColor','c','EdgeColor','none');
plot(v_center,f_fitted,'--b');
plot(v_center,f_expected,'--r');
hold off;

grid on;

title([name ': Inferred Measurement Noise Characteristics']);
xlabel('Measurement Noise Variable');
ylabel('Probability Density Function');
legend({'Actual Distribution','Normal Distribution (Fitted)','Normal Distribution (Expected)'});

end