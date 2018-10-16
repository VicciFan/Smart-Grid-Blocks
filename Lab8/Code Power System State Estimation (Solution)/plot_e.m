function h_fig = plot_e(e,name,id_fig)
% h_fig = plot_e(e,name,id_fig)
%
% INPUT
% - e       estimation error
% - name    state estimator name
% - id_fig  figure number
% 
% OUTPUT
% - h_fig   figure handle

h_fig = figure (id_fig);
clf;

[e_center,f_actual,f_fitted] = get_distribution(e(:),100);

hold on;
bar(e_center,f_actual,'FaceColor','c','EdgeColor','none');
plot(e_center,f_fitted,'--b');
hold off;

grid on;

title([name ': Estimation Error Characteristics']);
xlabel('Estimation Error Variable');
ylabel('Probability Density Function');
legend({'Actual Distribution','Normal Distribution (Fitted)'});

end