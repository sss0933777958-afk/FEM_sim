%% plot_combo_nrmse_overlay.m — overlay the 5 combo NRMSE curves on one axes
clear; clc; close all;
data_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit\fitting_trend';
fig_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
S = load(fullfile(data_dir,'combo_nrmse.mat'));   % R, NRMSE(5xnR), combos, labels, N

cols = lines(5);
f = figure('Color','w','Position',[100 100 820 520]); hold on;
leg = strings(1,5);
for c = 1:5
    plot(S.R, S.NRMSE(c,:), '-o','LineWidth',1.4,'MarkerSize',3.5, ...
         'Color',cols(c,:), 'MarkerFaceColor',cols(c,:));
    leg(c) = sprintf('%s: [%s]', S.labels{c}, strjoin(string(S.combos(c,:)),' '));
end
grid on; box on; set(gca,'FontSize',11); xlim([0 500]);
xlabel('Sampling radius  $R$  [$\mu$m]  (params from)','Interpreter','latex','FontSize',13);
ylabel('$\mathrm{NRMSE}_{\max}$  [\%]','Interpreter','latex','FontSize',13);
title(sprintf('Combo NRMSE on $R\\le50\\,\\mu$m ball ($N=%d$), $I$ in A', S.N), ...
      'Interpreter','latex','FontSize',12);
legend(leg, 'Location','northwest', 'Interpreter','none', 'FontSize',10);
exportgraphics(f, fullfile(fig_dir,'combo_nrmse_overlay.png'), 'Resolution',150);
fprintf('wrote combo_nrmse_overlay.png\n');
