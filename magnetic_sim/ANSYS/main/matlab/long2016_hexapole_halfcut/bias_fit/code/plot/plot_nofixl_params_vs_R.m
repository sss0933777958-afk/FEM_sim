%% plot_nofixl_params_vs_R.m
%  ============================================================================
%  讀 sweep_nofixl_vs_R.mat,畫 no_fix_l(18-param bias)模型的參數趨勢 vs 取樣半徑 R。
%  只兩個 panel:ell_hat 與 gB_hat(不畫 ||K̂||、不畫 NRMSE)。
%  y 軸標籤橫擺(Rotation 0),參照 sweep_alln_params_vs_R.png 的配色/標記但 2 格版。
%  ============================================================================
clear; clc; close all;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');
model = 'long2016_hexapole_halfcut';

data_dir = matlab_path(model, 'charge_fit', 'fitting_trend');
fig_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
S = load(fullfile(data_dir,'sweep_nofixl_vs_R.mat'));        % R_um, ell_R, gB_R, ...

f = figure('Color','w','Position',[100 100 880 440]);
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact'); %#ok<NASGU>

% --- panel 1: ell_hat [mm] ---
nexttile;
plot(S.R_um, S.ell_R*1e3, '-o','LineWidth',1.4,'MarkerSize',3,'Color',[0.2 0.45 0.85]);
grid on; box on; set(gca,'FontSize',10); xlim([45 505]);
yl1 = ylabel('$\hat\ell$ [mm]','Interpreter','latex','FontSize',13, ...
             'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle');
title('no\_fix\_l (18-param bias) fit, all in-ball nodes vs sampling radius $R$', ...
      'Interpreter','latex','FontSize',12);
set(gca,'XTickLabel',[]);                                    % 上格隱 x 數字(共用下格)

% --- panel 2: gB_hat [x1e-3] ---
nexttile;
plot(S.R_um, S.gB_R*1e3, '-s','LineWidth',1.4,'MarkerSize',3,'Color',[0.85 0.4 0.1]);
grid on; box on; set(gca,'FontSize',10); xlim([45 505]);
yl2 = ylabel('$\hat g_B$ [$\times10^{-3}$]','Interpreter','latex','FontSize',13, ...
             'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle');
xlabel('sampling radius $R$ [$\mu$m]','Interpreter','latex','FontSize',13);

drawnow;
exportgraphics(f, fullfile(fig_dir,'sweep_nofixl_params_vs_R.png'),'Resolution',150);
fprintf('wrote %s\n', fullfile(fig_dir,'sweep_nofixl_params_vs_R.png'));
