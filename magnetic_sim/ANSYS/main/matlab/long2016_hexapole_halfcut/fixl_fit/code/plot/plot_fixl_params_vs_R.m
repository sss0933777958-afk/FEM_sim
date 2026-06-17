%% plot_fixl_params_vs_R.m
%  ============================================================================
%  讀既有 sweep_alln_vs_R.mat(fix_l 全節點確定性 fit),重畫參數趨勢圖,
%  但**只留 ell_hat 與 gB_hat 兩個 panel,刪掉 ||K̂||_F 那格**。
%  樣式沿用原 sweep_alln_vs_R.m 的 3-panel 版(直立 y 標、藍圈/橘方、R* 標線)。
%  輸出覆寫原 sweep_alln_params_vs_R.png。
%  ============================================================================
clear; clc; close all;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');
model = 'long2016_hexapole_halfcut';

data_dir = matlab_path(model, 'charge_fit', 'fitting_trend');
fig_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
S = load(fullfile(data_dir,'sweep_alln_vs_R.mat'));         % R_um, ell_R, gB_R, R_best, ...

f = figure('Color','w','Position',[100 100 880 460]);
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact'); %#ok<NASGU>

nexttile;
plot(S.R_um, S.ell_R*1e3, '-o','LineWidth',1.4,'MarkerSize',3,'Color',[0.2 0.45 0.85]);
grid on; box on; set(gca,'FontSize',10);
ylabel('$\hat\ell$ [mm]','Interpreter','latex','FontSize',12);
title('Deterministic fit (all in-ball nodes) vs sampling radius $R$','Interpreter','latex','FontSize',12);

nexttile;
plot(S.R_um, S.gB_R*1e3, '-s','LineWidth',1.4,'MarkerSize',3,'Color',[0.85 0.4 0.1]);
grid on; box on; set(gca,'FontSize',10);
ylabel('$\hat g_B$ [$\times10^{-3}$]','Interpreter','latex','FontSize',12);
xlabel('sampling radius $R$ [$\mu$m]','Interpreter','latex','FontSize',12);

exportgraphics(f, fullfile(fig_dir,'sweep_alln_params_vs_R.png'),'Resolution',150);
fprintf('wrote %s (2-panel, K panel removed)\n', fullfile(fig_dir,'sweep_alln_params_vs_R.png'));
