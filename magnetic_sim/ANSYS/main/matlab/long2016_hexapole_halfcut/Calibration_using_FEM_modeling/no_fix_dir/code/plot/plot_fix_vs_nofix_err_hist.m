%% plot_fix_vs_nofix_err_hist.m -- 疊圖：fix-l vs 18-param bias 場誤差直方圖比較
% =========================================================================
%  把 fix_dir（fix-l 點電荷）與 no_fix_dir（18-param bias）兩模型 vs FEM（gap200）的
%  逐點逐激發向量差 |B_model-B_FEM| 直方圖**疊在同一張圖**比較。
%  資料來源 = 兩支直方圖腳本各自存的 err .mat（純 load，不 addpath 任何 code/function，
%  以避開 fix/no_fix 同名函式 select_ball/region_field_err 撞 path）。
%  風格：選項① 粗體框圖、半透明疊放、無 region_err、無 colorbar。
%  本副本輸出到 no_fix_dir/figures/（fix_dir/code/plot 有同檔輸出到 fix_dir/figures/）。
% =========================================================================

clear; clc;

VARIANT = 'gap200um_mueq';
CAL = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
       'long2016_hexapole_halfcut\Calibration_using_FEM_modeling'];
FIX_MAT   = fullfile(CAL,'fix_dir','data',    sprintf('field_err_hist_%s.mat', VARIANT));
NOFIX_MAT = fullfile(CAL,'no_fix_dir','data', sprintf('field_err_hist_%s.mat', VARIANT));
OUT_PNG   = fullfile(CAL,'no_fix_dir','figures', 'fix_vs_nofix_err_hist.png');   % <== 本副本輸出位置

%% ---- 載兩個 err（mT）------------------------------------------------------
F = load(FIX_MAT);   ef = F.err;   lf = sprintf('single-parameter  (median %.3f mT)', median(ef));   % err 已 mT
G = load(NOFIX_MAT); en = G.err;   ln = sprintf('18-param bias  (median %.3f mT)',     median(en));   % err 已 mT
fprintf('fix:   N=%d median=%.4f max=%.4f mT\n', numel(ef), median(ef), max(ef));
fprintf('nofix: N=%d median=%.4f max=%.4f mT\n', numel(en), median(en), max(en));

%% ---- 疊圖：統一尺度 + 0.01 mT/格（選項① 粗體框圖）-----------------------
XMAX = 0.8; edges = 0:0.005:XMAX;                 % 與單張統一尺度 + 0.005 mT/格（更細→更 smooth）
fig = figure('Color','w','Position',[100 100 820 580]);
ax  = axes(fig); hold(ax,'on');
histogram(ax, ef, edges, 'FaceColor',[0.85 0.33 0.10], 'FaceAlpha',0.55, 'EdgeColor','w');  % single-parameter（白邊線段，同單張直方圖）
histogram(ax, en, edges, 'FaceColor',[0.20 0.40 0.70], 'FaceAlpha',0.55, 'EdgeColor','w');  % bias

set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
box(ax,'on'); grid(ax,'off');
xlim(ax,[0 XMAX]);                                   % x=0 緊貼 y 軸
xt = get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));
yt = get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
xlabel(ax,'|B_{model} - B_{FEM}| (mT)','FontWeight','bold');
ylabel(ax,'Count','FontWeight','bold');
lg = legend(ax, {lf, ln}, 'Location','northeast','Interpreter','tex'); set(lg,'FontSize',14,'FontWeight','bold','Box','on');

exportgraphics(fig, OUT_PNG, 'Resolution',600);
fprintf('已輸出 %s\n', OUT_PNG);
