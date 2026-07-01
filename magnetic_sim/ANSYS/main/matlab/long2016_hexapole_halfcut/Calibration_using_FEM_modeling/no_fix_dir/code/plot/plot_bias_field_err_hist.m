%% plot_bias_field_err_hist.m -- 18-param bias model vs FEM: vector-difference error histogram
% =========================================================================
%  橫軸 = 誤差 = |B_model - B_FEM|（兩場逐點逐激發的「向量差大小」），縱軸 = count。
%  資料 = gap200um_mueq；模型 = 18 參數 bias（ell + e_hat(17)，逐激發 profile 電荷 g_j）；6 激發全合併。
%  場模型 = A*g_j（g_j=(A'A)^-1 A'b），與 Hall_sensor_base_fix_dir/main_Dmatrix 的場重建相同。
%
%  REUSES（不重寫；同 no_fix_dir/main.m）：load_coils_actuator(+variant), select_ball, fit_bias, make_Pc, build_A。
%  風格：選項① 粗體框圖。輸出 -> 本組 figures/bias_field_err_hist_gap200um_mueq.png。
% =========================================================================

clear; clc;

%% ---- config ----------------------------------------------------------------
VARIANT  = 'gap200um_mueq';   % FEM 變體
R_select = 150e-6;            % 取點半徑 [m]
ell0     = 0.5e-3;           % ell_hat 初值 [m]（fit_bias 在 SI）
NB       = 40;                % 直方圖 bin 數

%% ---- paths（同 no_fix_dir/main.m）-----------------------------------------
TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\no_fix_dir'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');     % ansys_path
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % mt_constants/import_ansys_data/filter_iron_nodes
addpath(fullfile(TREE,'code','function'));
model  = 'long2016_hexapole_halfcut';
figdir = fullfile(TREE,'figures');
if ~exist(figdir,'dir'); mkdir(figdir); end

cnst = mt_constants();
apdl_to_paper_idx = [1,3,6,5,2,4];

%% ---- 讀結果防呆 + 載 6 coil（gap200, 'all'）------------------------------
fprintf('讀結果：coil{1..6}/%s（dataset=all）；期望 gap200 graded（~656k），|B| 較 baseline 低 ~30%%。\n', VARIANT);
D = load_coils_actuator(model, cnst, apdl_to_paper_idx, 'all', VARIANT);   % [沿用，+variant]
[P, Bstack, npts] = select_ball(D, R_select);                              % [沿用] Bstack=3Np×6（mT）；P [m]

%% ---- 擬合 + profile 電荷 → 模型場 → 殘差 → 逐點逐激發向量差大小 ------------
[ell_hat, e_hat, J] = fit_bias(P, Bstack, D.Pc_base, ell0);   % [沿用]（P/ell 均公尺，well-scaled）
Pc  = make_Pc(e_hat, D.Pc_base);                              % [沿用]
A   = build_A(ell_hat, Pc, P);                                % [沿用] 3Np×6
g   = (A.' * A) \ (A.' * Bstack);                             % 逐欄 profile 電荷 g_j（6×6）
Res = A*g - Bstack;                                           % 殘差 = 模型場 - FEM 場（3Np×6）
err = vecnorm(reshape(Res, 3, npts*6), 2, 1).';              % 逐點逐激發向量差大小 [mT]（Bstack 已 mT）

errpct = 100 * sqrt(sum(err.^2) / sum(Bstack(:).^2));
fprintf('ell_hat=%.2f µm | numel(err)=%d (=npts*6=%d) | region err=%.3f%% | median=%.4f mT | max=%.4f mT\n', ...
        ell_hat*1e6, numel(err), npts*6, errpct, median(err), max(err));

%% ---- 存 err 供疊圖腳本載入（避免 fix/no_fix 同名函式撞 path）-----------------
model_label = '18-param bias';                                    %#ok<NASGU>
datadir = fullfile(TREE,'data'); if ~exist(datadir,'dir'); mkdir(datadir); end
save(fullfile(datadir, sprintf('field_err_hist_%s.mat', VARIANT)), ...
     'err','errpct','VARIANT','model_label');

%% ---- 繪圖：選項① 粗體框圖 --------------------------------------------------
XMAX = 0.8; edges = 0:0.005:XMAX;                                      % 統一尺度 + 0.005 mT/格（更細→更 smooth）
fig = figure('Color','w','Position',[100 100 760 560]);
ax  = axes(fig);
histogram(ax, err, edges, 'FaceColor',[0.20 0.40 0.70], 'EdgeColor','w');  % err 已 mT

set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
box(ax,'on'); grid(ax,'off');
xlim(ax,[0 XMAX]);                                   % x=0 緊貼 y 軸
xt = get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));
yt = get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
xlabel(ax,'|B_{model} - B_{FEM}| (mT)','FontWeight','bold');
ylabel(ax,'Count','FontWeight','bold');

txt = sprintf('N = %d (6 exc.)\nmedian = %.3f mT', numel(err), median(err));
text(ax, 0.97, 0.95, txt, 'Units','normalized', 'HorizontalAlignment','right', ...
     'VerticalAlignment','top', 'FontSize',14, 'FontWeight','bold', ...
     'BackgroundColor','w', 'EdgeColor','k', 'LineWidth',1.5, 'Margin',5);

%% ---- 輸出實檔 -------------------------------------------------------------
outpng = fullfile(figdir, sprintf('bias_field_err_hist_%s.png', VARIANT));
exportgraphics(fig, outpng, 'Resolution',600);
fprintf('已輸出 %s\n', outpng);
