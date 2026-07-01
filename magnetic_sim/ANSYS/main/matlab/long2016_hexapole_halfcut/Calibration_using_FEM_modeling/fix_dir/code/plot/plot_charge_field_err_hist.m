%% plot_charge_field_err_hist.m -- fix-l point-charge model vs FEM: vector-difference error histogram
% =========================================================================
%  橫軸 = 誤差 = |B_model - B_FEM|（兩場逐點逐激發的「向量差大小」），縱軸 = count。
%  資料 = gap200um_mueq；模型 = fix-l 點電荷（電荷鎖在極軸 ell*dhat、無 bias）；6 激發全合併。
%
%  REUSES（不重寫；同 fix_dir/main.m）：load_coils(+variant), select_ball, fit_KI_fixl, charge_residual。
%  風格：選項① 粗體框圖。輸出 -> 本組 figures/charge_field_err_hist_gap200um_mueq.png。
% =========================================================================

clear; clc;

%% ---- config ----------------------------------------------------------------
VARIANT  = 'gap200um_mueq';   % FEM 變體
R_select = 150e-6;            % 取點半徑 [m]
I_actual = 1;                 % 驅動電流 [A] = FEM 激發
NB       = 40;                % 直方圖 bin 數

%% ---- paths（同 fix_dir/main.m）--------------------------------------------
TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\fix_dir'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % mt_constants/import_ansys_data/filter_iron_nodes
addpath(fullfile(TREE,'code','function'));
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
figdir = fullfile(TREE,'figures');
if ~exist(figdir,'dir'); mkdir(figdir); end

cnst = mt_constants();
apdl_to_paper_idx = [1,3,6,5,2,4];
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);

%% ---- 讀結果防呆 + 載 6 coil（gap200, 'wp'）-------------------------------
fprintf('讀結果：coil{1..6}/%s（dataset=wp）；期望 gap200 |B| 較 baseline 低 ~30%%。\n', VARIANT);
C = load_coils(results_root, cnst, apdl_to_paper_idx, VARIANT);   % [沿用，+variant]（bfem 已 mT）
[coil, nmin] = select_ball(C, R_select);                          % [沿用] coil(k).p [m], .bfem [mT]

%% ---- 擬合 fix-l 模型 → 殘差 → 逐點逐激發向量差大小 -------------------------
[ell, ghat_I_B, K_bar, J] = fit_KI_fixl(coil, dhat, I_actual);    % 在 SI 公尺擬合；^Bg_I [mT/A]、ell [m]
freemask = true(6); freemask(1,1) = false;                        % K_bar(1,1) 固定 5/6
x = [ell*1e3; ghat_I_B; K_bar(freemask)];                        % 重組 packed 參數（ell 打包成 mm，同 fit_KI_fixl 順序）
r = charge_residual(x, coil, dhat, I_actual, freemask);           % [沿用] 殘差（逐 coil 區塊串接；mT）

err = []; off = 0; sumB2 = 0;
for k = 1:numel(coil)
    Mk  = size(coil(k).p, 1);                                     % 該 coil 取點數（各異）
    blk = r(off + (1:3*Mk));  off = off + 3*Mk;                   % 該 coil 區塊 [Bx;By;Bz]
    ek  = vecnorm(reshape(blk, Mk, 3), 2, 2);                     % 逐點向量差大小 [mT]
    err = [err; ek]; %#ok<AGROW>
    sumB2 = sumB2 + sum(coil(k).bfem.^2);
end
errpct = 100 * sqrt(sum(err.^2) / sumB2);
fprintf('ell=%.2f µm | ^Bg_I=%.4e mT/A | numel(err)=%d | region err=%.3f%% | median=%.4f mT | max=%.4f mT\n', ...
        ell*1e6, ghat_I_B, numel(err), errpct, median(err), max(err));

%% ---- 存 err 供疊圖腳本載入（避免 fix/no_fix 同名函式撞 path）-----------------
model_label = 'fix-\ell';                                          %#ok<NASGU>
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
outpng = fullfile(figdir, sprintf('charge_field_err_hist_%s.png', VARIANT));
exportgraphics(fig, outpng, 'Resolution',600);
fprintf('已輸出 %s\n', outpng);
