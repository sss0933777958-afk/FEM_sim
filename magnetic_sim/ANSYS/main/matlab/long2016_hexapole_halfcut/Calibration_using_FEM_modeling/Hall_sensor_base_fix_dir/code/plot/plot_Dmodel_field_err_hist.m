%% plot_Dmodel_field_err_hist.m -- D-matrix model vs FEM field: vector-difference error histogram
% =========================================================================
%  使用者要求：一張直方圖，橫軸 = 誤差 = |B_model - B_FEM|（兩場逐點逐激發的「向量差大小」），
%  縱軸 = 落在該誤差的數量 (count)。資料 = gap200um_mueq，模型 = D 矩陣模型（V2），6 激發全合併。
%
%  模型場（V2 D-matrix，main_Dmatrix.m 步驟 1-3 的重建）：
%     b_ij = S(p_i, e_hat, ell_hat) * Dv_j ，Dv 由 (A'A)^-1 A' Bstack 解出。
%     模型場矩陣 = A*Dv（3Np x 6）；殘差 Res = A*Dv - Bstack。
%     逐點逐激發向量差大小 err = ||Res 的每 3 列一組|| → 共 Np*6 個誤差值，全合併成直方圖。
%
%  REUSES（不重寫；同 main_Dmatrix.m）：
%     no_fix_dir/code/function : load_coils_actuator(+variant), select_ball, fit_bias, make_Pc, build_A
%  風格：選項① 粗體框圖（rules/figure-style.md）。輸出 -> 本組 figures/Dmodel_field_err_hist.png。
% =========================================================================

clear; clc;

%% ---- config ----------------------------------------------------------------
VARIANT  = 'gap200um_mueq';   % FEM 變體：gap200 2 段式 μ_eff（graded mesh）
R_select = 150e-6;            % 取點半徑 [m]：擬合/誤差只用此球內 air 節點
ell0     = 0.5e-3;           % ell_hat 初值 [m]（= ell_design；fit_bias 在 SI）
NB       = 40;                % 直方圖 bin 數（定案前可改）

%% ---- paths（沿用既有函式，同 main_Dmatrix.m）------------------------------
CAL  = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % mt_constants/import_ansys_data/filter_iron_nodes
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');  % ansys_path
addpath(fullfile(CAL,'no_fix_dir','code','function'));  % load_coils_actuator/select_ball/fit_bias/make_Pc/build_A
model = 'long2016_hexapole_halfcut';
figdir = fullfile(CAL,'Hall_sensor_base_fix_dir','figures');
if ~exist(figdir,'dir'); mkdir(figdir); end

cnst = mt_constants();
apdl_to_paper_idx = [1,3,6,5,2,4];

%% ---- 讀結果防呆（result-read-safety 層①）----------------------------------
fprintf('讀結果：ANSYS_data/%s/data/coil{1..6}/%s（dataset=all）\n', model, VARIANT);
fprintf('  物理意義：gap200 支撐座 2 段式 μ_eff（graded mesh）。\n');
fprintf('  期望指紋：|B|max 比 baseline 低 ~30%%（~0.7 T）；npts = R%.0fµm 球內 air 節點數。\n', R_select*1e6);

%% ---- 步驟1：載 6-coil FEM（gap200, actuator 框）→ 選 R 球 ------------------
D = load_coils_actuator(model, cnst, apdl_to_paper_idx, 'all', VARIANT);   % [沿用，+variant]
[P, Bstack, npts] = select_ball(D, R_select);                              % [沿用] Bstack=3Np×6（mT）；P [m]
fprintf('選出 N_p = %d air 節點（R ≤ %g µm）；|B|max(球內) = %.2f mT\n', ...
        npts, R_select*1e6, max(vecnorm(reshape(Bstack,3,[]),2,1)));

%% ---- 步驟2：擬合 ell_hat、e_hat（profile Dv；lsqnonlin，在 SI 公尺）--------
[ell_hat, e_hat, J] = fit_bias(P, Bstack, D.Pc_base, ell0);                % [沿用]（P/ell 均公尺）
Pc = make_Pc(e_hat, D.Pc_base);                                            % [沿用]
fprintf('擬合：ell_hat = %.2f µm | J = %.6e\n', ell_hat*1e6, J);

%% ---- 步驟3：profile Dv → 模型場 A*Dv → 殘差 Res ----------------------------
A   = build_A(ell_hat, Pc, P);     % [沿用] 3Np×6 stacked kernel（同交錯堆疊）
Dv  = (A.' * A) \ (A.' * Bstack);  % Dv=(ΣS'S)^-1(ΣS'b)，6×6
Res = A*Dv - Bstack;               % 殘差 = 模型場 - FEM 場（3Np×6）

%% ---- 步驟4：逐點逐激發「向量差大小」err = ||[εx;εy;εz]|| -------------------
R3  = reshape(Res, 3, npts*6);     % 每欄 = 一點一激發的 3 分量殘差
err = vecnorm(R3, 2, 1).';         % (npts*6)×1，單位 mT（Bstack 已 mT）

% 驗證：region err% 應 ≈ main 印的數字
errpct = 100 * sqrt(sum(err.^2) / sum(Bstack(:).^2));
fprintf('驗證：numel(err) = %d（= npts*6 = %d）| region err = %.3f%% | median = %.4f mT | max = %.4f mT\n', ...
        numel(err), npts*6, errpct, median(err), max(err));

%% ---- 繪圖：選項① 粗體框圖 --------------------------------------------------
XMAX = 0.8; edges = 0:0.005:XMAX;                                      % 統一尺度 + 0.005 mT/格（更細→更 smooth）
fig = figure('Color','w','Position',[100 100 760 560]);
ax  = axes(fig);
histogram(ax, err, edges, 'FaceColor',[0.20 0.40 0.70], 'EdgeColor','w');  % err 已 mT

set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
box(ax,'on'); grid(ax,'off');
xlim(ax,[0 XMAX]);                                   % x=0 緊貼 y 軸
xt = get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));   % x tick 減半
yt = get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));   % y tick 減半
xlabel(ax,'|B_{model} - B_{FEM}| (mT)','FontWeight','bold');
ylabel(ax,'Count','FontWeight','bold');

% 註記（N / median）
txt = sprintf('N = %d (6 exc.)\nmedian = %.3f mT', numel(err), median(err));
text(ax, 0.97, 0.95, txt, 'Units','normalized', 'HorizontalAlignment','right', ...
     'VerticalAlignment','top', 'FontSize',14, 'FontWeight','bold', ...
     'BackgroundColor','w', 'EdgeColor','k', 'LineWidth',1.5, 'Margin',5);

%% ---- 輸出實檔（rules/figure-output.md）------------------------------------
outpng = fullfile(figdir, 'Dmodel_field_err_hist.png');
exportgraphics(fig, outpng, 'Resolution',600);
fprintf('已輸出 %s\n', outpng);
