function plot_ell_gain_vs_R(USE_BIAS)
% plot_ell_gain_vs_R -- long2016 半切六極 paper 圖:ℓ̂ 與 ĝ_I 隨取樣半徑 R 的趨勢
% =========================================================================
%   graded FEM 資料 + current_base 校正(fix 單一 ℓ,USE_BIAS=false)。
%   掃描取樣半徑 R = 40~500 µm,每點:select_ball → fitting(fix) → solve_current。
%   上 panel = ℓ̂(R) [µm]、下 panel = ĝ_I(R) [mT/A](尺度差大 → 分兩 panel)。
%   前段(讀 .dat + build_actuator_data)只做一次;loop R 只重跑 select_ball/fitting/solve_current。
%   風格 = 選項①粗體框 @ paper scale(font 粗體、LineWidth3、box on、grid off、tick 減半、單位 ())。
%   輸出 → figures/paper_fig/Section2_E/ell_gain_vs_R.png(覆蓋迭代)。
% =========================================================================
    clc;
    if nargin < 1, USE_BIAS = false; end   % false = fix → ell_gain_vs_R.png；true = 18-param bias → ell_gain_vs_R_bias.png
    here   = fileparts(mfilename('fullpath'));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CAL,'function'));  addpath(fullfile(CAL,'common_path'));

    % ---- 前段 pipeline(只做一次):graded / current / actuator frame ----
    cfg = model_config('long2016_hexapole_halfcut','tip40um');
    raw = extract_ansys_data(cfg, 'all', 'graded');     % 印 Matched N(指紋核對)
    ad  = build_actuator_data(raw, cfg);   Pc_base = ad.Pc_base;
    F   = zeros(6, cfg.N_I);
    for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end

    % ---- 掃描 R(fix 模型;每 20µm 一點;從 ĝ_I 正值起 R=40µm;R≤20 病態非物理)----
    Rum = 40:20:500;                                    % 取樣半徑 [µm](固定 20µm 間距)
    nR  = numel(Rum);
    ell_R = nan(1,nR);  gI_R = nan(1,nR);  npts_R = zeros(1,nR);
    fprintf('\n  R[um]   npts    ell_hat[um]   gI_hat[mT/A]\n');
    for i = 1:nR
        [P, Bstack, npts] = cfg.select_ball(ad, Rum(i)*1e-6);
        npts_R(i) = npts;
        if npts < 3, fprintf('  %5d  %6d   (skip: npts<3)\n', Rum(i), npts); continue; end
        [e, l_hat, ~] = fitting(P, Bstack, Pc_base, 0.5e-3, USE_BIAS);   % fix / bias 由旗標切
        [~, gI_hat]   = solve_current(l_hat, e, Pc_base, P, Bstack, F);
        ell_R(i) = l_hat*1e6;  gI_R(i) = gI_hat;
        fprintf('  %5d  %6d   %10.2f   %11.4f\n', Rum(i), npts, ell_R(i), gI_R(i));
    end

    % ---- 畫圖(2×1 panel,選項①粗體框)----
    cL = [0.05 0.10 0.95];   cR = [0.85 0.10 0.10];     % ℓ̂ 亮藍 / ĝ_I 紅
    FS = 30;                                            % 對齊 circuit 圖字體
    fig = figure('Color','w','Position',[100 80 980 940]);
    t = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

    % 上:ℓ̂(R)
    ax1 = nexttile(t);  hold(ax1,'on');
    plot(ax1, Rum, ell_R, '-o', 'Color',cL, 'LineWidth',3.0, 'MarkerSize',8, 'MarkerFaceColor',cL);
    style_panel(ax1, FS);
    xlim(ax1,[40 500]);   set(ax1,'XTick',[40 100 200 300 400 500]);   % 含第一點 40 + 終點 500
    if ~USE_BIAS
        ylim(ax1,[750 870]);  set(ax1,'YTick',[770 800 830 860]);   % fix:4 內縮 tick(不含端點)
    else
        ylim(ax1,[750 870]);  set(ax1,'YTick',[770 800 830 860]);   % bias:ℓ̂ 765.8~864.4，4 內縮 tick
    end
    ax1.XTickLabel = {};                                % 上 panel x 數字隱藏(共用軸);無軸標題
    hold(ax1,'off');

    % 下:ĝ_I(R)
    ax2 = nexttile(t);  hold(ax2,'on');
    plot(ax2, Rum, gI_R, '-s', 'Color',cR, 'LineWidth',3.0, 'MarkerSize',8, 'MarkerFaceColor',cR);
    style_panel(ax2, FS);
    xlim(ax2,[40 500]);   set(ax2,'XTick',[40 100 200 300 400 500]);   % 含第一點 40 + 終點 500
    if ~USE_BIAS
        ylim(ax2,[7.7 8.5]);  set(ax2,'YTick',[7.8 8.0 8.2 8.4]);   % fix:4 內縮 tick(不含端點)
    else
        ylim(ax2,[9.5 10.5]);  set(ax2,'YTick',[9.7 9.9 10.1 10.3]);   % bias:ĝ_I 9.54~10.41，4 內縮 tick
    end
    ylabel(ax2, '$\mathbf{{}^{B}\hat{g}_{I}\;(mT/A)}$', 'Interpreter','latex', 'FontSize',36);   % 下 panel y 軸標題（標準數學字體 \mathbf CM）
    hold(ax2,'off');

    bstr = ''; if USE_BIAS, bstr = '_bias'; end
    out = fullfile(figdir, sprintf('ell_gain_vs_R%s.png', bstr));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('\nwrote %s\n', out);
end

% ============================================================================
function style_panel(ax, FS)
% 選項①粗體框 @ paper scale(tick 由呼叫端明設,不減半)。
    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.02 .02]);
    ax.Toolbar.Visible = 'off';
end
