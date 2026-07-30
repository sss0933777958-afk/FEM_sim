function plot_ell_gV_vs_R(USE_BIAS)
% plot_ell_gV_vs_R -- long2016 半切六極 paper 圖:ℓ̂ 與 ᴮĝ_V 隨取樣半徑 R 的趨勢
% =========================================================================
%   voltage_base(no bias = fix 單一 ℓ,USE_BIAS=false)。風格 = plot_ell_gain_vs_R
%   的 voltage 版:上 panel ℓ̂(R)[µm] 與 current 完全相同(fitting 不分電流/電壓),
%   下 panel 換成 ᴮĝ_V(R)[mT/mV]。V(sensor 電壓 6×6)只算一次、與 R 無關;
%   loop R 只重跑 select_ball → fitting(fix) → solve_voltage。
%   輸出 → figures/paper_fig/Section3_A/ell_gV_vs_R.png(覆蓋迭代)。
% =========================================================================
    clc;
    if nargin < 1, USE_BIAS = false; end   % false = fix(no bias)
    here   = fileparts(mfilename('fullpath'));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CAL,'function'));  addpath(fullfile(CAL,'common_path'));

    % ---- 前段 pipeline(只做一次):graded / actuator frame ----
    cfg = model_config('long2016_hexapole_halfcut','tip40um');
    VARIANT  = cfg.default_variant;
    V_METHOD = fieldget(cfg,'v_method','csv-tet');
    raw = extract_ansys_data(cfg, 'all', VARIANT);      % 印 Matched N(指紋核對)
    ad  = build_actuator_data(raw, cfg);   Pc_base = ad.Pc_base;

    % ---- V:sensor 電壓 6×6 [mV]（只算一次，與 R 無關）----
    SOFF_upper = 4.572e-3;  n_uniform = 1e4;  sensor_r = 0.15e-3;  axial_tol = 0.10e-3;
    [V, ~] = build_V_matrix(cfg, VARIANT, raw, cfg.S_hall, SOFF_upper, n_uniform, sensor_r, axial_tol, [], V_METHOD);

    % ---- 掃描 R(fix 模型;每 20µm 一點;R=40~500µm)----
    Rum = 40:20:500;   nR = numel(Rum);
    ell_R = nan(1,nR);  gV_R = nan(1,nR);  npts_R = zeros(1,nR);
    fprintf('\n  R[um]   npts    ell_hat[um]   gV_hat[mT/mV]\n');
    for i = 1:nR
        [P, Bstack, npts] = cfg.select_ball(ad, Rum(i)*1e-6);
        npts_R(i) = npts;
        if npts < 3, fprintf('  %5d  %6d   (skip: npts<3)\n', Rum(i), npts); continue; end
        [e, l_hat, ~] = fitting(P, Bstack, Pc_base, 0.5e-3, USE_BIAS);
        [~, gV_hat]   = solve_voltage(l_hat, e, Pc_base, P, Bstack, V);
        ell_R(i) = l_hat*1e6;  gV_R(i) = gV_hat;
        fprintf('  %5d  %6d   %10.2f   %11.5f\n', Rum(i), npts, ell_R(i), gV_R(i));
    end

    % ---- 畫圖(2×1 panel,選項①粗體框;風格同 plot_ell_gain_vs_R)----
    cL = [0.05 0.10 0.95];   cR = [0.85 0.10 0.10];     % ℓ̂ 亮藍 / ĝ_V 紅
    FS = 30;
    fig = figure('Color','w','Position',[100 80 980 940]);
    t = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

    % 上:ℓ̂(R)（與 current 相同）
    ax1 = nexttile(t);  hold(ax1,'on');
    plot(ax1, Rum, ell_R, '-o', 'Color',cL, 'LineWidth',3.0, 'MarkerSize',8, 'MarkerFaceColor',cL);
    style_panel(ax1, FS);
    xlim(ax1,[40 500]);   set(ax1,'XTick',[40 100 200 300 400 500]);
    ylim(ax1,[750 870]);  set(ax1,'YTick',[760 810 860]);   % 縱軸 3 根等距(間距50,涵蓋峰值)
    ax1.XTickLabel = {};
    hold(ax1,'off');

    % 下:ĝ_V(R)  [mT/mV]（自動 nice ticks）
    ax2 = nexttile(t);  hold(ax2,'on');
    plot(ax2, Rum, gV_R, '-s', 'Color',cR, 'LineWidth',3.0, 'MarkerSize',8, 'MarkerFaceColor',cR);
    style_panel(ax2, FS);
    xlim(ax2,[40 500]);   set(ax2,'XTick',[40 100 200 300 400 500]);
    gmin = min(gV_R);  gmax = max(gV_R);  pad = 0.10*(gmax-gmin);
    ylim(ax2,[gmin-pad gmax+pad]);   set(ax2,'YTick', ticks3([gmin gmax]));   % 縱軸 3 根
    ylabel(ax2, '$\mathbf{{}^{B}\hat{g}_{V}\;(mT/mV)}$', 'Interpreter','latex', 'FontSize',36);
    hold(ax2,'off');

    bstr = ''; if USE_BIAS, bstr = '_bias'; end     % false→ell_gV_vs_R.png / true→..._bias.png
    out = fullfile(figdir, sprintf('ell_gV_vs_R%s.png', bstr));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('\nwrote %s\n', out);
    fprintf('gV range: %.5f .. %.5f mT/mV\n', gmin, gmax);
end

% ============================================================================
function v = fieldget(s, f, d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

% ---- 3 根等距 nice-step tick（內縮,含在範圍內）----
function tk = ticks3(lm)
    s0 = (lm(2)-lm(1))/3;
    p  = 10^floor(log10(s0));  c = [1 2 2.5 3 4 5 10]*p;
    [~,i] = min(abs(c-s0));  s = c(i);
    ctr = round((lm(1)+lm(2))/2/s)*s;
    tk  = ctr + [-1 0 1]*s;
end

% ============================================================================
function style_panel(ax, FS)
    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.02 .02]);
    ax.Toolbar.Visible = 'off';
end
