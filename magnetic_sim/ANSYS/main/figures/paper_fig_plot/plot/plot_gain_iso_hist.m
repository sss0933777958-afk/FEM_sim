function plot_gain_iso_hist(USE_BIAS, R_FIT, R_EVAL, force, BINW)
% plot_gain_iso_hist -- 兩個六極設計的控制指標分布疊圖（Long Fei vs Zhi-Peng）
% =========================================================================
%   [MODIFIED 2026-08-21 使用者要求] **校正半徑與評估半徑分離**：
%     校正 = R <= R_FIT（預設 150 µm）的 **N_c 收斂點設計**（等測度網格降取樣）
%     評估 = R <= R_EVAL（預設 500 µm）的 **全部真實 .fld 格點**
%   （舊版兩者都用同一個 R_um。）兩個模型用同一組 (R_FIT, R_EVAL) → 可直接互比。
%
%     圖 1  Actuation volume 的立方根  C^(1/3)  [mT/A]
%           C = ∏ σ_k（σ = svd(S(p)·ᴮĤ_I) 的奇異值）→ 開三次方後回到 mT/A，
%           與 ĝ_I 同量綱、可直接比較兩個設計的「等效致動強度」。
%     圖 2  Isotropy  κ = σ₃/σ₁  [無因次]
%
%   ⚠ C 與 κ 是**點電荷模型導出的解析函數**，不是 raw FEM 場 —— 這裡的「節點」只用來
%     決定「在哪些位置求值」，並未讀取該點的 FEM 場值。
%   ⚠ **R_EVAL > R_FIT 時外圈是外推**：模型只在 R <= R_FIT 內被資料約束過。
%
%   ⚠ 兩個模型現在**工作區尺度相同**：Long Fei R_norm = 500 µm、
%     Zhi-Peng **R500** 變體 R_norm = 500 µm → R_EVAL = 500 對兩者都正好是工作球。
%     （舊版用 zhi_peng R594，兩者尺度不同、需另行換算，已不適用。）
%
%   ⚠ Zhi-Peng 的 'maxwell_split' 場（2026-08-20 改導磁係數後重解）六極已趨於等強，
%     K̄_I 非對角 30/30 全負、對角極差 12% —— 舊版檔頭寫的「±x 兩根強 2.3 倍、
%     8/30 非對角為正」是**舊場**（'maxwell'）的特性，已不適用。
%
%   資料量：兩個模型的六個 .fld 各 ~1.15 GB（201³ 格點），讀一輪約 1.6 分 →
%   結果存快取 data/gain_iso_hist_<model>_<variant>_fit<R>_eval<R>_<tag>.mat，第二次起秒級。
%
%   風格：①粗體框圖 + 疊圖直方圖慣例（nb=180、**共用 edges → 兩組 bin 寬相同**、
%   百分比縱軸）。[MODIFIED 2026-08-21] 照 err_hist 家族的定案：**不畫 mean 虛線、
%   長條不描黑邊、圖例只列系列名**；mean / CV / min / max 一律印在 console。
%   輸出 → figures/paper_fig/Section4_C/{gain_cbrt,iso}_hist_maxwell_<tag>.png
% =========================================================================
    clc;
    if nargin < 1 || isempty(USE_BIAS), USE_BIAS = false; end
    if nargin < 2 || isempty(R_FIT),    R_FIT    = 150;   end   % 校正半徑 [µm]
    if nargin < 3 || isempty(R_EVAL),   R_EVAL   = 500;   end   % 評估半徑 [µm]
    if nargin < 4 || isempty(force),    force    = false; end
    % [ADDED 2026-08-21] BINW = 直接指定 C^(1/3) 圖的 bin 寬 [mT/A]；[] = 沿用 nb=180。
    %   起因：兩個設計相距 23.1 mT/A（是 Design A 自身寬度的 39 倍），共用 edges 之下
    %   180 根裡有 157 根落在中間的空白 -> bin 寬被空白決定，A 的整個分布只佔 4 根、
    %   看起來像一根尖刺。給 BINW 可把格子切細（例 0.05 -> A 約 12 根）。
    %   ⚠ 只套用在 C^(1/3) 那張；κ 的值域是 [0,1]、尺度完全不同，仍用 nb=180。
    if nargin < 5, BINW = []; end
    % [MODIFIED 2026-08-21] BINW 可給 1x2 = [gain_bin, kappa_bin]，分別套到兩張圖；
    %   給純量 = 只套 C^(1/3)（舊行為）；[] = 兩張都用 nb=180。
    %   預設值由 Freedman-Diaconis 決定（h = 2*IQR*n^(-1/3)）：四組資料一致指向
    %   「每個分布切 16~20 根」。C 的 A/B 寬度差 4.6 倍、FD 各要 0.033/0.165，
    %   共用 edges 只能取一個 -> 折衷 0.05（A 12 根、B 雜訊 17%%）。
    %   kappa 兩組寬度幾乎相同、FD 給 0.0106/0.0109 -> 取 0.011（各 19 根、雜訊 10%%）。
    %   ⚠ 舊值 nb=180 對 kappa 是**過度解析 3.5 倍**（每根僅 26 點、雜訊 20%%），
    %     圖上的鋸齒是統計噪音而非結構。
    if isempty(BINW),      BW_G = 0.05;   BW_K = 0.011;
    elseif isscalar(BINW), BW_G = BINW;   BW_K = 0.011;
    else,                  BW_G = BINW(1); BW_K = BINW(2);
    end

    here   = fileparts(fileparts(mfilename('fullpath')));       % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section4_C');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    addpath(fullfile(CAL,'function'), fullfile(CAL,'utils'), fullfile(CAL,'common_path'));
    % [MODIFIED 2026-08-21] 移除 addpath(utils/long2016_hexapole_halfcut)：utils/ 已扁平化，
    %   conv_design_ws 在 function/（上一行已涵蓋），舊路徑不存在。

    tag = 'single';  if USE_BIAS, tag = 'eighteen'; end
    % {model, geom, variant, 圖例名, 顏色}
    %   [MODIFIED 2026-08-21] ① zhi_peng 改用 **R500** 幾何 + **maxwell_split** 場
    %   （2026-08-20 改導磁係數後重解的版本）；② 加 variant 欄位 —— 不明給的話
    %   conv_design_ws 會退回 cfg.default_variant，**靜默**用舊場。
    %   顏色照「顏色 = 模型」慣例：Long Fei 深藍、Zhi-Peng 紅。
    %   [MODIFIED 2026-08-21 使用者拍板] 圖例名改成 **Design A / Design B**（原 Long Fei /
    %   Zhi-Peng）；console 仍印 model 名，追溯不受影響。
    MD = { 'long2016_hexapole_halfcut', 'tip40um', '',              'Design A',  [0.05 0.10 0.95];
           'zhi_peng',                  'R500',    'maxwell_split', 'Design B',  [0.85 0.10 0.10] };

    nM = size(MD,1);   D = cell(1,nM);
    for a = 1:nM
        vt = MD{a,3};  if isempty(vt), vt = 'default'; end
        cf = fullfile(here, 'data', sprintf('gain_iso_hist_%s_%s_fit%03d_eval%03d_%s.mat', ...
                                            MD{a,1}, vt, R_FIT, R_EVAL, tag));
        if exist(cf,'file') && ~force
            D{a} = load(cf);   fprintf('由快取載入 %s\n', cf);
        else
            D{a} = compute_one(MD{a,1}, MD{a,2}, MD{a,3}, R_FIT, R_EVAL, USE_BIAS, here, CAL);
            S = D{a};   save(cf, '-struct', 'S');
            fprintf('已存 %s\n', cf);
        end
    end

    % ---- console 統計（使用者要求：mean 印在 console，不進圖）----
    fprintf('%s\n', repmat('=',1,78));
    fprintf('校正 R <= %d um（N_c 降取樣）｜評估 R <= %d um（真實 .fld 格點）｜%s\n', ...
            R_FIT, R_EVAL, tag);
    fprintf('%s\n', repmat('-',1,78));
    for a = 1:nM
        S = D{a};
        fprintf('[%-9s] %s / %s：N_c=%d (%d,%d,%d)  評估 %d 點  l_hat=%.1f um  g_I=%.4f mT/A  NMAE=%.2f%%\n', ...
                MD{a,4}, S.model, S.variant, S.Nc, S.tri, S.npts, S.l_hat*1e6, S.gI, S.NMAE);
        fprintf('    C^(1/3)  mean = %8.4f mT/A   CV = %5.2f%%   min = %8.4f   max = %8.4f\n', ...
                mean(S.Ccbrt), std(S.Ccbrt)/mean(S.Ccbrt)*100, min(S.Ccbrt), max(S.Ccbrt));
        fprintf('    kappa    mean = %8.4f         CV = %5.2f%%   min = %8.4f   max = %8.4f\n', ...
                mean(S.kap), std(S.kap)/mean(S.kap)*100, min(S.kap), max(S.kap));
    end
    fprintf('%s\n', repmat('-',1,78));
    fprintf('mean 比值（Zhi-Peng / Long Fei）：C^(1/3) %.3f 倍｜kappa %.3f 倍\n', ...
            mean(D{2}.Ccbrt)/mean(D{1}.Ccbrt), mean(D{2}.kap)/mean(D{1}.kap));
    fprintf('%s\n', repmat('=',1,78));

    % ---- 兩張疊圖 ----
    % [ADDED 2026-08-21 使用者拍板] 水平軸範圍可明給（[] = 由 xlim_pick 自動選）。
    %   C^(1/3) 定成 **[10, 50]**：自動選出來的 [5,55] 兩端各留 s/2 = 5，右邊那 4.1 個
    %   單位是空的、看起來「右邊留太多」。改成 10~50 後端點就是刻度、完全不留白。
    %   ⚠ 代價：Design B 有一小段尾巴 > 50 會落在視野外（比例印在 console）。
    %   ⚠ [10,50] 是**針對 R_EVAL=500 的資料**挑的；換評估半徑資料範圍就不同（R<=150
    %     的 C^(1/3) 只落在 12.2~16.6 / 26.9~35.6），沿用會空掉一大半 → 其餘半徑走自動。
    XR_GAIN = [];   if R_EVAL == 500, XR_GAIN = [10 50]; end
    XR_ISO  = [];                       % κ 的自動結果 [0,1] 已是最緊、不必 override
    % [MODIFIED 2026-08-21] 檔名帶評估半徑 _R<eval>，讓不同 R_EVAL 的圖並存
    %   （原本沒帶，R150 版會直接蓋掉 R500 版）。
    render_overlay({D{1}.Ccbrt, D{2}.Ccbrt}, MD(:,4), MD(:,5), ...
        '$\mathbf{\mathcal{C}^{1/3}\;(mT/A)}$', ...
        fullfile(figdir, sprintf('gain_cbrt_hist_maxwell_%s_R%d.png', tag, R_EVAL)), XR_GAIN, BW_G);
    render_overlay({D{1}.kap,   D{2}.kap},   MD(:,4), MD(:,5), ...
        '$\mathbf{\kappa}$', ...
        fullfile(figdir, sprintf('iso_hist_maxwell_%s_R%d.png', tag, R_EVAL)), XR_ISO, BW_K);
end

% ============================================================================
function S = compute_one(model, geom, variant, R_FIT, R_EVAL, USE_BIAS, here, CAL)
% 校正用 **R<=R_FIT 的收斂點設計 N_c**（降取樣，等測度網格）；
% 評估（畫直方圖）用 **R<=R_EVAL 的全部真實 .fld 格點**。
    fprintf('--- %s：載入 + 校正（N_c 降取樣 @R<=%d um）---\n', model, R_FIT);
    cfg = model_config(model, geom);
    if isempty(variant), variant = cfg.default_variant; end
    raw = extract_maxwell_data(cfg, 'all', variant);
    ad  = build_actuator_data(raw, cfg);
    [P, ~, np] = cfg.select_ball(ad, R_EVAL*1e-6);    % 評估集（真實格點）

    F = zeros(6, cfg.N_I);
    for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end

    % ---- 收斂點設計 N_c（per model）----
    %   zhi_peng 的六極不等強 → K̄_I 非對角恆有正值，故關掉那道閘（ki_gate=false），
    %   只保留「ℓ̂ 穩定 ∧ ĝ_I 穩定 ∧ 對角全正 ∧ 對角占優」。
    %   與 plot_svd_polar 同一把尺（使用者拍板 2026-08-15）：非 long2016 時 K̄_I 的
    %   物理結構條件不適用（六極不等強）→ **完全不納入判準**（ki_req=false），
    %   收斂點只由 ℓ̂ 與 ĝ_I 決定。
    is_l2016 = strcmp(model,'long2016_hexapole_halfcut');
    sg = struct('model',model, 'geom',geom, 'variant',variant, ...
                'ki_gate',is_l2016, 'ki_req',is_l2016);
    % ---- 讀 main.m 產的收斂設計校正結果（不再自己重跑階梯 + 校正）--------
    %   [MODIFIED 2026-08-23 使用者拍板] 校正與收斂判準已搬回 main.m
    %   （conv_design_ws / conv_design_sensor 只負責決定內插點位置與取場），
    %   繪圖端改成**接收 main 產完的結果** -> 圖與結果 PDF 保證出自同一次校正。
    %   ⚠ 該組合必須先跑過 main.m（GRID_NRPT='auto'）；找不到就報錯，不猜。
    %   ⚠ 同一組合可能有多顆 convN 檔（舊實驗留下的，例如 long2016 R150
    %     eighteen 就有 convN6/convN80/convN88）-> 只認 conv_auto==true 那顆。
    md_ = fullfile(CAL, 'data', model, '.mat');
    tg_ = 'single';   if USE_BIAS, tg_ = 'eighteen'; end
    dd_ = dir(fullfile(md_, sprintf('calib_current_%s_convN*_R%03d_%s.mat', ...
                                    variant, round(R_FIT), tg_)));
    if numel(dd_) > 1
        ok_ = false(1, numel(dd_));
        for k_ = 1:numel(dd_)
            f_ = fullfile(md_, dd_(k_).name);   w_ = whos('-file', f_);
            if ismember('conv_auto', {w_.name})
                r_ = load(f_, 'conv_auto');   ok_(k_) = logical(r_.conv_auto);
            end
        end
        dd_ = dd_(ok_);
    end
    assert(numel(dd_) == 1, ['找到 %d 顆收斂校正 .mat（需恰好 1 顆）。請先跑 ' ...
           'main.m：MODEL=''%s''、R_select=%ge-6、USE_BIAS=%d、GRID_NRPT=''auto''。'], ...
           numel(dd_), model, R_FIT, USE_BIAS);
    cal = load(fullfile(md_, dd_(1).name));
    e = cal.e;   l_hat = cal.l_hat;   KI = cal.KI_bar;   gI = cal.gI_hat;
    tri = cal.GRID_NRPT;   Nc = cal.npts;   rm = struct('NMAE', cal.NMAE);
    fprintf('  N_c = %d 點 (%d,%d,%d) @R<=%d um；評估集 = %d 個真實格點 @R<=%d um\n', ...
            Nc, tri, R_FIT, np, R_EVAL);

    Pc   = make_Pc(e, cfg.Pc_base);
    Hhat = gI * KI;                          % ᴮĤ_I [mT/A]
    C    = zeros(np,1);   kap = zeros(np,1);
    for i = 1:np
        d  = P(i,:)/l_hat - Pc.';            % 6×3
        Sm = (d ./ (vecnorm(d,2,2).^3)).';   % 3×6
        sv = svd(Sm * Hhat);
        C(i) = prod(sv);   kap(i) = sv(3)/sv(1);
    end
    S = struct('model',model, 'geom',geom, 'variant',variant, ...
               'R_FIT',R_FIT, 'R_EVAL',R_EVAL, 'USE_BIAS',USE_BIAS, ...
               'npts',np, 'Nc',Nc, 'tri',tri, ...
               'l_hat',l_hat, 'gI',gI, 'KI',KI, 'NMAE',rm.NMAE, ...
               'C',C, 'Ccbrt',C.^(1/3), 'kap',kap);
end

% ============================================================================
function render_overlay(vals, names, cols, xlab, out, xrset, BINW)
% 兩組資料的疊圖直方圖（各自正規化成百分比 → 比較的是分布形狀）
% [MODIFIED 2026-08-21 使用者拍板，比照 err_hist 家族]：
%   ① 長條**不描黑邊**（EdgeColor 'none'）——細長條逐根描邊會糊成一團黑
%   ② **不畫 mean 虛線**；③ 圖例**只列系列名**（mean / CV 改印 console）
%   ④ 縱軸換成 figure-style 2026-08-20 定案：內部刻度標數字、起點與終點都不標、
%      [0,T] 平分 N+1 段、步長為 0.1 的倍數且數字為整數或 0.N
%   ⑤ 兩組**共用同一組 edges → bin 寬必然相同**（figure-style 硬條件：
%      面積 = 100% × bin 寬，寬度不同就不能比高度）
    if nargin < 6, xrset = []; end
    if nargin < 7, BINW  = []; end
    FS = 28;   ALPH = 0.60;   nb = 180;
    allv = [vals{1}(:); vals{2}(:)];
    if isempty(BINW)
        edg = linspace(min(allv), max(allv), nb+1);
    else
        % [ADDED 2026-08-21] 明給 bin 寬：自 min 起等寬鋪到蓋過 max。仍是兩組共用
        %   同一組 edges -> 面積 = 100%% x BINW 的可比性不變。
        edg = min(allv) : BINW : (min(allv) + ceil((max(allv)-min(allv))/BINW)*BINW);
        nb  = numel(edg) - 1;
    end
    % **兩組共用同一組 edges** → 兩組 bin 寬必然相同（面積可比的硬條件）
    ctr  = (edg(1:end-1) + edg(2:end))/2;
    fprintf('  bin：共用 %d 個、寬 %.4g（兩組同寬；由兩組合併後的 %.4g ~ %.4g 均分）\n', ...
            nb, edg(2)-edg(1), min(allv), max(allv));
    for q = 1:2                                   % [ADDED] 各組實際佔幾根 / 峰值多高
        wv = vals{q}(:);   hc = histcounts(wv, edg);
        fprintf(['    %-10s 全寬 %.4g -> %.1f 根（非空 %d），峰值 %.1f%%' char(10)], ...
                names{q}, max(wv)-min(wv), (max(wv)-min(wv))/(edg(2)-edg(1)), ...
                sum(hc>0), max(hc)/numel(wv)*100);
    end

    fig = figure('Color','w','Position',[100 100 1180 860]);
    ax  = axes(fig);   hold(ax,'on');
    h = gobjects(1,2);   pk = [];
    for a = 1:2
        p = histcounts(vals{a}, edg) / numel(vals{a}) * 100;
        h(a) = bar(ax, ctr, p, 1, 'FaceColor',cols{a}, 'FaceAlpha',ALPH, 'EdgeColor','none');
        pk = [pk p]; %#ok<AGROW>
    end

    box(ax,'on');  grid(ax,'off');
    % [MODIFIED 2026-08-21 使用者拍板] tick 朝外（TickDir 'out'，同 plot_sigma_hist）
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5, ...
           'TickLength',[.015 .015],'TickDir','out');
    if isempty(xrset)
        [xr, xt] = xlim_pick(min(allv), max(allv));
    else
        xr = xrset;   xt = xticks_in(xr);            % 使用者指定範圍 → 只挑內部刻度
        for a = 1:2                                  % 落在視野外的樣本比例（誠實回報）
            f = 100*mean(vals{a}(:) < xr(1) | vals{a}(:) > xr(2));
            if f > 0, fprintf('  ⚠ %s 有 %.2f%% 的點落在視野 [%g, %g] 之外\n', ...
                              names{a}, f, xr(1), xr(2)); end
        end
    end
    xlim(ax, xr);   set(ax,'XTick',xt);
    [yr, yt] = ylim_from_zero(max(pk), 4);
    ylim(ax, yr);   set(ax,'YTick', yt);   ytop = yr(2);

    % x 起訖：只標數字、不畫 tick mark（figure-style）
    for xv = xr
        text(ax, xv, -0.022*ytop, sprintf('%g', xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    xlabel(ax, xlab, 'Interpreter','latex', 'FontSize',36);
    ylabel(ax, '$\mathbf{Percentage\;(\%)}$', 'Interpreter','latex', 'FontSize',36);

    % 圖例：框內右上角、一行一系列、只列系列名
    lg = legend(ax, h, names, 'Interpreter','tex', 'Location','northeast', 'NumColumns',1);
    lg.FontSize = 24;  lg.FontWeight = 'bold';
    lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;  lg.Color = 'w';
    ax.Toolbar.Visible = 'off';   hold(ax,'off');

    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
    close(fig);
end

% ============================================================================
function [xr, xt] = xlim_pick(lo, hi)
% 橫軸（figure-style 2026-08-19/20/21）。硬條件，優先序由高到低：
%   ① 內部刻度**數量奇數**（3 或 5）  ② 刻度**等距**
%   ③ 每個刻度與兩個端點都是**乾淨的數**（整數或一位小數 0.N）
%   ④ 兩端留白盡量小（端點另以 text 標數字、不畫 tick）
% ⚠ 別走回頭路：舊版「先由資料定 x0=floor(lo/s)*s、x1=ceil(hi/s)*s，再取內部刻度」
%   常常無解 —— 它要求 (x1−x0)/s 剛好是 4 或 6，而 x0/x1 又被資料釘死。2026-08-21 實測
%   C^(1/3)（12.2~50.9）在那個寫法下完全無解，掉進 fallback 後端點變成 8.3426 / 54.7651、
%   刻度互相疊字。
% 現行做法 = 「**先定刻度、再由刻度反推 xlim**」：
%   刻度取 o + m*s（o ∈ {0, s/2} 兩種相位），xlim = [t_first − s/2, t_last + s/2]
%   → 兩端留白恆為 s/2（對稱、且是最小可能值）。挑「總留白最小」者，同分再挑
%   「單邊留白最大值最小」者（避免一邊貼很緊、另一邊空一大塊 —— 使用者反映的
%   「右邊留太多」正是這種：舊版 [10,60] 右邊空 9.1，新版 [5,55] 只空 4.1）。
%   ① 端點取 **s/2 或 s/4 的整數倍**（都保證乾淨），由資料 floor/ceil 而來 → 貼著資料
%   ② 刻度取 o + k*s（o ∈ {0, s/2} 兩種相位）落在 (x0,x1) **嚴格內部**者
%   ③ 端點各再試「往外推一格」一次 —— 只用 ①② 常常湊不出奇數根（實測 R500 的
%      12.2~50.9 在基礎版全是偶數 4 或 8 根），往外推可以翻轉奇偶
%   排序：**整數刻度優先** → 留白小 → 根數接近 5
% [MODIFIED 2026-08-21 使用者反映「兩邊留太寬」] 端點格加入 **s/4**：
%   只用 s/2 時，「整數刻度」會把端點鎖在 s/2 的倍數上 —— C^(1/3) R150（14.36~40.77）
%   因此只能取 [10,45]（兩邊各空 4.4 / 4.2）。加了 s/4 之後可取 [12.5, 42.5]，
%   刻度仍是整數 20/30/40，但留白縮到 1.86 / 1.73（少 58%）。
    rng_ = max(hi-lo, realmin);
    cand = [1 2 2.5 5 10];
    xr = [];   xt = [];   best = [inf inf inf];
    for k = (floor(log10(rng_))-2) : (floor(log10(rng_))+1)
        for c = cand
            s = c*10^k;
            for hq = [s/2, s/4]                            % 端點格（粗 / 細兩種都試）
                b0 = floor(lo/hq)*hq;    b1 = ceil(hi/hq)*hq;
                for o = [0, s/2]                           % 刻度相位
                    for ext = 1:3
                        x0 = b0;   x1 = b1;
                        if ext == 2, x0 = b0 - hq; elseif ext == 3, x1 = b1 + hq; end
                        t = (ceil((x0-o)/s + 1e-9) : floor((x1-o)/s - 1e-9))*s + o;
                        t = t(t > x0+1e-9 & t < x1-1e-9);
                        n = numel(t);
                        if mod(n,2)~=1 || n<3 || n>5, continue; end
                        if ~(isclean(x0) && isclean(x1) && isclean(t)), continue; end
                        isInt = double(any(abs(t - round(t)) > 1e-9));   % 0 = 全整數（優先）
                        sc = [isInt, (lo-x0)+(x1-hi), abs(n-5)];
                        if lexlt(sc, best), best = sc;   xr = [x0 x1];   xt = t; end
                    end
                end
            end
        end
    end
    if ~isempty(xr), return; end
    s  = rng_/6;   xr = [lo-0.1*rng_, hi+0.1*rng_];   xt = lo + (1:5)*s;   % 保底
end

% ============================================================================
function tf = lexlt(a, b)      % 字典序比較（容忍浮點）
    for i = 1:numel(a)
        if a(i) < b(i) - 1e-9, tf = true;  return; end
        if a(i) > b(i) + 1e-9, tf = false; return; end
    end
    tf = false;
end

% ============================================================================
function xt = xticks_in(xr)
% 給定的 [x0 x1] 內挑內部刻度：等距、**奇數個**（3 或 5）、每個都是乾淨的數。
%   兩種相位都試（s 的整數倍 / 再偏半格），取最接近 5 根的。端點本身不放 tick
%   —— 端點值由呼叫端以 text 標，照 figure-style「端點只標數字、不畫 tick mark」。
    span = xr(2) - xr(1);
    cand = [1 2 2.5 5 10];
    xt = [];   bestd = inf;
    for k = (floor(log10(span))-2) : (floor(log10(span))+1)
        for c = cand
            s = c*10^k;
            for o = [0, s/2]
                t = (ceil((xr(1)-o)/s + 1e-9) : floor((xr(2)-o)/s - 1e-9))*s + o;
                t = t(t > xr(1)+1e-9 & t < xr(2)-1e-9);
                n = numel(t);
                if mod(n,2)~=1 || n<3 || n>5 || ~isclean(t), continue; end
                if abs(n-5) < bestd, bestd = abs(n-5);  xt = t; end
            end
        end
    end
    if isempty(xt), xt = xr(1) + (1:3)*(span/4); end          % 保底
end

% ============================================================================
function tf = isclean(v)
% 「乾淨的數」＝整數或一位小數 0.N（figure-style 2026-08-20：端點不可是 0.123 這種尾數）
    tf = all(abs(v*10 - round(v*10)) < 1e-9);
end

% ============================================================================
function [lim, tk] = ylim_from_zero(maxv, N)
% 縱軸（figure-style 2026-08-20 定案，與 plot_err_hist_shell 同一支）：
%   • 自 0 起；**終點（軸上緣）與起點 0 都不標數字、也不畫 tick**（端點值由水平軸負責）
%   • **內部刻度要標數字**，每個必須是整數或一位小數 0.N
%   • 「N 根 tick」= 不含端點的內部刻度數；[0,T] 平分 N+1 段 → tk=(1:N)*s、T=(N+1)*s
%   • s 必為 0.1 的倍數且 (N+1)*s >= 1.08*maxv；在 [smin, 1.15*smin] 內挑最漂亮的
%     （整數 > 0.5 的倍數 > 0.2 的倍數 > 其餘；同分取最小 s，填充率最高）
    if nargin < 2 || isempty(N), N = 4; end
    smin = 1.08*maxv/(N+1);
    k0 = max(1, ceil(smin/0.1 - 1e-9));   k1 = max(k0, ceil(1.15*smin/0.1));
    best = k0;   bs = -1;
    for k = k0:k1
        if     mod(k,10) == 0, sc = 3;               % 整數
        elseif mod(k,5)  == 0, sc = 2;               % 0.5 的倍數
        elseif mod(k,2)  == 0, sc = 1;               % 0.2 的倍數
        else,                  sc = 0;
        end
        if sc > bs, bs = sc;  best = k; end
    end
    s   = best*0.1;
    tk  = round((1:N)*s*10)/10;
    lim = [0 (N+1)*s];
end

% ============================================================================
function Pc = make_Pc(e17, Pc_base)
    if isempty(e17) || all(e17(:) == 0), Pc = Pc_base;  return; end
    E = zeros(3, 6);
    E(:,1) = e17(1:3);     E(:,2) = e17(4:6);
    E(:,3) = e17(7:9);     E(:,4) = e17(10:12);
    E(:,5) = e17(13:15);
    E(1,6) = e17(16);      E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);
    Pc = Pc_base + E;
end
