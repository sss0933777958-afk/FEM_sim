function plot_gain_iso_hist(USE_BIAS, R_um, force)
% plot_gain_iso_hist -- 兩個六極設計的控制指標分布疊圖（Long Fei vs Zhi-Peng）
% =========================================================================
%   對每個模型各自做一次校正（**R <= R_um 的 .fld 格點**，預設 150 µm），
%   再在**各自那批取樣點**上逐點算兩個控制指標，畫成疊圖直方圖：
%
%     圖 1  Actuation volume 的立方根  C^(1/3)  [mT/A]
%           C = ∏ σ_k（σ = svd(S(p)·ᴮĤ_I) 的奇異值）→ 開三次方後回到 mT/A，
%           與 ĝ_I 同量綱、可直接比較兩個設計的「等效致動強度」。
%     圖 2  Isotropy  κ = σ₃/σ₁  [無因次]
%
%   ⚠ 兩個模型的取樣球半徑相同（R_um），但**工作區尺度不同**：
%       Long Fei  R_norm = 500 µm；Zhi-Peng R_norm = 594 µm。
%     所以同一個 R 對兩者代表的「佔工作區比例」不同（30% vs 25%）。
%
%   ⚠ Zhi-Peng 的已知特性（首次載入實測，2026-08-15）：±x 兩根極的場是另外四根的
%     2.3 倍（極板外緣半徑 25.29 vs 34.76 mm → 磁路長度不同），K̄_I 有 8/30 個非對角
%     為正、NMAE 4.73%。點電荷模型對「六極不等強」配得不好，數值僅供比較用。
%
%   資料量：Zhi-Peng 的六個 .fld 各 ~1.15 GB（201³ 格點），讀一輪約 1.6 分 →
%   結果存快取 data/gain_iso_hist_<model>_R<R>_<tag>.mat，第二次起秒級。
%
%   風格：①粗體框圖 + 疊圖直方圖慣例（nb=180、共用 edges、mean 虛線、百分比縱軸、
%   圖例框外一行一系列、統計併入該行）。
%   輸出 → figures/paper_fig/Section4_C/{gain_cbrt,iso}_hist_maxwell.png
% =========================================================================
    clc;
    if nargin < 1 || isempty(USE_BIAS), USE_BIAS = false; end
    if nargin < 2 || isempty(R_um),     R_um     = 150;   end
    if nargin < 3 || isempty(force),    force    = false; end

    here   = fileparts(fileparts(mfilename('fullpath')));       % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section4_C');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    addpath(fullfile(CAL,'function'), fullfile(CAL,'utils'), fullfile(CAL,'common_path'));
    % sphere_grid_sample 目前放在 model 專屬子夾（已參數化成可跨 model，但檔案位置未動）
    addpath(fullfile(CAL,'utils','long2016_hexapole_halfcut'));

    tag = 'single';  if USE_BIAS, tag = 'eighteen'; end
    % {model, geom, 圖例名, 顏色}
    MD = { 'long2016_hexapole_halfcut', 'tip40um', 'Long Fei',  [0.05 0.10 0.95];
           'zhi_peng',                  '',        'Zhi-Peng',  [0.85 0.10 0.10] };

    nM = size(MD,1);   D = cell(1,nM);
    for a = 1:nM
        cf = fullfile(here, 'data', sprintf('gain_iso_hist_%s_R%03d_%s.mat', MD{a,1}, R_um, tag));
        if exist(cf,'file') && ~force
            D{a} = load(cf);   fprintf('由快取載入 %s\n', cf);
        else
            D{a} = compute_one(CAL, MD{a,1}, MD{a,2}, R_um, USE_BIAS, here);
            S = D{a};   save(cf, '-struct', 'S');
            fprintf('已存 %s\n', cf);
        end
        S = D{a};
        fprintf('[%-10s] N_c=%4d (%d,%d,%d)  評估 %4d 點  l_hat=%7.1f um  g_I=%7.4f mT/A  NMAE=%5.2f%%\n', ...
                MD{a,3}, S.Nc, S.tri, S.npts, S.l_hat*1e6, S.gI, S.NMAE);
        fprintf('             C^(1/3): %6.3f ~ %6.3f mT/A (mean %6.3f)   kappa: %.3f ~ %.3f (mean %.3f)\n', ...
                min(S.Ccbrt), max(S.Ccbrt), mean(S.Ccbrt), min(S.kap), max(S.kap), mean(S.kap));
    end

    % ---- 兩張疊圖 ----
    render_overlay({D{1}.Ccbrt, D{2}.Ccbrt}, MD(:,3), MD(:,4), ...
        '$\mathbf{\mathcal{C}^{1/3}\;(mT/A)}$', 'mT/A', '%.3f', ...
        fullfile(figdir, sprintf('gain_cbrt_hist_maxwell_%s.png', tag)));
    render_overlay({D{1}.kap,   D{2}.kap},   MD(:,3), MD(:,4), ...
        '$\mathbf{\kappa}$', '', '%.3f', ...
        fullfile(figdir, sprintf('iso_hist_maxwell_%s.png', tag)));
end

% ============================================================================
function S = compute_one(CAL, model, geom, R_um, USE_BIAS, here)
% 校正用 **R<=R_um 的收斂點設計 N_c**（降取樣，等測度網格）；
% 評估（畫直方圖）用 **R<=R_um 的全部真實 .fld 格點**。
    fprintf('--- %s：載入 + 校正（N_c 降取樣）---\n', model);
    cfg = model_config(model, geom);
    raw = extract_maxwell_data(cfg, 'all', cfg.default_variant);
    ad  = build_actuator_data(raw, cfg);
    [P, ~, np] = cfg.select_ball(ad, R_um*1e-6);      % 評估集（真實格點）

    F = zeros(6, cfg.N_I);
    for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end

    % ---- 收斂點設計 N_c（per model）----
    %   zhi_peng 的六極不等強 → K̄_I 非對角恆有正值，故關掉那道閘（ki_gate=false），
    %   只保留「ℓ̂ 穩定 ∧ ĝ_I 穩定 ∧ 對角全正 ∧ 對角占優」。
    %   與 plot_svd_polar 同一把尺（使用者拍板 2026-08-15）：非 long2016 時 K̄_I 的
    %   物理結構條件不適用（六極不等強）→ **完全不納入判準**（ki_req=false），
    %   收斂點只由 ℓ̂ 與 ĝ_I 決定。
    is_l2016 = strcmp(model,'long2016_hexapole_halfcut');
    sg = struct('model',model, 'geom',geom, 'ki_gate',is_l2016, 'ki_req',is_l2016);
    [tri, Nc] = conv_design(R_um*1e-6, cfg, F, USE_BIAS, here, cfg.R_norm, sg);
    assert(all(tri > 0), '%s：R=%d um 找不到收斂點設計', model, R_um);

    og = struct('frame','actuator', 'NRPT',tri, 'model',model, 'geom',geom);
    [gx,gy,gz,gB] = sphere_grid_sample(R_um*1e-6, [], og);
    Pc_ = [gx gy gz];
    Bc  = zeros(3*size(Pc_,1), size(gB,3));
    for j = 1:size(gB,3), Bc(:,j) = reshape(gB(:,:,j).', [], 1); end

    [e, l_hat] = fitting(Pc_, Bc, cfg.Pc_base, cfg.R_norm, USE_BIAS);
    [KI, gI, ~, rm] = solve_current(l_hat, e, cfg.Pc_base, Pc_, Bc, F);
    fprintf('  N_c = %d 點 (%d,%d,%d)；評估集 = %d 個真實格點\n', Nc, tri, np);

    Pc   = make_Pc(e, cfg.Pc_base);
    Hhat = gI * KI;                          % ᴮĤ_I [mT/A]
    C    = zeros(np,1);   kap = zeros(np,1);
    for i = 1:np
        d  = P(i,:)/l_hat - Pc.';            % 6×3
        Sm = (d ./ (vecnorm(d,2,2).^3)).';   % 3×6
        sv = svd(Sm * Hhat);
        C(i) = prod(sv);   kap(i) = sv(3)/sv(1);
    end
    S = struct('model',model, 'geom',geom, 'R_um',R_um, 'USE_BIAS',USE_BIAS, ...
               'npts',np, 'Nc',Nc, 'tri',tri, ...
               'l_hat',l_hat, 'gI',gI, 'KI',KI, 'NMAE',rm.NMAE, ...
               'C',C, 'Ccbrt',C.^(1/3), 'kap',kap);
end

% ============================================================================
function render_overlay(vals, names, cols, xlab, unit, fmt, out)
% 兩組資料的疊圖直方圖（各自正規化成百分比 → 比較的是分布形狀）
    % 配色沿用 err_hist_conv/shell 那組（使用者指定）：FaceAlpha 0.60 + **細黑邊**
    FS = 28;   ALPH = 0.60;   nb = 180;
    allv = [vals{1}(:); vals{2}(:)];
    edg  = linspace(min(allv), max(allv), nb+1);
    ctr  = (edg(1:end-1) + edg(2:end))/2;

    fig = figure('Color','w','Position',[100 100 1180 860]);
    ax  = axes(fig);   hold(ax,'on');
    h = gobjects(1,2);   pk = [];
    for a = 1:2
        p = histcounts(vals{a}, edg) / numel(vals{a}) * 100;
        h(a) = bar(ax, ctr, p, 1, 'FaceColor',cols{a}, 'FaceAlpha',ALPH, ...
                   'EdgeColor','k', 'LineWidth',0.3);
        pk = [pk p]; %#ok<AGROW>
    end
    for a = 1:2      % mean 虛線畫在長條之上（中性色，不隨系列色走）
        xline(ax, mean(vals{a}), '--', 'Color',[0 0 0]*(a==1) + [0 0.60 0]*(a==2), ...
              'LineWidth',2.6, 'HandleVisibility','off');
    end

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.015 .015]);
    [xr, xt] = xlim_pick(min(allv), max(allv));
    xlim(ax, xr);   set(ax,'XTick',xt);
    % 縱軸：自 0 起、上緣只留 8% 裕度、刻度取整數（figure-style）
    [yr, yt] = ylim_from_zero(max(pk));
    ylim(ax, yr);   set(ax,'YTick', yt);   ytop = yr(2);

    % x 起訖：只標數字、不畫 tick mark（figure-style）
    for xv = xr
        text(ax, xv, -0.022*ytop, sprintf('%g', xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    xlabel(ax, xlab, 'Interpreter','latex', 'FontSize',36);
    ylabel(ax, '$\mathbf{Percentage\;(\%)}$', 'Interpreter','latex', 'FontSize',36);

    % 圖例：一行一系列、標記方框、統計併入該行（2026-08-15 慣例）
    us = '';  if ~isempty(unit), us = [' ' unit]; end
    lb = cell(1,2);
    for a = 1:2
        lb{a} = sprintf(['%s: mean = ' fmt '%s, CV = %.2f%%'], names{a}, ...
                        mean(vals{a}), us, std(vals{a})/mean(vals{a})*100);
    end
    lg = legend(ax, h, lb, 'Interpreter','tex', 'Location','northoutside', 'NumColumns',1);
    lg.FontSize = 24;  lg.FontWeight = 'bold';
    lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;
    ax.Toolbar.Visible = 'off';   hold(ax,'off');

    drawnow;
    axp = get(ax,'Position');   lgp = get(lg,'Position');   lgh = lgp(4);
    GAPN = 0.020;   newTop = 1 - lgh - GAPN - 0.006;
    axp(4) = newTop - axp(2);   set(ax,'Position',axp);
    set(lg, 'Position', [axp(1), newTop + GAPN, axp(3), lgh]);

    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
    close(fig);
end

% ============================================================================
function [xr, xt] = xlim_pick(lo, hi)
% 橫軸：貼著資料選 nice 步長（**不強迫從 0 起**），取 3~7 個內部刻度；
%   端點由呼叫端以 text 補數字。cand 不含 2.5 → 盡量落在整數刻度（figure-style）。
    cand = [1 2 5 10];
    rng_ = max(hi-lo, realmin);
    for k = (floor(log10(rng_))-1) : (floor(log10(rng_))+1)
        for c = cand
            s  = c*10^k;
            x0 = floor(lo/s)*s;   x1 = ceil(hi/s)*s;
            t  = (round(x0/s)+1 : round(x1/s)-1) * s;
            if numel(t) >= 3 && numel(t) <= 7
                xr = [x0 x1];   xt = t;   return;
            end
        end
    end
    s  = rng_/5;   xr = [lo-0.1*rng_, hi+0.1*rng_];   xt = lo + (1:4)*s;
end

% ============================================================================
function [lim, tk] = ylim_from_zero(maxv)
% 自 0 起的縱軸：上緣只留 8% 裕度；刻度間距由資料範圍取 nice（不隨上緣壓縮而變小）
    cand = [1 2 2.5 3 4 5 10];
    x = maxv/4;   k = floor(log10(max(x,realmin)));
    s = cand(find(cand*10^k >= x, 1)) * 10^k;
    top = 1.08 * maxv;
    n = floor(top/s);   if n < 1, n = 1; end
    lim = [0 top];   tk = (1:n)*s;
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
