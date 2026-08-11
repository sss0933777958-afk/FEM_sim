function plot_ell_nmin_vs_full(QTY, RSEL, TOLPC, KWIN, force)
% plot_ell_nmin_vs_full -- paper 圖：減量取樣 vs 全取樣的 l_hat(R)，並算兩者的 RMSPE
% =========================================================================
%   兩條曲線（同一組 10 um 內插格、同一條校正管線、single 模型）：
%     ① **減量取樣**：每個 R 只取 N_min(R) 個隨機點（10 組排列取中位數）
%     ② **全取樣**  ：每個 R 取球內**全部**格點
%   縱軸 l_hat [um]、水平軸取樣半徑 R [um]；另算 RMSPE = rms((①-②)/②)。
%
%   N_min(R) 的定義（使用者拍板 2026-08-12）：對每組排列，找第一個 n 使
%   **其後 KWIN 點**的相對變化率都 <= TOLPC%；10 組取中位數。
%   ⚠ K=10 是必要的：K=5 會在震盪期產生偽陽性（n=42/57/67 這種不可能的值）。
%
%   🔑 **本圖不做任何新擬合**，全部從既有快取取值：
%     ① data/nmin_vs_R_maxwell10_single_random.mat（10 組 × 每個 n 的 l_hat）
%     ② data/ell_gain_sweep_maxwell_interp.mat（每個 R 用完整 10 um 池的 l_hat）
%     兩者是同一套網格與管線，可直接逐點相比。
%
%   風格 = 選項①粗體框 @ paper scale；圖例照 figure-style「圖例標準樣式」。
%   輸出 → figures/paper_fig/Section2_E/ell_nmin_vs_full_maxwell10_single.png
% =========================================================================
    clc;
    % [MODIFIED 2026-08-12] 起點拉回 R=100、視窗改 K=7（使用者拍板）。
    %   K=10 時 R=100 判不出來 —— 原因是**技術性的**：n 網格 log 間距、固定 45 點，池越小
    %   跨度越短，安定區（n>~1500）之後只剩 8.4 個格點，湊不出要檢查的 10 點。
    %   K=7 讓 21/21 個 R 都有解，且 10 組中位數**沒有一個**被偽陽（<200）污染。
    % [ADDED 2026-08-12] QTY='ell'（預設，從快取取值、零成本）| 'gain'（ĝ_I，需在
    %   那 21 個 (R, N_min) 組合上重跑擬合 —— 快取沒存 ĝ_I；結果另存 nmin_gain_*.mat）。
    if nargin < 1 || isempty(QTY),   QTY   = 'ell';      end
    if nargin < 2 || isempty(RSEL),  RSEL  = 100:20:500; end
    if nargin < 3 || isempty(TOLPC), TOLPC = 0.2;        end
    if nargin < 4 || isempty(KWIN),  KWIN  = 7;          end
    if nargin < 5 || isempty(force), force = false;      end
    switch lower(QTY)
        case 'ell',  ylab = '$\mathbf{\hat{\ell}\;(micro\;meter)}$';   fmt = '%8.3f';
        case 'gain', ylab = '$\mathbf{{}^{B}\hat{g}_{I}\;(mT/A)}$';    fmt = '%8.4f';
        otherwise,   error('QTY 必為 ''ell'' | ''gain''');
    end

    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    S = load(fullfile(here,'data','nmin_vs_R_maxwell10_single_random.mat'));   % 減量（10 組）
    F = load(fullfile(here,'data','ell_gain_sweep_maxwell_interp.mat'));       % 全取樣

    nR = numel(RSEL);
    ellA = nan(1,nR);  ellB = nan(1,nR);  Nuse = nan(1,nR);  pool = nan(1,nR);
    fprintf('   R    N_used   pool     佔比%%    ell_減量   ell_全取樣    誤差%%\n');
    for t = 1:nR
        i = find(S.RLIST == RSEL(t));
        assert(~isempty(i), '快取沒有 R=%d', RSEL(t));
        n = S.nall{i};   E = S.ellall{i};   pool(t) = S.pool(i);

        nd = nan(1, size(E,1));                       % 各組的 N_min
        for d = 1:size(E,1)
            nd(d) = flat_one(n, E(d,:), TOLPC, KWIN);
        end
        Nmed = median(nd, 'omitnan');
        [~, k] = min(abs(n - Nmed));                  % 取最接近中位數的 n 格
        Nuse(t) = n(k);
        ellA(t) = median(E(:,k), 'omitnan');          % 減量：該 n 的 10 組中位數
        ellB(t) = F.ell_R(F.Rum == RSEL(t));          % 全取樣
    end

    if strcmpi(QTY,'gain')                             % ĝ_I 需重跑擬合（見檔頭 [ADDED]）
        [ellA, ellB] = gain_curves(here, RSEL, Nuse, S, F, force);
    end
    for t = 1:nR
        fprintf(['  %3d  %7d  %7d  %6.2f   ' fmt '   ' fmt '   %+7.3f\n'], ...
                RSEL(t), Nuse(t), pool(t), Nuse(t)/pool(t)*100, ellA(t), ellB(t), ...
                (ellA(t)-ellB(t))/ellB(t)*100);
    end
    rel   = (ellA - ellB) ./ ellB;
    RMSPE = sqrt(mean(rel.^2)) * 100;
    fprintf('\n  RMSPE = %.4f %%   最大絕對誤差 %.3f %%   點數佔比 %.2f ~ %.2f %%\n', ...
            RMSPE, max(abs(rel))*100, min(Nuse./pool)*100, max(Nuse./pool)*100);

    % ---- 畫圖（單 panel，選項①粗體框）----
    FS = 36;  LWBOX = 4.0;
    cA = [0.05 0.10 0.95];   cB = [0.85 0.10 0.10];
    fig = figure('Color','w','Position',[100 100 1160 830]);
    ax  = axes(fig);  hold(ax,'on');
    h1 = plot(ax, RSEL, ellA, '-o', 'Color',cA, 'LineWidth',3.0, 'MarkerSize',8, ...
              'MarkerFaceColor',cA, 'Clipping','off');
    h2 = plot(ax, RSEL, ellB, '--s','Color',cB, 'LineWidth',3.0, 'MarkerSize',8, ...
              'MarkerFaceColor',cB, 'Clipping','off');
    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX,'TickLength',[.015 .015]);
    ax.Toolbar.Visible = 'off';

    XL = [RSEL(1) RSEL(end)];   xlim(ax, XL);
    set(ax,'XTick', inner_ticks(XL(1), XL(2), 3));
    [YL, YT] = axlim_auto(min([ellA ellB]), max([ellA ellB]), [3 5]);
    ylim(ax, YL);   set(ax,'YTick',YT);

    % [REMOVED 2026-08-12] 原本每點標 N_min 的數字已拿掉（使用者要求）——
    %   點數佔比改由獨立的火柴棒圖 plot_nmin_ratio_stem.m 呈現。
    yoff = YL(1) - 0.022*diff(YL);                    % 端點只標數字、不畫 tick
    for xv = XL
        text(ax, xv, yoff, sprintf('%g',xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top','FontSize',FS,'FontWeight','bold','Clipping','off');
    end
    xlabel(ax, '$\mathbf{Sampling\;range\;(micro\;meter)}$', 'Interpreter','latex', 'FontSize',FS);
    ylabel(ax, ylab, 'Interpreter','latex', 'FontSize',FS);

    lg = legend(ax, [h1 h2], ...
         {sprintf('Reduced sampling (RMSPE = %.2f%%)', RMSPE), 'Full sampling'}, ...
         'Interpreter','tex', 'Location','northoutside', 'NumColumns',2);
    lg.FontSize = 24;  lg.FontWeight = 'bold';
    lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;
    hold(ax,'off');  drawnow;
    axp = get(ax,'Position');  lgh = get(lg,'Position');  lgh = lgh(4);
    GAPN = 0.022;  newTop = 1 - lgh - GAPN - 0.006;
    axp(4) = newTop - axp(2);  set(ax,'Position',axp);
    set(lg, 'Position', [axp(1), newTop + GAPN, axp(3), lgh]);

    out = fullfile(figdir, sprintf('%s_nmin_vs_full_maxwell10_single.png', lower(QTY)));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('\nwrote %s\n', out);
end

% ============================================================================
function [gA, gB] = gain_curves(here, RSEL, Nuse, S, F, force)
% ĝ_I 的兩條曲線：減量（每個 R 在 n=Nuse 抽 M 組、取中位數）與全取樣（既有 sweep）。
%   ⚠ 這是唯一需要重跑擬合的部分 —— 快取只存了 l_hat，沒存 ĝ_I。
%   點數都很小（222~2807），21 個 R × M 組，成本遠低於原本的 n-sweep。
    gcache = fullfile(here, 'data', 'nmin_gain_maxwell10_single.mat');
    gB = arrayfun(@(R) F.gI_R(F.Rum == R), RSEL);           % 全取樣（既有）
    if exist(gcache,'file') && ~force
        C = load(gcache);
        if isequal(C.RSEL(:).', RSEL(:).') && isequal(C.Nuse(:).', Nuse(:).')
            gA = C.gA;   fprintf('loaded gain cache %s\n', gcache);   return;
        end
        fprintf('  [warn] gain 快取的 (R, N_min) 與本次不同 → 重算\n');
    end

    solver_path();
    cfg = model_config('long2016_hexapole_halfcut','tip40um');
    raw = extract_maxwell_data(cfg, 'all', cfg.default_variant);
    ad  = build_actuator_data(raw, cfg);
    ad  = interp_grid_sample(ad, cfg, max(RSEL)*1e-6, 10e-6);
    rr  = sqrt(ad.r2);
    Fm  = zeros(6, cfg.N_I);
    for j = 1:cfg.N_I, Fm(cfg.apdl_to_paper_idx(j), j) = 1; end

    MREP = size(S.ellall{1}, 1);
    gA = nan(1, numel(RSEL));
    fprintf('\n  [gain] 逐 R 重跑（%d 組取中位數）\n', MREP);
    for t = 1:numel(RSEL)
        sel = find(rr <= RSEL(t)*1e-6);   np = numel(sel);
        gd  = nan(1, MREP);
        for d = 1:MREP
            rng(d-1);   ord = sel(randperm(np));            % 與 l_hat 用同一組 seed/排列
            idx = ord(1:Nuse(t));
            P   = ad.Pa(idx,:);
            Bs  = zeros(3*Nuse(t), cfg.N_I);
            for j = 1:cfg.N_I
                Bs(:,j) = reshape(ad.Ba(idx,:,j).', [], 1);
            end
            [e, l_hat] = fitting(P, Bs, ad.Pc_base, 0.5e-3, false);
            [~, g]     = solve_current(l_hat, e, ad.Pc_base, P, Bs, Fm);
            gd(d) = g;
        end
        gA(t) = median(gd, 'omitnan');
        fprintf('   R=%3d  n=%5d  gI_med=%.4f  [%.4f ~ %.4f]\n', RSEL(t), Nuse(t), gA(t), min(gd), max(gd));
    end
    save(gcache, 'RSEL','Nuse','gA');
    fprintf('saved gain cache %s\n', gcache);
end

% ============================================================================
function solver_path()
    APDL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    CAL  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    warning('off','MATLAB:rmpath:DirNotFound');
    rmpath(fullfile(APDL,'function'));  rmpath(fullfile(APDL,'common_path'));
    warning('on','MATLAB:rmpath:DirNotFound');
    addpath(fullfile(CAL,'function'));  addpath(fullfile(CAL,'common_path'));
end

% ============================================================================
function adq = interp_grid_sample(ad, cfg, Rmax, h)
    g = -Rmax:h:Rmax;
    [X, Y, Z] = ndgrid(g, g, g);
    Pq = [X(:), Y(:), Z(:)];
    rq = vecnorm(Pq, 2, 2);
    in = rq <= Rmax;   Pq = Pq(in,:);   rq = rq(in);
    Pm  = (ad.R_act.' * Pq.').';
    air = filter_iron_nodes(Pm(:,1), Pm(:,2), Pm(:,3) + cfg.SPH_OFST, cfg);
    Pq  = Pq(air,:);   rq = rq(air);
    src = ad.r2 < (Rmax + 60e-6)^2;
    Xs  = ad.Pa(src,:);
    N_I = size(ad.Ba,3);
    Bq  = zeros(size(Pq,1), 3, N_I);
    Fi  = scatteredInterpolant(Xs(:,1), Xs(:,2), Xs(:,3), ad.Ba(src,1,1), 'linear', 'none');
    for j = 1:N_I
        for c = 1:3
            Fi.Values = ad.Ba(src,c,j);
            Bq(:,c,j) = Fi(Pq);
        end
    end
    ok = ~any(any(isnan(Bq),3),2);
    fprintf('  [interp] 格距 %.0f um：球內 %d 點 → 濾鐵後 %d → 有效 %d\n', ...
            h*1e6, nnz(in), size(Pq,1), nnz(ok));
    adq = ad;
    adq.Pa = Pq(ok,:);   adq.r2 = rq(ok).^2;   adq.Ba = Bq(ok,:,:);
end

% ============================================================================
function k = flat_one(nlist, ell, TOLPC, KWIN)
% 單一條 l_hat(n) 的平坦段起點：第一個 a，使其後 KWIN 點的相對變化率都 <= TOLPC%。
    dch = [NaN, abs(diff(ell))./ell(2:end)*100];
    k = NaN;
    for a = 1:numel(nlist)-KWIN
        if all(dch(a+1 : a+KWIN) <= TOLPC), k = nlist(a);  return; end
    end
end

% ============================================================================
function tk = inner_ticks(lo, hi, n)
    cand = [1 2 2.5 3 4 5 10];
    s = (hi-lo)/(n+1);   k = floor(log10(s));
    [~, i] = min(abs(cand*10^k - s));   s = cand(i)*10^k;
    ctr = round(((lo+hi)/2)/s)*s;
    tk  = ctr + (-(n-1)/2 : (n-1)/2)*s;
    tk  = tk(tk > lo & tk < hi);
end

% ============================================================================
function [lim, tk] = axlim_auto(lo, hi, nlist)
% 奇數個等距 tick，兩端留白 = tick 間距（figure-style）。
    cand = [1 2 2.5 3 4 5 10];
    mid  = (lo+hi)/2;   rng_ = max(hi-lo, realmin);
    best = {};   bestSpan = inf;
    for n = nlist
        k0 = floor(log10(rng_/(n+1)));
        for k = k0:(k0+4)
            hit = false;
            for c = cand
                s = c*10^k;
                t = round(mid/s)*s + (-(n-1)/2 : (n-1)/2)*s;
                L = [t(1)-s, t(end)+s];   clr = 0.15*s;
                if lo >= L(1)+clr && hi <= L(2)-clr
                    if (n+1)*s < bestSpan, bestSpan = (n+1)*s;  best = {L, t}; end
                    hit = true;  break;
                end
            end
            if hit, break; end
        end
    end
    if isempty(best)
        n = nlist(1);  s = rng_/(n+1);
        t = mid + (-(n-1)/2 : (n-1)/2)*s;   best = {[t(1)-s, t(end)+s], t};
    end
    lim = best{1};   tk = best{2};
end
