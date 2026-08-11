function plot_nmin_vs_R(RLIST, TOLPC, NN, NCAP, KTAIL, ORDER, MREP, force)
% plot_nmin_vs_R -- paper 圖：每個取樣半徑 R 所需的最小取樣點數 N_min
% =========================================================================
%   水平軸 = 取樣半徑 R [um]、縱軸 = N_min（達到指定帶寬所需的最小點數）。
%
%   判準（使用者拍板 2026-08-11）：
%     ① 模型 = **single（1 參數）** —— 它是點數軸上較嚴的那個
%        （R=150：single 要 1588 點才進 0.26%，eighteen 20 點就進了）。
%     ② 帶寬 = **±TOLPC%（預設 0.26）**，基準 = 尾段**最大 KTAIL 個 n 的 l_hat 平均**。
%     ③ N_min = 最小的 n*，使得**所有 n >= n* 都留在帶內**（不是第一次進入帶內
%        —— 震盪期會偶然穿過，那樣會嚴重低估）。
%
%   [MODIFIED 2026-08-11] 基準由「最大 n 的單一值」改為「尾段 KTAIL 個的平均」。
%     為什麼：舊法的基準本身在擺動（R=150 尾段 n>=1588 仍在 874.0~876.3 游走，
%     峰對峰 0.26%）。基準隨機落在擺動的高點/低點 → 整條帶上下平移 → N_min 被推大/推小，
%     純粹是抽樣運氣。實測 R=380/400/420 跳 2018/10045/4502 即此故。
%     更嚴重的是**舊法在尾端自我指涉**：最後一點的偏差恆為 0、必在帶內，判準失去鑑別力
%     —— 所以舊法算出「N_min 佔池 64%」（R=120）那種值其實沒有意義，不只是不可信。
%     改用尾段平均後最後一點的偏差不再恆為 0；若連最大 n 都在帶外，回 NaN 並警告
%     （代表該 R 的尾段擺幅 >= 帶寬，判不出 N_min），**不給假數字**。
%
%   取樣：10 um 內插格、`rng(0)` 打亂一次、取前 n 個（巢狀累積，不重抽）。
%     ⚠ 效率關鍵：**10 um 格只在 R=max(RLIST) 內插一次**，各 R 只取 |p| <= R 的子集；
%       否則每個 R 重建一次 scatteredInterpolant 會慢數十倍。
%     ⚠ NCAP：每個 R 的最大 n 上限（預設 5e4）。R=150 的池 14,147 < NCAP，故先前
%       N_min = 1588 的結果可逐位重現。若某個 R 的 N_min 逼近 NCAP，代表基準本身
%       未收斂 → 該點不可信，須加大 NCAP 重跑（腳本會警告）。
%
%   RLIST 預設 100:20:500 —— 自 **R = 100 um** 起，因為那是「l_hat 與 g_I 同時穩定」
%   的最小半徑（見 ell_vs_npts_maxwell_interp.png：R<100 是 g_I 崩掉，非點數不足）。
%
%   風格 = 選項①粗體框 @ paper scale；**不放軸標題**（同 ell_vs_npts_* 系列）。
%   輸出 → figures/paper_fig/Section2_E/nmin_vs_R_maxwell10_single.png
% =========================================================================
    clc;
    if nargin < 1 || isempty(RLIST), RLIST = 100:20:500; end
    if nargin < 2 || isempty(TOLPC), TOLPC = 0.26;       end   % 帶寬 [%]
    if nargin < 3 || isempty(NN),    NN    = 40;         end   % 每個 R 測幾個點數
    if nargin < 4 || isempty(NCAP),  NCAP  = 5e4;        end   % 每個 R 的最大 n 上限
    if nargin < 5 || isempty(KTAIL), KTAIL = 5;          end   % 前瞻視窗長度（後 K 點）
    % [ADDED 2026-08-12] 取點順序 + 重複抽樣數：
    %   'random'（預設）＝每個 R 把池隨機打亂、巢狀取前 n → n 點**散布在整個 R 球**，
    %     取樣範圍**始終是 R** ⇒ 這才回答「在範圍 R 內最少取幾點」。
    %     ⚠ 單一組排列的 N_min 帶抽樣運氣（實測相鄰 R 可差 15 倍）→ 必須跑 MREP 組取**中位數**。
    %   'radial' ＝依半徑由內往外取。⚠ 該序下「少取點」等於「縮小範圍」，
    %     且曲線與 R 無關 ⇒ N_min(R) 必為水平線、且其值不是 R 的答案（誤差 R=500 達 9.4%）。
    if nargin < 6 || isempty(ORDER), ORDER = 'random';   end
    if nargin < 7 || isempty(MREP),  MREP  = 5;          end   % 隨機序的重複抽樣組數
    if nargin < 8 || isempty(force), force = false;      end

    HGRID = 10e-6;   N0 = 20;
    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    if strcmpi(ORDER,'radial')
        cachef = fullfile(here,'data','nmin_vs_R_maxwell10_single_radial_tol020.mat');   % 徑向：曲線與判準無關，共用一顆
    else
        cachef = fullfile(here,'data','nmin_vs_R_maxwell10_single_random.mat');   % 不含 M：可續抽
    end

    if strcmpi(ORDER,'random')
        % ===== 隨機序：快取可**續抽**（已有 seed 0..d0-1 就只補 d0..MREP-1）=====
        [RLIST, pool, nall, ellall] = random_draws(cachef, RLIST, MREP, NN, N0, NCAP, HGRID, force);
        [Nmin, Nlo, Nhi] = median_flat(nall, ellall, TOLPC, KTAIL);
        fprintf('\n  判準：後 %d 點變化率 <= %.2f%%；%d 組取中位數\n', KTAIL, TOLPC, MREP);
        for i = 1:numel(RLIST)
            fprintf('  R=%4d  pool=%7d   N_med=%7g   [min %g, max %g]\n', ...
                    RLIST(i), pool(i), Nmin(i), Nlo(i), Nhi(i));
        end
        % 從「判得出來」的第一個 R 開始畫（使用者要求：從點數夠的範圍開始）
        k0 = find(~isnan(Nmin), 1);
        if isempty(k0), error('所有 R 都判不出 N_flat'); end
        RLIST = RLIST(k0:end);  Nmin = Nmin(k0:end);  pool = pool(k0:end);
    elseif exist(cachef,'file') && ~force
        % 徑向：快取存的是**完整曲線**，載入後重算判準（曲線與 R 無關，各 R 只是截斷位置不同）。
        S = load(cachef);   pool = S.pool;   RLIST = S.RLIST;
        fprintf('loaded cache %s（重算判準：後 %d 點變化率 <= %.2f%%）\n', cachef, KTAIL, TOLPC);
        [Nmin, reffm] = flat_start(S.nall{1}, S.ellall{1}, S.ellref, pool, TOLPC, KTAIL);
        for i = 1:numel(RLIST)
            fprintf('  R=%4d  pool=%7d   N_flat=%8g   r_eff=%6.1f um\n', RLIST(i), pool(i), Nmin(i), reffm(i));
        end
    else
        solver_path();
        cfg = model_config('long2016_hexapole_halfcut','tip40um');
        raw = extract_maxwell_data(cfg, 'all', cfg.default_variant);
        ad  = build_actuator_data(raw, cfg);
        ad  = interp_grid_sample(ad, cfg, max(RLIST)*1e-6, HGRID);   % **只做一次**
        rr  = sqrt(ad.r2);

        % =====================================================================
        % [ADDED 2026-08-11] 徑向排序模式（使用者拍板）：由內往外取點。
        %   判準改為 **逐點變化率**：|Δl_hat|/l_hat <= TOLPC，且其後不再超出。
        %   （不用「偏離終值」判準 —— 那在徑向序下是自我指涉的：終點就是 r=R 的答案。）
        %
        %   🔑 效率關鍵：徑向序下 l_hat(r) **與 R 無關** —— R=150 的序列就是 R=500 序列的
        %   前段（巢狀）。所以只跑**一次** R=max 的完整掃描，各 R 只把曲線截斷到 n <= N(R)
        %   再套判準，不必跑 21 次。
        % =====================================================================
        if strcmpi(ORDER,'random')
            % ===== 隨機序：每個 R 抽 MREP 組獨立排列，各算一次 N_flat，取中位數 =====
            nR = numel(RLIST);   pool = zeros(1,nR);
            nall = cell(1,nR);   ellall = cell(1,nR);
            fprintf('\n    R    pool   n_max   ');
            for d = 1:MREP, fprintf('  d%d   ', d); end
            fprintf('  ->  中位數\n');
            for i = 1:nR
                sel = find(rr <= RLIST(i)*1e-6);   pool(i) = numel(sel);
                nmax  = min(pool(i), NCAP);
                nlist = unique(round(logspace(log10(N0), log10(nmax), NN)));
                E = nan(MREP, numel(nlist));
                for d = 1:MREP
                    rng(d-1);   ord = sel(randperm(pool(i)));   % 每組獨立排列（seed 0..M-1）
                    for k = 1:numel(nlist)
                        idx = ord(1:nlist(k));                  % 巢狀累積：只增不換
                        P   = ad.Pa(idx,:);
                        Bs  = zeros(3*nlist(k), cfg.N_I);
                        for j = 1:cfg.N_I
                            Bs(:,j) = reshape(ad.Ba(idx,:,j).', [], 1);
                        end
                        [~, l_hat] = fitting(P, Bs, ad.Pc_base, 0.5e-3, false);   % single
                        E(d,k) = l_hat*1e6;
                    end
                end
                nall{i} = nlist;   ellall{i} = E;
                nd = arrayfun(@(d) flat_one(nlist, E(d,:), TOLPC, KTAIL), 1:MREP);
                fprintf('  %4d  %6d  %6d  ', RLIST(i), pool(i), nmax);
                fprintf('%6g ', nd);   fprintf('  ->  %7g\n', median(nd,'omitnan'));
            end
            [Nmin, Nlo, Nhi] = median_flat(nall, ellall, TOLPC, KTAIL);
            save(cachef, 'RLIST','Nmin','Nlo','Nhi','pool','nall','ellall', ...
                         'TOLPC','NN','NCAP','KTAIL','MREP','HGRID','ORDER');
            fprintf('saved cache %s\n', cachef);
        else
        Rmax = max(RLIST);
        nfull = unique(round(logspace(log10(N0), log10(numel(rr)), NN)));
        [~, ordr] = sort(rr, 'ascend');                 % 由內往外（確定性）
        ellf = nan(size(nfull));   refff = nan(size(nfull));
        fprintf('\n  [徑向全掃描 @R=%d]      n     ell_hat[um]   r_eff[um]\n', Rmax);
        for k = 1:numel(nfull)
            idx = ordr(1:nfull(k));
            P   = ad.Pa(idx,:);
            Bs  = zeros(3*nfull(k), cfg.N_I);
            for j = 1:cfg.N_I
                Bs(:,j) = reshape(ad.Ba(idx,:,j).', [], 1);
            end
            [~, l_hat] = fitting(P, Bs, ad.Pc_base, 0.5e-3, false);   % single
            ellf(k)  = l_hat*1e6;
            refff(k) = sqrt(max(ad.r2(idx)))*1e6;
            fprintf('                     %7d   %10.3f   %7.1f\n', nfull(k), ellf(k), refff(k));
        end
        dch = [NaN, abs(diff(ellf))./ellf(2:end)*100];   % 逐點變化率 [%]

        nR = numel(RLIST);   pool = zeros(1,nR);
        for i = 1:nR, pool(i) = nnz(rr <= RLIST(i)*1e-6); end
        [Nmin, reffm] = flat_start(nfull, ellf, refff, pool, TOLPC, KTAIL);   % 平坦段起點
        fprintf('\n    R    pool    N_flat(後%d點Δ<=%.2f%%)   r_eff\n', KTAIL, TOLPC);
        for i = 1:nR
            fprintf('  %4d  %6d   %14g        %6.1f\n', RLIST(i), pool(i), Nmin(i), reffm(i));
        end
        ellall = {ellf};  nall = {nfull};  ellref = refff;  tailpp = dch;   % 存全掃描供診斷
        save(cachef, 'RLIST','Nmin','pool','reffm','ellall','nall','ellref','tailpp', ...
                     'TOLPC','NN','NCAP','KTAIL','HGRID');
        fprintf('saved cache %s\n', cachef);
        end   % ORDER 分支
    end
    if false   % --- 以下為舊的 random-order 分支，保留備查、目前不執行 ---
        nR = numel(RLIST);
        Nmin = nan(1,nR);  pool = zeros(1,nR);  ellref = nan(1,nR);  tailpp = nan(1,nR);
        ellall = cell(1,nR);  nall = cell(1,nR);      % [ADDED] 存每個 R 的完整 l_hat(n) 供診斷
        fprintf('\n    R    pool    n_max   ell_ref[um]  尾段擺幅%%   N_min(+-%.2f%%)\n', TOLPC);
        for i = 1:nR
            sel  = find(rr <= RLIST(i)*1e-6);
            pool(i) = numel(sel);
            nmax = min(pool(i), NCAP);
            nlist = unique(round(logspace(log10(N0), log10(nmax), NN)));

            rng(0);   ord = sel(randperm(pool(i)));        % 該 R 的池打亂一次
            ell = nan(size(nlist));
            for k = 1:numel(nlist)
                idx = ord(1:nlist(k));                     % 前 n 個 → 巢狀累積
                P   = ad.Pa(idx,:);
                Bs  = zeros(3*nlist(k), cfg.N_I);
                for j = 1:cfg.N_I
                    Bs(:,j) = reshape(ad.Ba(idx,:,j).', [], 1);
                end
                [~, l_hat] = fitting(P, Bs, ad.Pc_base, 0.5e-3, false);   % single
                ell(k) = l_hat*1e6;
            end
            ellall{i} = ell;   nall{i} = nlist;      % [ADDED] 診斷用
            % 基準 = 尾段 KTAIL 個的平均（消掉基準自身的擺動；見檔頭 [MODIFIED]）
            kt = min(KTAIL, numel(nlist));
            ellref(i) = mean(ell(end-kt+1:end));
            tailpp(i) = (max(ell(end-kt+1:end)) - min(ell(end-kt+1:end))) / ellref(i) * 100;

            dev = abs(ell - ellref(i))/ellref(i)*100;
            bad = find(dev > TOLPC, 1, 'last');
            flag = '';
            if isempty(bad)
                Nmin(i) = nlist(1);
            elseif bad < numel(nlist)
                Nmin(i) = nlist(bad+1);
            else
                Nmin(i) = NaN;                    % 連最大 n 都在帶外 → 判不出來
                flag = '  <-- 連 n_max 都越界，尾段擺幅 >= 帶寬，N_min 無定義';
            end
            if ~isnan(Nmin(i)) && Nmin(i) > 0.5*nmax
                flag = '  <-- N_min > 半個池，該 R 的池不夠大，值不可用';
            end
            fprintf('  %4d  %6d  %6d   %9.3f   %8.3f   %7g%s\n', ...
                    RLIST(i), pool(i), nmax, ellref(i), tailpp(i), Nmin(i), flag);
        end
        save(cachef, 'RLIST','Nmin','pool','ellref','tailpp','ellall','nall', ...
                     'TOLPC','NN','NCAP','KTAIL','HGRID');
        fprintf('saved cache %s\n', cachef);
    end

    % ---- 畫圖（單 panel，選項①粗體框；無軸標題）----
    FS = 36;  LWBOX = 4.0;  col = [0.05 0.10 0.95];
    fig = figure('Color','w','Position',[100 100 1120 820]);
    ax  = axes(fig);  hold(ax,'on');
    h1 = plot(ax, RLIST, Nmin, '-o', 'Color',col, 'LineWidth',3.0, ...
              'MarkerSize',8, 'MarkerFaceColor',col, 'Clipping','off');
    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX,'TickLength',[.015 .015]);
    ax.Toolbar.Visible = 'off';

    XL = [RLIST(1) RLIST(end)];   xlim(ax, XL);
    XT = inner_ticks(XL(1), XL(2), 3);   set(ax,'XTick',XT);
    if all(isnan(Nmin)), error('所有 R 都判不出 N_flat'); end
    ylo = min(Nmin,[],'omitnan');   yhi = max(Nmin,[],'omitnan');
    if yhi <= ylo                    % [ADDED] 全部相同（徑向序下常見）→ 人為撐開，否則 ylim 報錯
        ylo = ylo*0.8;   yhi = yhi*1.2;
    end
    [YL, YT] = axlim_auto(ylo, yhi, [3 5]);
    ylim(ax, YL);   set(ax,'YTick',YT);

    yoff = YL(1) - 0.022*diff(YL);                  % 端點只標數字、不畫 tick
    for xv = XL
        text(ax, xv, yoff, sprintf('%g',xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top','FontSize',FS,'FontWeight','bold','Clipping','off');
    end

    % [MODIFIED 2026-08-11] 使用者要求加回軸標題（照 figure-style：LaTeX \mathbf、FS 36、單位括號）
    xlabel(ax, '$\mathbf{Sampling\;range\;(micro\;meter)}$', 'Interpreter','latex', 'FontSize',FS);
    ylabel(ax, '$\mathbf{Minimum\;number\;of\;points}$',     'Interpreter','latex', 'FontSize',FS);

    if strcmpi(ORDER,'radial'), otxt = 'radial order';
    else,                       otxt = sprintf('median of %d draws', MREP);   end
    lg = legend(ax, h1, {sprintf('Single parameter, %s, next %d steps |\\Delta| \\leq %.2f%%', otxt, KTAIL, TOLPC)}, ...
                'Interpreter','tex', 'Location','northoutside', 'NumColumns',1);
    lg.FontSize = 24;  lg.FontWeight = 'bold';
    lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;
    hold(ax,'off');   drawnow;
    axp = get(ax,'Position');  lgh = get(lg,'Position');  lgh = lgh(4);
    GAPN = 0.022;  newTop = 1 - lgh - GAPN - 0.006;
    axp(4) = newTop - axp(2);  set(ax,'Position',axp);
    set(lg, 'Position', [axp(1), newTop + GAPN, axp(3), lgh]);

    out = fullfile(figdir, 'nmin_vs_R_maxwell10_single.png');
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('\nwrote %s\n', out);
end

% ============================================================================
function [RLIST, pool, nall, ellall] = random_draws(cachef, RLIST, MREP, NN, N0, NCAP, HGRID, force)
% 隨機序的抽樣主體，**支援續抽**：快取已有 d0 組（seed 0..d0-1）就只補 d0..MREP-1。
%   每組：把該 R 的池用 seed 打亂 → 巢狀取前 n → 擬合 → 一條 l_hat(n)。
%   回傳 nall{i}=該 R 的 n 清單、ellall{i}= MREP × numel(n) 的 l_hat 矩陣。
    nR = numel(RLIST);   pool = zeros(1,nR);
    nall = cell(1,nR);   ellall = cell(1,nR);   d0 = 0;

    if exist(cachef,'file') && ~force
        S = load(cachef);
        if isequal(S.RLIST(:).', RLIST(:).')
            nall = S.nall;   ellall = S.ellall;   pool = S.pool;
            d0   = size(ellall{1}, 1);
            fprintf('loaded cache %s（已有 %d 組）\n', cachef, d0);
            if d0 >= MREP
                for i = 1:nR, ellall{i} = ellall{i}(1:MREP,:); end
                return;
            end
            fprintf('  → 續抽 seed %d..%d（共補 %d 組）\n', d0, MREP-1, MREP-d0);
        else
            fprintf('  [warn] 快取的 RLIST 與本次不同 → 全部重算\n');   d0 = 0;
        end
    end

    solver_path();
    cfg = model_config('long2016_hexapole_halfcut','tip40um');
    raw = extract_maxwell_data(cfg, 'all', cfg.default_variant);
    ad  = build_actuator_data(raw, cfg);
    ad  = interp_grid_sample(ad, cfg, max(RLIST)*1e-6, HGRID);
    rr  = sqrt(ad.r2);

    for i = 1:nR
        sel = find(rr <= RLIST(i)*1e-6);   pool(i) = numel(sel);
        if d0 == 0
            nall{i} = unique(round(logspace(log10(N0), log10(min(pool(i),NCAP)), NN)));
            ellall{i} = nan(0, numel(nall{i}));
        end
        nlist = nall{i};
        E = ellall{i};
        for d = d0+1 : MREP
            rng(d-1);   ord = sel(randperm(pool(i)));      % seed = d-1
            row = nan(1, numel(nlist));
            for k = 1:numel(nlist)
                idx = ord(1:nlist(k));                     % 巢狀累積：只增不換
                P   = ad.Pa(idx,:);
                Bs  = zeros(3*nlist(k), cfg.N_I);
                for j = 1:cfg.N_I
                    Bs(:,j) = reshape(ad.Ba(idx,:,j).', [], 1);
                end
                [~, l_hat] = fitting(P, Bs, ad.Pc_base, 0.5e-3, false);   % single
                row(k) = l_hat*1e6;
            end
            E(d,:) = row;
        end
        ellall{i} = E;
        fprintf('  R=%4d  pool=%7d  已完成 %d 組\n', RLIST(i), pool(i), MREP);
    end
    save(cachef, 'RLIST','pool','nall','ellall','NN','NCAP','HGRID');
    fprintf('saved cache %s（%d 組）\n', cachef, MREP);
end

% ============================================================================
function k = flat_one(nlist, ell, TOLPC, KWIN)
% 單一條 l_hat(n) 曲線的「平坦段起點」：第一個 k，使其後 KWIN 點的相對變化率都 <= TOLPC%。
    dch = [NaN, abs(diff(ell))./ell(2:end)*100];
    k = NaN;
    for a = 1:numel(nlist)-KWIN
        if all(dch(a+1 : a+KWIN) <= TOLPC), k = nlist(a);  return; end
    end
end

% ============================================================================
function [Nmed, Nlo, Nhi] = median_flat(nall, ellall, TOLPC, KWIN)
% 各 R：對 MREP 組排列各算一次 flat_one，回中位數與全距。
    nR = numel(nall);   Nmed = nan(1,nR);  Nlo = nan(1,nR);  Nhi = nan(1,nR);
    for i = 1:nR
        E  = ellall{i};                                  % MREP × numel(nlist)
        nd = arrayfun(@(d) flat_one(nall{i}, E(d,:), TOLPC, KWIN), 1:size(E,1));
        Nmed(i) = median(nd, 'omitnan');
        Nlo(i)  = min(nd);   Nhi(i) = max(nd);
    end
end

% ============================================================================
function [Nmin, reffm] = flat_start(nfull, ellf, refff, pool, TOLPC, KWIN)
% 「平坦段的第一個點」（使用者拍板 2026-08-12）：
%   對候選點 k，看**其後 KWIN 個點**（k+1..k+KWIN，不含 k 自己）的相對變化率
%   |Δl_hat|/l_hat；若這 KWIN 個全部 <= TOLPC%，則 k 就是平坦段起點。取最小的 k。
%   ⚠ 這是「前瞻視窗」判準 —— 只看接下來幾點，不要求「其後永遠不超出」，
%     所以大 R 時不會因為外圈的系統性下滑而變成無解。
    dch  = [NaN, abs(diff(ellf))./ellf(2:end)*100];      % 逐點變化率 [%]
    nR   = numel(pool);
    Nmin = nan(1,nR);   reffm = nan(1,nR);
    for i = 1:nR
        m  = find(nfull <= pool(i));                     % 截斷到該 R 的池
        for k = 1:numel(m)-KWIN
            w = dch(m(k+1 : k+KWIN));
            if all(w <= TOLPC)
                Nmin(i)  = nfull(m(k));
                reffm(i) = refff(m(k));
                break;
            end
        end
    end
end

% ============================================================================
function tk = inner_ticks(lo, hi, n)
% n 個等距內部刻度（不含端點；n 取奇數）。
    cand = [1 2 2.5 3 4 5 10];
    s = (hi-lo)/(n+1);   k = floor(log10(s));
    [~, i] = min(abs(cand*10^k - s));   s = cand(i)*10^k;
    ctr = round(((lo+hi)/2)/s)*s;
    tk  = ctr + (-(n-1)/2 : (n-1)/2)*s;
    tk  = tk(tk > lo & tk < hi);
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
% 把場內插到 h 間距均勻格、只留球內（與 plot_ell_gain_vs_R 同名 local 同演算法）。
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
    fprintf('  [interp] 格距 %.0f um：球內 %d 點 → 濾鐵後 %d → 有效 %d（源點 %d）\n', ...
            h*1e6, nnz(in), size(Pq,1), nnz(ok), nnz(src));
    adq = ad;
    adq.Pa = Pq(ok,:);   adq.r2 = rq(ok).^2;   adq.Ba = Bq(ok,:,:);
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
                s   = c*10^k;
                T   = { round(mid/s)*s + (-(n-1)/2 : (n-1)/2)*s };
                if lo >= 0 && lo < s, T = [{(0:n-1)*s}, T]; end %#ok<AGROW>
                for it = 1:numel(T)
                    t = T{it};   L = [t(1)-s, t(end)+s];   clr = 0.15*s;
                    if lo >= L(1)+clr && hi <= L(2)-clr
                        if (n+1)*s < bestSpan, bestSpan = (n+1)*s;  best = {L, t}; end
                        hit = true;  break;
                    end
                end
                if hit, break; end
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
