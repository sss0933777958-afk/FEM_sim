function plot_conv_npts(force)
% plot_conv_npts -- R=150um：l_hat、g_I 隨取樣點數的收斂（兩張圖，判到收斂就收手）
% =========================================================================
%   沿「按比例輪流 +1」的等測度網格階梯（起點 [1 2 3]、權重 1 : 3 : 3pi）逐個設計擬合，
%   同時追兩件事：
%
%   圖 A  l_hat vs 點數                     → ell_vs_npts_conv_maxwell_R150.png
%         含 N=0 錨點（值 = 初值 l0 = 500 um，畫在 log 軸最左端、端點標 0）。
%         收斂判準：第一個設計 i，其後**連續 10 步**「相對前一點」變化率皆 < 0.5%
%         （視窗 KELL；10 段變化率 ⇒ 涵蓋 v(i)..v(i+10) 共 11 個設計點）。
%         [MODIFIED 2026-08-17] 使用者拍板：**圖上不再畫收斂點虛線**（判定照舊）。
%
%   圖 B  [REMOVED 2026-08-17] K_I 的 Frobenius 相對變化 vs 點數（使用者拍板不需要）。
%         縱軸曾是 ||K_i - K_{i-1}||_F / ||K_{i-1}||_F * 100 [%]，收斂判準為
%         「第一個 i，其後連續 10 步該值 < 0.5%」（視窗 KKI）。
%         ⚠ **計算保留**（console 報告 + 快取 S.fro1/fro2 + sweep 早停的六判準之一），
%         只是不再出圖 —— 拿掉計算會改變早停時機、縮短快取涵蓋的設計範圍。
%
%   圖 C  g_I_hat vs 點數                   → gain_vs_npts_conv_maxwell_R150.png
%         [ADDED 2026-08-13] 收斂判準與 l_hat 同一把尺（後 10 步相對變化 < 0.5%）。
%         [MODIFIED 2026-08-17] 使用者拍板：**圖上不再畫收斂點虛線**（判定照舊）。
%         無 N=0 錨點（g_I 沒有初值）；小 N 端近簡併值離譜，只畫 N >= 20。
%
%   ⚠ K_I 的號誌條件（off-diag 全負、對角全正）**不列入本腳本的判準**：本腳本圖 B 用的是
%     量化指標。號誌狀態仍算出來存進快取供追溯。
%     ⚠ 注意這與 plot_full_vs_conv_vs_R.m 現行的判準**不同**：那支已改成
%     「l_hat 穩定 ∧ g_I 穩定 ∧ K_I 符合物理（號誌閘門）」三者交集，因為 K_I 的
%     Frobenius 判準對 single 不可達（詳見該檔 conv_fit 註解）。兩支要不要統一未定。
%
%   早停：每加一個設計就檢查兩模型的兩個判準是否都成立；成立後**再多跑 5 個設計當緩衝**
%   即停止。先前 200 個設計全跑的版本，九成時間花在對收斂點毫無貢獻的大設計上。
%
%   風格①粗體框圖：FS 36 粗體、box on、grid off、log 橫軸、刻度等比且兩端留 15% 淨空、
%   曲線首末點貼齊框邊、起訖數字以 text 補、圖例照 figure-style 標準樣式。
%   配色：single 藍、eighteen 紅（顏色 = 模型）。
% =========================================================================
    clc;
    if nargin < 1 || isempty(force), force = false; end

    R      = 150e-6;                 % 取樣半徑 [m]
    l0     = 0.5e-3;                 % l_hat 初值 [m]（= N=0 錨點的值）
    TOL    = 0.005;                  % 兩個判準共用的容忍度：相對變化 < 0.5%
    % [MODIFIED 2026-08-13] 視窗 5 → 10（使用者拍板）。理由：K=5..7 對三個量的收斂點
    %   **完全相同**（實測 ell 8/6、gain 528/10、K_I 1200/88 逐字不變）；K=10 是最小
    %   會動的值，且只動 eighteen 的 gain（10 → 80，離全格點 +0.49% → −0.10%）。
    %   ⚠ single 的 gain 不論 K 多長都停在 528：它從 528 之後每步都 <0.5%、但持續同向
    %   漂移（累積 +1.02%），逐步變化率判準對慢漂移沒有鑑別力。
    KELL   = 10;                     % l_hat 與 g_I 判準的視窗（後 10 步）
    KKI    = 10;                     % K_I Frobenius 判準的視窗（後 10 步）
    BUFFER = 5;                      % 各判準都成立後再多跑幾個設計當緩衝
    NDMAX  = 120;                    % 保險上限（正常會早停，用不到）
    MODEL  = 'long2016_hexapole_halfcut';   GEOM = 'tip40um';

    here   = fileparts(fileparts(mfilename('fullpath')));       % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    cachef = fullfile(here, 'data', 'conv_npts_maxwell_R150.mat');

    %% ---- 計算（快取）------------------------------------------------------
    if exist(cachef,'file') && ~force
        S = load(cachef);   fprintf('由快取載入 %s\n', cachef);
    else
        S = sweep(MODEL, GEOM, R, l0, TOL, KELL, KKI, BUFFER, NDMAX);
        save(cachef, '-struct', 'S');   fprintf('已存 %s\n', cachef);
    end

    % [ADDED] 收斂索引一律用**當下的判準**從快取資料重新判定，不沿用快取裡存的舊值
    %   ⇒ 之後調 TOL 或視窗長度都不必重跑擬合（擬合結果與判準是分開的兩件事）。
    S.iell1 = first_stable(S.ell1, TOL, KELL);   S.iell2 = first_stable(S.ell2, TOL, KELL);
    S.igI1  = first_stable(S.gI1,  TOL, KELL);   S.igI2  = first_stable(S.gI2,  TOL, KELL);
    S.ifro1 = first_below(S.fro1, TOL*100, KKI); S.ifro2 = first_below(S.fro2, TOL*100, KKI);

    %% ---- console 報告 -----------------------------------------------------
    fprintf('\n  N        : %s\n', num2str(S.N,    '%8d'));
    fprintf('  l_hat 1p : %s um\n', num2str(S.ell1,'%8.1f'));
    fprintf('  l_hat 18p: %s um\n', num2str(S.ell2,'%8.1f'));
    fprintf('  g_I   1p : %s mT/A\n', num2str(S.gI1,'%8.3f'));
    fprintf('  g_I   18p: %s mT/A\n', num2str(S.gI2,'%8.3f'));
    fprintf('  dK/K 1p  : %s %%\n', num2str([NaN S.fro1(2:end)],'%8.2f'));
    fprintf('  dK/K 18p : %s %%\n', num2str([NaN S.fro2(2:end)],'%8.2f'));
    rep('l_hat  single  ', S.N, S.iell1, S.TRI);
    rep('l_hat  eighteen', S.N, S.iell2, S.TRI);
    rep('g_I    single  ', S.N, S.igI1,  S.TRI);
    rep('g_I    eighteen', S.N, S.igI2,  S.TRI);
    rep('K_I    single  ', S.N, S.ifro1, S.TRI);
    rep('K_I    eighteen', S.N, S.ifro2, S.TRI);

    %% ---- 圖 A：l_hat（含 N=0 錨點）---------------------------------------
    % [MODIFIED 2026-08-14] 使用者拍板：**只畫到「收斂點之後再 3 個設計」**、橫軸改**線性**。
    %   理由：收斂點在 N=6/8，後面到 N=1776 全是平的，log 軸把九成篇幅給了沒有資訊的尾段。
    %   截到收斂後 3 個設計（索引 5 → N=14）＝「跳一下就穩住」再走三步確認，畫面最緊湊。
    %   （曾試 +10 → N=88、+7 → N=30，右側都偏空。）圖 B / 圖 C 不動（仍 log、完整階梯）。
    ic  = min(numel(S.N), max([S.iell1, S.iell2]) + 3);
    NA  = [0,          S.N(1:ic)];
    A1  = [l0*1e6,     S.ell1(1:ic)];
    A2  = [l0*1e6,     S.ell2(1:ic)];
    fprintf('\n[圖A] 截到索引 %d（N=%d）＝收斂點後 3 個設計\n', ic, S.N(ic));
    % [MODIFIED 2026-08-17] 使用者拍板：本圖**不畫收斂點虛線**（Nc 傳 NaN 即不畫）。
    %   收斂點仍照常判定並印在 console / 存進快取，只是不畫在圖上。
    mk('ell', NA, A1, A2, NaN, NaN, ...
       '$\mathbf{\hat{\ell}\;(micro\;meter)}$', true, figdir, [], 'linear');

    %% ---- 圖 B：K_I 的 Frobenius 相對變化 —— [REMOVED 2026-08-17] -----------
    % 使用者拍板：這張圖不需要了，`kifrob_vs_npts_conv_maxwell_R150.png` 已刪除。
    % ⚠ Frobenius 的**計算**刻意保留（`S.fro1/fro2`、`S.ifro1/ifro2`、console 報告、
    %   sweep 早停的六判準之一）—— 拿掉會改變早停時機、讓快取涵蓋的設計範圍縮短。
    % 要恢復這張圖，把下面兩行取消註解即可：
    %   m = 2:numel(S.N);
    %   mk('kifrob', S.N(m), S.fro1(m), S.fro2(m), idx2N(S.N,S.ifro1), idx2N(S.N,S.ifro2), ...
    %      '$\mathbf{\|\Delta \bar{K}_{I}\|_{F}/\|\bar{K}_{I}\|_{F}\;(\%)}$', false, figdir);

    %% ---- 圖 C：g_I_hat（無 N=0 錨點：g_I 沒有初值）------------------------
    % [ADDED 2026-08-13] 判準加入 g_I 後補這張。**從階梯第一個設計（N=6）畫起**，
    %   eighteen 的收斂點 N=10 才進得了圖。y 刻度使用者指定為整數 → 明給 [6 8 10]。
    % [MODIFIED 2026-08-14] **截到「收斂點之後再 3 個設計」**（索引 31 → N=600）；
    %   橫軸維持**對數**、主刻度以 **10^n 指數形式**呈現（使用者拍板，'logexp'）。
    %   （曾試線性軸：N=6~30 的 9 個點會壓在最左緣、single 的震盪看不清 → 還原回對數。）
    icg = min(numel(S.N), max([S.igI1, S.igI2]) + 3);
    fprintf('[圖C] 截到索引 %d（N=%d）＝收斂點後 3 個設計\n', icg, S.N(icg));
    % [MODIFIED 2026-08-17] 同圖 A：**不畫收斂點虛線**。
    mk('gain', S.N(1:icg), S.gI1(1:icg), S.gI2(1:icg), NaN, NaN, ...
       '$\mathbf{{}^{B}\hat{g}_{I}\;(mT/A)}$', false, figdir, [6 8 10], 'logexp');
end

% ============================================================================
function S = sweep(MODEL, GEOM, R, l0, TOL, KELL, KKI, BUFFER, NDMAX)
% 逐設計擬合，早停。回傳 N、l_hat、K_I、Frobenius 變化率與四個收斂索引。
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    addpath(fullfile(CAL,'function'), fullfile(CAL,'common_path'), fullfile(CAL,'utils'));
    addpath(fullfile(CAL,'utils','long2016_hexapole_halfcut'));

    cfg = model_config(MODEL, GEOM);
    F   = zeros(6, cfg.N_I);
    for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end

    % 設計階梯：按比例輪流 +1
    w = [1, 3, 3*pi];   t = [1 2 3];   TRI = zeros(NDMAX,3);
    for q = 1:NDMAX, TRI(q,:) = t;  [~, j] = min(t ./ w);  t(j) = t(j) + 1;  end

    [N, ell1, ell2, fro1, fro2, gI1, gI2] = deal(nan(1,NDMAX));
    K1 = nan(6,6,NDMAX);   K2 = nan(6,6,NDMAX);
    [ki1, ki2] = deal(false(1,NDMAX));                  % 號誌狀態（存查、不入判準）
    done_at = NaN;

    for q = 1:NDMAX
        [x,y,z,B] = sphere_grid_sample(R, [], struct('frame','actuator','NRPT',TRI(q,:)));
        P  = [x y z];   N(q) = size(P,1);
        Bs = zeros(3*N(q), size(B,3));
        for j = 1:size(B,3), Bs(:,j) = reshape(B(:,:,j).', [], 1); end

        [ell1(q), K1(:,:,q), ki1(q), gI1(q)] = fit_one(P, Bs, cfg.Pc_base, l0, false, F);
        [ell2(q), K2(:,:,q), ki2(q), gI2(q)] = fit_one(P, Bs, cfg.Pc_base, l0, true,  F);
        if q > 1
            fro1(q) = frob_rel(K1(:,:,q), K1(:,:,q-1));
            fro2(q) = frob_rel(K2(:,:,q), K2(:,:,q-1));
        end
        fprintf('  N=%6d (%2d,%2d,%3d)：1p l=%7.1f g=%7.3f dK/K=%7.3f%% ki=%d｜18p l=%7.1f g=%7.3f dK/K=%7.3f%% ki=%d\n', ...
                N(q), TRI(q,:), ell1(q), gI1(q), fro1(q), ki1(q), ell2(q), gI2(q), fro2(q), ki2(q));

        % --- 早停：六個判準都成立後再跑 BUFFER 個設計 ---
        %   [MODIFIED 2026-08-13] 加入 g_I 的收斂判準（與 l_hat 同一把尺）。
        if isnan(done_at)
            i1 = first_stable(ell1(1:q), TOL, KELL);
            i2 = first_stable(ell2(1:q), TOL, KELL);
            g1 = first_stable(gI1(1:q),  TOL, KELL);
            g2 = first_stable(gI2(1:q),  TOL, KELL);
            j1 = first_below(fro1(1:q), TOL*100, KKI);
            j2 = first_below(fro2(1:q), TOL*100, KKI);
            if all(~isnan([i1 i2 g1 g2 j1 j2])), done_at = q; end
        end
        if ~isnan(done_at) && q >= done_at + BUFFER, break; end
    end
    nq = q;

    S.N = N(1:nq);   S.ell1 = ell1(1:nq);   S.ell2 = ell2(1:nq);
    S.gI1 = gI1(1:nq);     S.gI2 = gI2(1:nq);
    S.fro1 = fro1(1:nq);   S.fro2 = fro2(1:nq);
    S.K1 = K1(:,:,1:nq);   S.K2 = K2(:,:,1:nq);
    S.ki1 = ki1(1:nq);     S.ki2 = ki2(1:nq);     S.TRI = TRI(1:nq,:);
    S.iell1 = first_stable(S.ell1, TOL, KELL);   S.iell2 = first_stable(S.ell2, TOL, KELL);
    S.igI1  = first_stable(S.gI1,  TOL, KELL);   S.igI2  = first_stable(S.gI2,  TOL, KELL);
    S.ifro1 = first_below(S.fro1, TOL*100, KKI); S.ifro2 = first_below(S.fro2, TOL*100, KKI);
    S.R = R;  S.TOL = TOL;  S.KELL = KELL;  S.KKI = KKI;  S.BUFFER = BUFFER;
end

% ============================================================================
function [ell_um, K, kiok, gI] = fit_one(P, Bstack, Pc_base, l0, USE_BIAS, F)
% 擬合一個設計 → l_hat [um]、K_I、g_I_hat [mT/A]；kiok 是號誌狀態。
    try
        [e, l]  = fitting(P, Bstack, Pc_base, l0, USE_BIAS);
        [K, gI] = solve_current(l, e, Pc_base, P, Bstack, F);
        ell_um = l*1e6;
        M = ~eye(6);   [~, im] = max(abs(K), [], 2);
        kiok = all(diag(K) > 0) && all(K(M) < 0) && isequal(im(:).', 1:6);
    catch ME
        warning('fit_one:degenerate','N=%d 擬合失敗（%s）', size(P,1), ME.identifier);
        ell_um = NaN;   K = nan(6);   kiok = false;   gI = NaN;
    end
end

% ============================================================================
function r = frob_rel(Ka, Kb)
% ||Ka - Kb||_F / ||Kb||_F * 100 [%]
    if any(~isfinite(Ka(:))) || any(~isfinite(Kb(:))), r = NaN;  return; end
    r = 100 * norm(Ka - Kb, 'fro') / norm(Kb, 'fro');
end

% ============================================================================
function i0 = first_stable(v, tol, K)
% 第一個 i，使其後連續 K 個「相對前一點」的變化率都 < tol。未達判準回 NaN。
    i0  = NaN;
    if numel(v) < K+1, return; end
    rel = abs(diff(v) ./ v(1:end-1));
    for i = 1:numel(rel)-K+1
        if all(rel(i:i+K-1) < tol), i0 = i;  return; end
    end
end

% ============================================================================
function i0 = first_below(f, tolpct, K)
% 第一個 i，使 f(i..i+K-1) 都 < tolpct（f 已是百分比、f(1)=NaN）。未達判準回 NaN。
    i0 = NaN;
    for i = 2:numel(f)-K+1
        seg = f(i:i+K-1);
        if all(isfinite(seg)) && all(seg < tolpct), i0 = i;  return; end
    end
end

% ============================================================================
function Nc = idx2N(N, i0)
    if isnan(i0), Nc = NaN; else, Nc = N(i0); end
end

% ============================================================================
function rep(name, N, i0, TRI)
    if isnan(i0)
        fprintf('  收斂點（%s）：未達判準\n', name);
    else
        fprintf('  收斂點（%s）：N = %d，設計 (%d,%d,%d)\n', name, N(i0), TRI(i0,:));
    end
end

% ============================================================================
function mk(tag, N, v1, v2, Nc1, Nc2, ylab, zero_anchor, figdir, YTfix, XSCALE)
% 一張圖：橫軸 N、single 藍 / eighteen 紅、收斂點以同色虛線標出。
%   YTfix（選填）：明給 y 刻度（等距）；上下界由「外側刻度再留一個間距」推得。
%   不給則走 axlim_auto 自動取刻度。
%   XSCALE（選填）：'log'（預設，主刻度標一般數字 10 / 100）
%                  | 'logexp'（對數軸、主刻度標 **10^n 指數形式**）
%                  | 'linear'（線性軸；N=0 錨點直接畫在 0，內部刻度取奇數個等距）
%     三者的端點數字都仍以 text 補、首末點貼框。
    if nargin < 10, YTfix  = []; end
    if nargin < 11 || isempty(XSCALE), XSCALE = 'log'; end
    isLin = strcmpi(XSCALE, 'linear');
    isExp = strcmpi(XSCALE, 'logexp');
    FS = 36;   LW = 3.0;   MS = 10;
    c1 = [0.05 0.10 0.95];   c2 = [0.85 0.10 0.10];
    fig = figure('Color','w','Position',[100 100 1120 820]);
    ax  = axes(fig);   hold(ax,'on');

    % log 軸畫不了 N=0 → 錨點擺在最左端（往左一個十年），端點數字標 0
    % 線性軸沒這個問題，0 就畫在 0
    Np = N;
    if zero_anchor && any(N == 0) && ~isLin
        pos = N(N > 0);   Np(N == 0) = min(pos)/10;
    end

    h1 = plot(ax, Np, v1, '-o', 'Color',c1, 'LineWidth',LW, 'MarkerSize',MS, ...
              'MarkerFaceColor',c1, 'Clipping','off');
    h2 = plot(ax, Np, v2, '-s', 'Color',c2, 'LineWidth',LW, 'MarkerSize',MS, ...
              'MarkerFaceColor','w', 'Clipping','off');

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.02 .02]);
    ax.Toolbar.Visible = 'off';

    XL = [min(Np) max(Np)];   xlim(ax, XL);
    lb0 = sprintf('%g',XL(1));   lb1 = sprintf('%g',XL(2));
    if zero_anchor && any(N == 0), lb0 = '0'; end
    if isLin
    % [ADDED 2026-08-14] 線性橫軸：內部刻度取**奇數個等距**（不含端點，figure-style），
    %   端點數字另以 text 補、首末資料點貼齊左右框邊；線性軸不放次刻度。
    set(ax, 'XScale', 'linear', 'XMinorTick', 'off');
    XT = ticks_lin_inner(XL(1), XL(2), 3);
    set(ax, 'XTick', XT);
    set(ax, 'XTickLabel', arrayfun(@(v) sprintf('%g', v), XT, 'UniformOutput', false));
    else
    set(ax,'XScale','log');
    % [MODIFIED 2026-08-13] **標準對數軸**（使用者要求「對數刻度畫得像樣」）：
    %   主刻度 = 整十年（10 / 100 / 1000），標一般數字（不用 10^3 指數形式）；
    %   次刻度 = 每個十年內的 2..9，短線無數字 —— 科學圖表的慣例畫法。
    %   端點（起點/終點）仍照 figure-style 用 text 補數字 ⇒ 會與太靠近的整十年刻度撞字，
    %   故依「字元數 × 每字元軸寬」估算所需間距，撞到的整十年刻度**拿掉**
    %   （例：6~1776 的軸，1000 距右端僅 10% 軸寬、放不下「1000」+「1776」→ 拿掉 1000；
    %    10 距左端 9%，但「10」+「6」只需 ~5%，故保留）。次刻度不受影響，軸仍是標準對數樣貌。
    e0  = floor(log10(XL(1)));   e1 = ceil(log10(XL(2)));
    dec = 10.^(e0:e1);           dec = dec(dec >= XL(1) & dec <= XL(2));
    sp  = log10(XL(2)/XL(1));    CW  = 0.024;      % 每字元約佔軸寬比例（FS36 粗體）
    keep = true(size(dec));
    for q = 1:numel(dec)
        % isExp 時標籤是 10^{n}，顯示寬度不隨指數變（上標較小）→ 一律以 3 個字元估
        if isExp, nd = 3; else, nd = numel(sprintf('%g', dec(q))); end
        keep(q) = log10(dec(q)/XL(1))/sp >= 0.5*CW*(nd+numel(lb0)) + 0.010 && ...
                  log10(XL(2)/dec(q))/sp >= 0.5*CW*(nd+numel(lb1)) + 0.010;
    end
    XT = dec(keep);
    set(ax, 'XTick', XT);
    if isExp
        % [ADDED 2026-08-14] 主刻度以 10^n 指數形式呈現（使用者要求）。用 tex interpreter
        %   做上標即可 —— 刻度數字仍是 Helvetica 粗體，符合 figure-style（不套 latex）。
        set(ax, 'TickLabelInterpreter', 'tex');
        set(ax, 'XTickLabel', arrayfun(@(v) sprintf('10^{%d}', round(log10(v))), ...
                                       XT, 'UniformOutput', false));
    else
        set(ax, 'XTickLabel', arrayfun(@(v) sprintf('%g', v), XT, 'UniformOutput', false));
    end
    mv = [];
    for e = e0:e1, mv = [mv, (2:9)*10^e]; end %#ok<AGROW>
    set(ax, 'XMinorTick', 'on');
    ax.XAxis.MinorTickValues = mv(mv >= XL(1) & mv <= XL(2));
    end

    if isempty(YTfix)
        [YL, YT] = axlim_auto(min([v1 v2]), max([v1 v2]), [3 5]);
    else
        YT = YTfix;   s = YT(2) - YT(1);   YL = [YT(1)-s, YT(end)+s];
        assert(min([v1 v2]) > YL(1) && max([v1 v2]) < YL(2), '資料超出指定刻度範圍');
    end
    ylim(ax, YL);   set(ax,'YTick',YT);

    % 收斂點虛線（同色；兩模型同值時只畫一條深灰）
    if ~isnan(Nc1) && ~isnan(Nc2) && Nc1 == Nc2
        xline(ax, Nc1, '--', 'Color',[0.25 0.25 0.25], 'LineWidth',2.5, 'HandleVisibility','off');
    else
        if ~isnan(Nc1), xline(ax, Nc1, '--', 'Color',c1, 'LineWidth',2.5, 'HandleVisibility','off'); end
        if ~isnan(Nc2), xline(ax, Nc2, '--', 'Color',c2, 'LineWidth',2.5, 'HandleVisibility','off'); end
    end

    % 水平軸起訖：只標數字、不畫 tick mark（起點若是 N=0 錨點則標 0）
    yoff = YL(1) - 0.022*diff(YL);
    for q = 1:2
        lb = {lb0, lb1};
        text(ax, XL(q), yoff, lb{q}, 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    xlabel(ax, '$\mathbf{Number\;of\;points}$', 'Interpreter','latex', 'FontSize',FS);
    ylabel(ax, ylab, 'Interpreter','latex', 'FontSize',FS);

    lg = legend(ax, [h1 h2], {'Single parameter', 'Eighteen parameters'}, ...
                'Interpreter','tex', 'Location','northoutside', 'NumColumns',2);
    lg.FontSize = 24;  lg.FontWeight = 'bold';  lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;
    hold(ax,'off');

    drawnow;
    axp = get(ax,'Position');   lgp = get(lg,'Position');
    lgw = lgp(3);   lgh = lgp(4);
    GAPN = 0.022;   newTop = 1 - lgh - GAPN - 0.006;
    axp(4) = newTop - axp(2);   set(ax,'Position',axp);
    if lgw < 0.70*axp(3)
        set(lg, 'Position', [axp(1) + (axp(3)-lgw)/2, newTop + GAPN, lgw, lgh]);
    else
        set(lg, 'Position', [axp(1), newTop + GAPN, axp(3), lgh]);
    end

    out = fullfile(figdir, sprintf('%s_vs_npts_conv_maxwell_R150.png', tag));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function tk = ticks_lin_inner(lo, hi, n)
% 線性軸：n 個（奇數）等距 nice 刻度，嚴格落在 (lo,hi) **內部**、不含端點。
%   端點數字由呼叫端以 text 補（figure-style：曲線首末點貼框、端點只標數字不畫 tick mark）。
%   n 是**目標**個數：若該範圍湊不出 n 個 nice 刻度，依序退到 n+2 / n-2 / n+4（**都維持奇數**）。
%   例：(0,30) 湊不出 3 個 nice 刻度（10/20 只有 2 個），退到 5 個 → [5 10 15 20 25]。
    cand = [1 2 2.5 3 4 5 10];
    for nn = [n, n+2, n-2, n+4]
        if nn < 1, continue; end
        k0 = floor(log10(max(hi-lo, realmin)/(nn+1)));
        for k = k0:(k0+3)
            for c = cand
                s = c*10^k;
                t = s*ceil((lo + 0.5*s)/s) : s : (hi - 0.5*s);
                if numel(t) == nn, tk = t;  return; end
            end
        end
    end
    tk = lo + (1:n)/(n+1)*(hi-lo);        % 退路：等分（極端範圍才會用到）
end

% ============================================================================
function [lim, tk] = axlim_auto(lo, hi, nlist)
% 奇數個等距 tick、兩端留白 = tick 間距（house style：各腳本自帶一份）。
    cand = [1 2 2.5 3 4 5 10];
    mid  = (lo+hi)/2;   rng_ = max(hi-lo, realmin);
    best = {};   bestSpan = inf;
    for n = nlist
        k0 = floor(log10(rng_/(n+1)));
        for k = k0:(k0+4)
            hit = false;
            for c = cand
                s   = c*10^k;
                ctr = round(mid/s)*s;
                T = { ctr + (-(n-1)/2 : (n-1)/2)*s };
                if lo >= 0 && lo < s, T = [{(0:n-1)*s}, T]; end %#ok<AGROW>
                for it = 1:numel(T)
                    t   = T{it};
                    L_  = [t(1)-s, t(end)+s];
                    clr = 0.15*s;
                    if lo >= L_(1)+clr && hi <= L_(2)-clr
                        if (n+1)*s < bestSpan, bestSpan = (n+1)*s;  best = {L_, t}; end
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
