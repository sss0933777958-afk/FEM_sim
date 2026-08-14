function plot_full_vs_conv_vs_R(force)
% plot_full_vs_conv_vs_R -- 「全格點校正」vs「收斂點固定設計校正」隨取樣半徑 R 的對照
% =========================================================================
%   兩張圖，每張上下兩個子圖（上 l_hat、下 g_I_hat），橫軸都是取樣半徑 R [um]：
%
%     圖 1  single   (1-param)  ：全格點  vs  **該 R 自己的收斂點設計**
%     圖 2  eighteen (18-param) ：全格點  vs  **該 R 自己的收斂點設計**
%
%   減量設計 = **每個 R 各自**沿「按比例輪流 +1」階梯搜尋，取三判準的交集
%   （各自持續 5 步）：
%     (a) l_hat   相對前一個設計的變化 < 0.5%（自足判準）
%     (b) g_I_hat 相對前一個設計的變化 < 0.5%（同一把尺）
%     (c) K_I_bar 符合物理：off-diag 全負 + 對角全正 + 對角占優
%   [MODIFIED 2026-08-13] (c) 原為「K_I_bar 的 Frobenius 相對變化 < 0.5%」，
%   已證實對 single 不可達（見 conv_fit 內註解），改為物理有效性閘門；
%   同時補上 (b)，否則 g_I 無人約束（eighteen 會惡化到 3.14%）。
%
%   四個子圖各算一個 **RMSPE**（減量 vs 全格點）：
%       RMSPE = sqrt( sum (v_conv - v_full)^2 / sum v_full^2 ) * 100  [%]
%   數值寫進圖例（照 figure-style「統計數值寫進圖例文字」）。
%
%   資料來源
%     全格點：既有快取 data/ell_gain_sweep_maxwell{,_bias}.mat（plot_ell_gain_vs_R 產）
%     減量  ：sphere_grid_sample 逐 R 重算（點數固定、只有球半徑變）
%
%   R = 60:20:500，起點與 *_vs_R_conv_maxwell.png 那組統一
%   （R<=50 時六顆電荷在小球內幾乎不可分辨，g_I 會翻負，實測 R=50 給 -8.1 mT/A）。
%
%   風格①粗體框圖：FS 36 粗體、box on、grid off、tiledlayout(2,1) 共用橫軸
%   （上子圖不印 x 刻度數字）、線性橫軸、曲線首末點貼齊框邊、起訖數字以 text 補。
%   配色：全格點 = 黑實心圓、減量 = 該模型色（single 藍 / eighteen 紅）空心方。
%
%   輸出 → figures/paper_fig/Section2_E/full_vs_conv_vs_R_{single,eighteen}_maxwell.png
% =========================================================================
    clc;
    if nargin < 1 || isempty(force), force = false; end

    R_um  = 20:20:500;               % 取樣半徑 [um]
                                     % [MODIFIED 2026-08-14] 使用者要求「取樣範圍自小往大、
                                     % 每 20um 跑一次流程」→ 起點降到 **40**（= 全格點快取
                                     % ell_gain_sweep_maxwell*.mat 的最小 R；R=20 那顆快取沒有、
                                     % R=0 球內 0 點無法擬合）。
                                     % ⚠ 已知：R=40/60 小球內六顆電荷近簡併，conv 可能判不到
                                     % 收斂（回 NaN）；全格點基準本身也只有 28 / 106 點，
                                     % 且 R=40 的 single 全格點 g_I = 27.6 mT/A（正常 ~9.5）。
                                     % [歷史] 2026-08-13 曾因此把起點由 60 提到 80。
    l0    = 0.5e-3;
    TOL   = 0.005;                   % 兩判準共用：相對變化 < 0.5%
    % [MODIFIED 2026-08-13] 視窗 5 → 10（使用者拍板；R=150 的實測顯示 K=5..7 三個量的
    %   收斂點完全相同，K=10 是最小會動的值）。連帶把 NDMAX 從 120 拉到 150：判準變嚴、
    %   每個判準都要多 5 步才成立，大 R 端本來就吃緊。
    KWIN  = 10;                      % 三判準共用：後 10 步
    NDMAX = 150;                     % 每個 R 最多走幾個設計（早停）
    MODEL = 'long2016_hexapole_halfcut';   GEOM = 'tip40um';

    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    cachef = fullfile(here, 'data', 'full_vs_conv_vs_R_maxwell.mat');

    %% ---- 計算（快取）------------------------------------------------------
    if exist(cachef,'file') && ~force
        S = load(cachef);
        if isequal(S.R_um(:).', R_um(:).')
            fprintf('由快取載入 %s\n', cachef);
        else
            % [ADDED 2026-08-14] **增量補跑**：R 清單改變時，舊快取已有的 R 直接沿用，
            %   只算新增的那幾個（全部重跑要數小時，九成花在大 R 端的大設計上）。
            fprintf('快取 R 清單與設定不符 → 只補跑新增的 R\n');
            S = sweep(here, MODEL, GEOM, R_um, l0, TOL, KWIN, NDMAX, S);
            save(cachef, '-struct', 'S');   fprintf('已更新 %s\n', cachef);
        end
    else
        S = sweep(here, MODEL, GEOM, R_um, l0, TOL, KWIN, NDMAX);
        save(cachef, '-struct', 'S');   fprintf('已存 %s\n', cachef);
    end

    %% ---- RMSPE（減量 vs 全格點）------------------------------------------
    % 未收斂的 R 會是 NaN（conv_fit 已 warn），RMSPE 只取兩邊皆有限的點
    rp = @(a,b) 100*sqrt(sum((a(isfinite(a)&isfinite(b))-b(isfinite(a)&isfinite(b))).^2) ...
                       / sum(b(isfinite(a)&isfinite(b)).^2));
    % [ADDED 2026-08-14] **本檔的兩張圖與 RMSPE 只取 R>=80**（維持既有論文圖不變）：
    %   快取雖已擴到 R=40（供 plot_ell_gain_2panel 用），但 R=40 的全格點 single
    %   g_I = 27.6 mT/A（正常 ~9.5，小球內六顆電荷近簡併）會把 y 軸拉爆、
    %   也會污染 RMSPE。R=40/60 的 conv-single 本來就是 NaN（判不到收斂）。
    pm   = S.R_um >= 80;
    r_e1 = rp(S.ell_c1(pm), S.ell_f1(pm));   r_g1 = rp(S.gI_c1(pm), S.gI_f1(pm));
    r_e2 = rp(S.ell_c2(pm), S.ell_f2(pm));   r_g2 = rp(S.gI_c2(pm), S.gI_f2(pm));
    fprintf('\n  RMSPE（減量 vs 全格點）\n');
    fprintf('    single   ：l_hat %.3f%%   g_I %.3f%%   (N = %d-%d)\n', r_e1, r_g1, S.N1(1), S.N1(2));
    fprintf('    eighteen ：l_hat %.3f%%   g_I %.3f%%   (N = %d-%d)\n', r_e2, r_g2, S.N2(1), S.N2(2));

    %% ---- 兩張圖 -----------------------------------------------------------
    % g_I 的 y 刻度使用者指定為整數（single 9/10、eighteen 10/11），不走自動取刻度
    mk('single',   S.R_um, S.ell_f1, S.ell_c1, S.gI_f1, S.gI_c1, ...
       [0.05 0.10 0.95], S.N1, r_e1, r_g1, figdir, [9 10]);
    mk('eighteen', S.R_um, S.ell_f2, S.ell_c2, S.gI_f2, S.gI_c2, ...
       [0.85 0.10 0.10], S.N2, r_e2, r_g2, figdir, [10 11]);
end

% ============================================================================
function S = sweep(here, MODEL, GEOM, R_um, l0, TOL, KWIN, NDMAX, old)
% 全格點曲線讀既有快取；減量曲線逐 R 重算（設計固定、只有球半徑變）。
%   old（選填）：舊快取。其 R 已算過者直接沿用，只跑新增的 R（判準相同才可沿用）。
    if nargin < 9, old = []; end
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    addpath(fullfile(CAL,'function'), fullfile(CAL,'common_path'), fullfile(CAL,'utils'));
    addpath(fullfile(CAL,'utils','long2016_hexapole_halfcut'));

    cfg = model_config(MODEL, GEOM);
    F   = zeros(6, cfg.N_I);
    for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end

    % --- 全格點：既有 R-sweep 快取（plot_ell_gain_vs_R 產，R = 40:20:500）---
    A = load(fullfile(here,'data','ell_gain_sweep_maxwell.mat'));       % single
    B = load(fullfile(here,'data','ell_gain_sweep_maxwell_bias.mat'));  % eighteen
    assert(isequal(A.Rum(:), B.Rum(:)), '兩顆全格點快取的 R 不一致');
    % [MODIFIED 2026-08-14] 全格點快取缺某些 R 時**不再中止**，該點的 full 曲線填 NaN。
    %   （R=20 不在全格點快取裡；本檔的兩張圖只畫 R>=80，而 plot_ell_gain_2panel 只用 conv。）
    [tf, loc] = ismember(R_um, A.Rum);
    nR0 = numel(R_um);
    [S.ell_f1, S.gI_f1, S.ell_f2, S.gI_f2, S.npts_f] = deal(nan(1,nR0));
    S.ell_f1(tf) = A.ell_R(loc(tf));   S.gI_f1(tf) = A.gI_R(loc(tf));
    S.ell_f2(tf) = B.ell_R(loc(tf));   S.gI_f2(tf) = B.gI_R(loc(tf));
    S.npts_f(tf) = A.npts_R(loc(tf));
    if ~all(tf)
        fprintf('⚠ 全格點快取缺 R = %s um → 該點 full 曲線填 NaN（conv 不受影響）\n', ...
                num2str(R_um(~tf)));
    end

    % --- 減量：**每個 R 各自沿階梯搜尋收斂點**（雙判準交集）---
    w = [1, 3, 3*pi];   tt = [1 2 3];   TRI = zeros(NDMAX,3);
    for q = 1:NDMAX, TRI(q,:) = tt;  [~, jj] = min(tt ./ w);  tt(jj) = tt(jj) + 1;  end
    nR = numel(R_um);
    [S.ell_c1, S.gI_c1, S.ell_c2, S.gI_c2, S.n_c1, S.n_c2] = deal(nan(1,nR));
    S.tri_c1 = zeros(nR,3);   S.tri_c2 = zeros(nR,3);
    S.fb_c1  = false(1,nR);   S.fb_c2  = false(1,nR);   % [ADDED] true = 只用 l_hat 收斂點的 fallback
    % [ADDED 2026-08-14] 舊快取可沿用的條件：判準（TOL/KWIN）一致
    reuse = ~isempty(old) && isfield(old,'R_um') && ...
            isfield(old,'TOL') && isequal(old.TOL,TOL) && ...
            isfield(old,'KWIN') && isequal(old.KWIN,KWIN);
    for a = 1:nR
        R = R_um(a)*1e-6;
        % [MODIFIED 2026-08-14] 兩個模型**分開**沿用；舊值是 NaN 的重算（現在有 fallback 了）
        got1 = false;   got2 = false;
        if reuse
            j = find(old.R_um == R_um(a), 1);
            if ~isempty(j)
                if ~isnan(old.ell_c1(j))
                    S.ell_c1(a) = old.ell_c1(j);  S.gI_c1(a) = old.gI_c1(j);
                    S.n_c1(a)   = old.n_c1(j);    S.tri_c1(a,:) = old.tri_c1(j,:);
                    if isfield(old,'fb_c1'), S.fb_c1(a) = old.fb_c1(j); end
                    got1 = true;
                end
                if ~isnan(old.ell_c2(j))
                    S.ell_c2(a) = old.ell_c2(j);  S.gI_c2(a) = old.gI_c2(j);
                    S.n_c2(a)   = old.n_c2(j);    S.tri_c2(a,:) = old.tri_c2(j,:);
                    if isfield(old,'fb_c2'), S.fb_c2(a) = old.fb_c2(j); end
                    got2 = true;
                end
            end
        end
        if got1 && got2
            fprintf('  R=%3d um | 沿用舊快取（conv 1p N=%d / 18p N=%d）\n', R_um(a), S.n_c1(a), S.n_c2(a));
            continue;
        end
        if ~got1
            [S.ell_c1(a), S.gI_c1(a), S.n_c1(a), S.tri_c1(a,:), S.fb_c1(a)] = ...
                conv_fit(R, TRI, cfg, l0, false, F, TOL, KWIN);
        end
        if ~got2
            [S.ell_c2(a), S.gI_c2(a), S.n_c2(a), S.tri_c2(a,:), S.fb_c2(a)] = ...
                conv_fit(R, TRI, cfg, l0, true,  F, TOL, KWIN);
        end
        fprintf('  R=%3d um | full %6d : 1p l=%6.1f g=%6.3f | 18p l=%6.1f g=%6.3f\n', ...
                R_um(a), S.npts_f(a), S.ell_f1(a), S.gI_f1(a), S.ell_f2(a), S.gI_f2(a));
        fprintf('           conv 1p N=%5d (%d,%d,%2d) l=%6.1f g=%6.3f | 18p N=%5d (%d,%d,%2d) l=%6.1f g=%6.3f\n', ...
                S.n_c1(a), S.tri_c1(a,1),S.tri_c1(a,2),S.tri_c1(a,3), S.ell_c1(a), S.gI_c1(a), ...
                S.n_c2(a), S.tri_c2(a,1),S.tri_c2(a,2),S.tri_c2(a,3), S.ell_c2(a), S.gI_c2(a));
    end
    S.R_um = R_um;   S.TRI = TRI;   S.TOL = TOL;   S.KWIN = KWIN;
    S.N1 = [min(S.n_c1) max(S.n_c1)];   S.N2 = [min(S.n_c2) max(S.n_c2)];
    % [ADDED] 全部欄位強制成**橫列**。快取取出來的 ell_R/gI_R 是直行、這裡算的是橫列，
    %   不統一的話 (a-b) 會被廣播成 nR x nR 矩陣：RMSPE 變成向量、fprintf 跟著循環輸出，
    %   最後 [ellF ellC] 直接 horzcat 失敗（踩過）。
    for f = {'ell_f1','gI_f1','ell_f2','gI_f2','npts_f','ell_c1','gI_c1','ell_c2','gI_c2', ...
             'n_c1','n_c2','R_um'}
        S.(f{1}) = S.(f{1})(:).';
    end
end

% ============================================================================
function [ell_um, gI, n, tri, isfb] = conv_fit(R, TRI, cfg, l0, USE_BIAS, F, TOL, KWIN)
% 沿階梯走，找**三判準的交集**（各自都要持續 KWIN 步）：
%   (a) l_hat 穩定    ：「相對前一個設計」的變化 < TOL
%   (b) g_I_hat 穩定  ：同一把尺（相對前一個設計的變化 < TOL）
%   (c) K_I_bar 符合物理：off-diag 全負 + 對角全正 + 對角占優
% 回該設計的 l_hat / g_I / 實得點數 / 三元組。三判準都成立即早停。
%
% ⚠ 為什麼要有 (b)（2026-08-13 加）：判準必須約束**每一個要報的量**。只用 (a)+(c) 時，
%   eighteen 的 K_I_bar 在階梯第一個設計（N=6）就號誌全對 → (c) 形同虛設、判準退化成只剩
%   (a)，g_I 完全沒人管：實測 g_I 對全格點的 RMSPE 由 0.40% 惡化到 3.14%（R=180~260 的
%   g_I 飄到 11.3 mT/A，全格點只有 10.5~10.8）。
%
% ⚠ 為什麼 K_I_bar 不用數值容差（2026-08-13 實測後改）：
%   以 R=150 全 1771 格點的 K_I_bar 為參考，等測度網格的 K_I_bar 在
%   N=204 差 10.7%、N=1200 差 5.0%、N=1776 差 4.1%，誤差只以 N^-0.41 下降，
%   外推到 0.5% 需約 22 萬點 —— 對 single 而言數值收斂不可達。成因是模型設定不良
%   （殘差有結構 → 最小平方解隨取樣測度改變）；佐證：eighteen 用同一套內插取樣器、
%   同樣 88 點，g_I 對全格點只差 0.01%。故 K_I_bar 改判「號誌結構是否物理」。
% ⚠ l_hat 刻意用**自足**判準（比前一個設計）而非比全格點：否則圖上的 RMSPE
%   會變成自我實現，失去獨立驗證的意義。
    nD = size(TRI,1);
    [ell, np, gIv] = deal(nan(1,nD));
    ok = false(1,nD);   i0 = NaN;
    for q = 1:nD
        try
            [x,y,z,B] = sphere_grid_sample(R, [], struct('frame','actuator','NRPT',TRI(q,:)));
            P  = [x y z];   np(q) = size(P,1);
            Bs = zeros(3*np(q), size(B,3));
            for j = 1:size(B,3), Bs(:,j) = reshape(B(:,:,j).', [], 1); end
            [e, l] = fitting(P, Bs, cfg.Pc_base, l0, USE_BIAS);
            [K, gIv(q)] = solve_current(l, e, cfg.Pc_base, P, Bs, F);
            ell(q) = l*1e6;
            ok(q)  = ki_physical(K);
        catch
            ell(q) = NaN;   gIv(q) = NaN;   ok(q) = false;
        end
        ie = first_stable(ell(1:q), TOL, KWIN);
        ig = first_stable(gIv(1:q), TOL, KWIN);
        ik = first_true(ok(1:q), KWIN);
        if ~isnan(ie) && ~isnan(ig) && ~isnan(ik), i0 = max([ie ig ik]);  break; end
    end
    % [ADDED 2026-08-14] fallback（使用者拍板）：三判準交集達不到時，**只要 l_hat 有收斂
    %   就用 l_hat 的收斂點**，g_I 取同一個設計的值（不是另外找 g_I 的收斂點——它根本沒有）。
    %   適用 R<=60：實測 single 的 g_I 與 K_I 號誌在 55 個設計內都判不到，但 l_hat 都收斂。
    %   ⚠ 這種點回報的 g_I **未經自身收斂驗證**，isfb=true 供上層標示 / 排除。
    isfb = false;
    if isnan(i0)
        ie = first_stable(ell, TOL, KWIN);
        if isnan(ie)
            warning('conv_fit:noconv','R=%.0f um (USE_BIAS=%d) %d 個設計內 l_hat 也未收斂', ...
                    R*1e6, USE_BIAS, nD);
            ell_um = NaN;  gI = NaN;  n = NaN;  tri = [0 0 0];  return
        end
        warning('conv_fit:ellonly', ...
                'R=%.0f um (USE_BIAS=%d) 三判準交集未達 → 改用 l_hat 收斂點 N=%d（g_I 未經收斂驗證）', ...
                R*1e6, USE_BIAS, np(ie));
        i0 = ie;   isfb = true;
    end
    ell_um = ell(i0);   gI = gIv(i0);   n = np(i0);   tri = TRI(i0,:);
end

% ============================================================================
function tf = ki_physical(K)
% K_I_bar 的物理有效性：自激發（對角）為正、交叉項（off-diag）全負、且對角占優。
    od = K(~eye(6));
    [~, am] = max(abs(K), [], 2);
    tf = all(od < 0) && all(diag(K) > 0) && isequal(am(:).', 1:6);
end

% ============================================================================
function i0 = first_true(v, K)
% 第一個 i，使 v(i..i+K-1) 全為 true。
    i0 = NaN;
    for i = 1:numel(v)-K+1
        if all(v(i:i+K-1)), i0 = i;  return; end
    end
end

% ============================================================================
function i0 = first_stable(v, tol, K)
% 第一個 i，使其後連續 K 個「相對前一點」的變化率都 < tol。
    i0 = NaN;
    if numel(v) < K+1, return; end
    rel = abs(diff(v) ./ v(1:end-1));
    for i = 1:numel(rel)-K+1
        if all(rel(i:i+K-1) < tol), i0 = i;  return; end
    end
end

% ============================================================================
function mk(tag, R_um, ellF, ellC, gF, gC, cM, ~, rmse_e, rmse_g, figdir, gYT)
% 一張圖、上下兩子圖（上 l_hat、下 g_I），共用橫軸 R。
%   兩個子圖**各有自己的圖例**（照 figure-style 標準樣式），RMSPE 寫進該子圖的圖例文字，
%   不用浮動 text。座標軸手動定位（不走 tiledlayout）——圖例要各自貼齊自己的框。
    FS = 36;   LW = 3.0;   MS = 9;   GAPN = 0.045;   % GAPN 要夠大：圖例外框不可壓到座標框
    cF = [0.15 0.15 0.15];                       % 全格點 = 黑
    fig = figure('Color','w','Position',[100 40 1180 1360]);

    XL = [min(R_um) max(R_um)];   XT = [200 300 400];
    AXX = 0.160;   AXW = 0.795;   AXH = 0.350;   % 兩軸共用左緣/寬/高

    % ---- 上：l_hat ----
    ax1 = axes(fig, 'Position', [AXX 0.560 AXW AXH]);   hold(ax1,'on');
    h1 = plot(ax1, R_um, ellF, '-o', 'Color',cF, 'LineWidth',LW, 'MarkerSize',MS, ...
              'MarkerFaceColor',cF, 'Clipping','off');
    h2 = plot(ax1, R_um, ellC, '-s', 'Color',cM, 'LineWidth',LW, 'MarkerSize',MS, ...
              'MarkerFaceColor','w', 'Clipping','off');
    style_ax(ax1, FS, XL, XT);   ax1.XTickLabel = {};          % 共用軸：上子圖不印 x 數字
    [YL1, YT1] = axlim_auto(min([ellF ellC],[],'omitnan'), max([ellF ellC],[],'omitnan'), [3 5]);
    ylim(ax1, YL1);   set(ax1,'YTick',YT1);
    ylabel(ax1, '$\mathbf{\hat{\ell}\;(micro\;meter)}$', 'Interpreter','latex', 'FontSize',FS);
    lg1 = mklg(ax1, [h1 h2], rmse_e);
    hold(ax1,'off');

    % ---- 下：g_I ----
    ax2 = axes(fig, 'Position', [AXX 0.090 AXW AXH]);   hold(ax2,'on');
    h3 = plot(ax2, R_um, gF, '-o', 'Color',cF, 'LineWidth',LW, 'MarkerSize',MS, ...
              'MarkerFaceColor',cF, 'Clipping','off');
    h4 = plot(ax2, R_um, gC, '-s', 'Color',cM, 'LineWidth',LW, 'MarkerSize',MS, ...
              'MarkerFaceColor','w', 'Clipping','off');
    style_ax(ax2, FS, XL, XT);
    % [MODIFIED 2026-08-13] g_I 用使用者指定的**整數**刻度（single 9/10、eighteen 10/11）。
    %   g_I 的資料跨距只有 ~0.5-0.7 mT/A：single 區間內只含一個整數、eighteen 一個都不含，
    %   所以刻度由外部給定，上下界改用「外側刻度再留 0.15」而非自動取刻度。
    YT2 = gYT;   YL2 = [YT2(1)-0.15, YT2(end)+0.15];
    assert(min([gF gC],[],'omitnan') > YL2(1) && max([gF gC],[],'omitnan') < YL2(2), ...
           'g_I 資料超出指定刻度範圍');
    ylim(ax2, YL2);   set(ax2,'YTick',YT2);
    ylabel(ax2, '$\mathbf{{}^{B}\hat{g}_{I}\;(mT/A)}$', 'Interpreter','latex', 'FontSize',FS);
    xlabel(ax2, '$\mathbf{Sampling\;range\;(micro\;meter)}$', 'Interpreter','latex', 'FontSize',FS);
    lg2 = mklg(ax2, [h3 h4], rmse_g);
    % 水平軸起訖：只標數字、不畫 tick mark
    yoff = YL2(1) - 0.030*diff(YL2);
    for xv = XL
        text(ax2, xv, yoff, sprintf('%g', xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    hold(ax2,'off');

    place_lg(lg1, ax1, GAPN);   place_lg(lg2, ax2, GAPN);

    out = fullfile(figdir, sprintf('full_vs_conv_vs_R_%s_maxwell.png', tag));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function lg = mklg(ax, h, rmspe)
% figure-style 標準圖例：northoutside、黑粗外框、FS24 粗體、兩欄（系列｜含統計值的系列）。
    lg = legend(ax, h, {'All grid points', ...
                sprintf('Reduce sampling (RMSPE = %.2f%%)', rmspe)}, ...
                'Interpreter','tex', 'Location','northoutside', 'NumColumns',2);
    lg.FontSize = 24;   lg.FontWeight = 'bold';
    lg.Box = 'on';      lg.EdgeColor = 'k';   lg.LineWidth = 2.5;
end

% ============================================================================
function place_lg(lg, ax, gap)
% 置於該座標框正上方；寬度自動貼合內容（<70% 框寬則置中、否則左右緣切齊框）。
    drawnow;
    axp = get(ax,'Position');   lgp = get(lg,'Position');
    lgw = lgp(3);   lgh = lgp(4);
    if lgw < 0.70*axp(3)
        set(lg, 'Position', [axp(1)+(axp(3)-lgw)/2, axp(2)+axp(4)+gap, lgw, lgh]);
    else
        set(lg, 'Position', [axp(1), axp(2)+axp(4)+gap, axp(3), lgh]);
    end
end

% ============================================================================
function style_ax(ax, FS, XL, XT)
    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.02 .02]);
    xlim(ax, XL);   set(ax,'XTick',XT);
    ax.Toolbar.Visible = 'off';
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
