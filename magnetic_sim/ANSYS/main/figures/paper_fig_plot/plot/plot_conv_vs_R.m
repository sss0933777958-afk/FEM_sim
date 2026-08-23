function plot_conv_vs_R(force)
% plot_conv_vs_R -- 「用收斂點的點數校正，拿真實 FEM 格點驗證」隨取樣半徑 R 的變化
% =========================================================================
%   對每一個取樣半徑 R 做三件事：
%     ① **沿用 plot_full_vs_conv_vs_R.m 已找到的該 R 收斂設計**（讀其快取
%        full_vs_conv_vs_R_maxwell.mat 的 tri_c1/tri_c2）。這樣兩支腳本的 N_c
%        保證一致，也省掉重跑一次階梯搜尋（那是最貴的一段）。
%        判準＝三者交集：l_hat 穩定 ∧ g_I 穩定 ∧ K_I 符合物理（各持續 KWIN=10 步）。
%     ② 用該設計的減量取樣點校正 → l_hat、e、G、g_I_hat；
%     ③ **拿該 R 內全部真實 .fld 格點當評估集**（不是拿來校正的內插查詢點），
%        用校正出來的 (l_hat, e, G) 去預測，算殘差。
%
%   [MODIFIED 2026-08-13] 原本此腳本自己跑階梯、判準只有「l_hat + K_I 號誌」且 KWIN=7，
%   與現行判準不一致（少了 g_I 條件）→ 改成讀 full_vs_conv 的設計。
%
%   四張圖（橫軸都是 R [um]，single 藍 / eighteen 紅）：
%     ell   : l_hat [um]                                （收斂點校正值）
%     gain  : g_I_hat [mT/A]                            （同上）
%     rms   : sqrt(J/N) [mT]  —— J = 殘差平方和、N = 殘差項數（點數 x 3 分量 x 6 激發）
%     nmae  : mean|res| / mean|B|  [%]                  —— 兩者皆對全部殘差項取平均
%
%   ⚠ 評估時**用校正階段解出來的 G**（不在評估集上重解），所以這是真正的外推驗證：
%     幾何參數與電荷強度都被檢驗，而不是只檢驗幾何。
%
%   風格①粗體框圖：FS 36 粗體、box on、grid off、線性橫軸、刻度奇數等距、
%   曲線首末點貼齊左右框邊、起訖數字以 text 補、圖例照 figure-style 標準樣式。
%
%   輸出 → figures/paper_fig/Section2_E/{ell,gain,rms,nmae}_vs_R_conv_maxwell.png
% =========================================================================
    clc;
    if nargin < 1 || isempty(force), force = false; end

    l0    = 0.5e-3;                  % l_hat 初值 [m]
    MODEL = 'long2016_hexapole_halfcut';   GEOM = 'tip40um';
    % R 掃描範圍與逐 R 的收斂設計都來自 full_vs_conv_vs_R_maxwell.mat（見檔頭）

    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    cachef = fullfile(here, 'data', 'conv_vs_R_maxwell.mat');

    %% ---- 計算（快取）------------------------------------------------------
    if exist(cachef,'file') && ~force
        S = load(cachef);   fprintf('由快取載入 %s\n', cachef);
    else
        S = sweep_R(MODEL, GEOM, l0, here);
        save(cachef, '-struct', 'S');   fprintf('已存 %s\n', cachef);
    end

    fprintf('\n  R [um]   : %s\n', num2str(S.R_um,   '%7d'));
    fprintf('  N_c 1p   : %s\n',   num2str(S.Nc1,    '%7d'));
    fprintf('  N_c 18p  : %s\n',   num2str(S.Nc2,    '%7d'));
    fprintf('  l_hat 1p : %s um\n',  num2str(S.ell1, '%7.1f'));
    fprintf('  l_hat 18p: %s um\n',  num2str(S.ell2, '%7.1f'));
    fprintf('  gain  1p : %s mT/A\n',num2str(S.gI1,  '%7.3f'));
    fprintf('  gain  18p: %s mT/A\n',num2str(S.gI2,  '%7.3f'));
    fprintf('  minJ  1p : %s mT\n',  num2str(S.rmsc1,'%7.4f'));   % sqrt(minJ/(3·Nc·6))
    fprintf('  minJ  18p: %s mT\n',  num2str(S.rmsc2,'%7.4f'));
    fprintf('  grid  1p : %s mT\n',  num2str(S.rms1, '%7.3f'));   % 外推到全格點（未出圖）
    fprintf('  grid  18p: %s mT\n',  num2str(S.rms2, '%7.3f'));
    fprintf('  NMAE  1p : %s %%\n',  num2str(S.nmae1,'%7.2f'));
    fprintf('  NMAE  18p: %s %%\n',  num2str(S.nmae2,'%7.2f'));
    fprintf('  評估格點 : %s\n',     num2str(S.Neval,'%7d'));

    %% ---- 四張圖 -----------------------------------------------------------
    % [MODIFIED 2026-08-14] 曲線自 R=0 起（使用者指定），起始值 = 0。
    % ⚠ l_hat 與 g_I 的 R-sweep **只保留 plot_ell_gain_2panel.m 那一版**（左右兩格、
    %   同一張圖）—— 本檔不再產 ell_vs_R_conv / gain_vs_R_conv（重複呈現，使用者拍板）。
    % [MODIFIED 2026-08-13] rms 圖改畫**校正當下的 min J**（不是外推到全格點的殘差）：
    %   sqrt(min J / (3·N_c·6))。min J 來自 fitting 的 variable projection 閉式解。
    %   分母用**純量項數**（3 分量 × N_c 點 × 6 激發），這樣單位才是 mT、才是每個殘差
    %   分量的 RMS。外推到全格點的那組仍算在 S.rms1/S.rms2（存快取，未出圖）。
    mk('rms',  S.rmsc1, S.rmsc2, '$\mathbf{\sqrt{J/N}\;(mT)}$',                  figdir, S.R_um, 0);
    mk('nmae', S.nmae1, S.nmae2, '$\mathbf{NMAE\;(\%)}$',                        figdir, S.R_um, 0);
end

% ============================================================================
function S = sweep_R(MODEL, GEOM, l0, here)
% 逐 R：取 full_vs_conv 已定的收斂設計 → 用該設計的減量點校正 → 在真實格點上評估。
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    addpath(fullfile(CAL,'function'), fullfile(CAL,'common_path'), fullfile(CAL,'utils'));
    addpath(fullfile(CAL,'utils','long2016_hexapole_halfcut'));

    cfg = model_config(MODEL, GEOM);
    F   = zeros(6, cfg.N_I);
    for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end

    % 評估集來源：全部 .fld 格點（一次載入，之後按 R 取子集）
    raw = extract_maxwell_data(cfg, 'all', cfg.default_variant);
    ad  = build_actuator_data(raw, cfg);

    % 逐 R 的收斂設計：直接讀 plot_full_vs_conv_vs_R 的快取（判準與 N_c 由那支決定）
    df = fullfile(here, 'data', 'full_vs_conv_vs_R_maxwell.mat');
    assert(exist(df,'file')==2, '找不到 %s —— 先跑 plot_full_vs_conv_vs_R(true)', df);
    D    = load(df);
    R_um = D.R_um(:).';   nR = numel(R_um);
    TRI  = {D.tri_c1, D.tri_c2};

    [S.Nc1,S.Nc2,S.ell1,S.ell2,S.gI1,S.gI2,S.rms1,S.rms2,S.nmae1,S.nmae2,S.Neval, ...
     S.rmsc1,S.rmsc2] = deal(nan(1,nR));                % rmsc = sqrt(min J / (3·N_c·6))
    S.tri1 = zeros(nR,3);   S.tri2 = zeros(nR,3);

    for a = 1:nR
        R = R_um(a)*1e-6;
        [P0, B0, n0] = cfg.select_ball(ad, R);          % 評估集：該 R 內全部真實格點
        S.Neval(a) = n0;

        for m = 1:2
            tri = TRI{m}(a,:);
            if any(tri == 0), continue; end             % 該 R 在 full_vs_conv 未達判準
            % [MODIFIED 2026-08-23] sphere_grid_sample 已併入 conv_design_ws。
            [P, Bs] = conv_design_ws(tri(1), tri(2), tri(3), R, ...
                                     struct('frame','actuator'));
            np = size(P,1);
            r = fit_pack(P, Bs, cfg.Pc_base, l0, m==2, F);
            [rms, nmae] = eval_on_grid(P0, B0, r, cfg.Pc_base);
            rmsc = sqrt(r.J / numel(Bs));               % numel(Bs) = 3·N_c·6（純量項數）
            if m == 1
                S.Nc1(a)=np; S.ell1(a)=r.l*1e6; S.gI1(a)=r.gI;
                S.rms1(a)=rms; S.nmae1(a)=nmae; S.tri1(a,:)=tri;  S.rmsc1(a)=rmsc;
            else
                S.Nc2(a)=np; S.ell2(a)=r.l*1e6; S.gI2(a)=r.gI;
                S.rms2(a)=rms; S.nmae2(a)=nmae; S.tri2(a,:)=tri;  S.rmsc2(a)=rmsc;
            end
        end
        fprintf(['  R=%3d um｜評估格點 %6d｜1p: N=%5d (%d,%d,%2d) l=%6.1f g=%6.3f ' ...
                 'rms=%.4f nmae=%.2f%%｜18p: N=%5d (%d,%d,%2d) l=%6.1f g=%6.3f rms=%.4f nmae=%.2f%%\n'], ...
                R_um(a), n0, S.Nc1(a), S.tri1(a,1),S.tri1(a,2),S.tri1(a,3), S.ell1(a), S.gI1(a), S.rms1(a), S.nmae1(a), ...
                S.Nc2(a), S.tri2(a,1),S.tri2(a,2),S.tri2(a,3), S.ell2(a), S.gI2(a), S.rms2(a), S.nmae2(a));
    end
    S.R_um = R_um;
end

% ============================================================================
function r = fit_pack(P, Bstack, Pc_base, l0, USE_BIAS, F)
% 擬合一個設計，回傳打包好的校正結果（l, e, G, gI, K, J）。
%   J = fitting 收斂後的 **min J**（variable projection：電荷強度以閉式解
%   gj = (SᵀS)\(SᵀB_j) 投影掉，只在 l_hat[+e] 上搜尋）。
%   長度＝numel(Bstack)=3·N_c·6，故 J 已含 3 分量 × 6 激發。
    try
        [e, l, J] = fitting(P, Bstack, Pc_base, l0, USE_BIAS);
        [K, gI, G] = solve_current(l, e, Pc_base, P, Bstack, F);
        r = struct('l',l, 'e',e, 'G',G, 'gI',gI, 'K',K, 'J',J);
    catch ME
        warning('fit_pack:degenerate','N=%d 擬合失敗（%s）', size(P,1), ME.identifier);
        r = struct('l',NaN, 'e',[], 'G',[], 'gI',NaN, 'K',[], 'J',NaN);
    end
end

% ============================================================================
function [rms, nmae] = eval_on_grid(P0, B0, r, Pc_base)
% 拿校正出來的 (l, e, G) 去預測**真實 FEM 格點**的場，算殘差。
%   rms  = sqrt(J/N)              J = 殘差平方和、N = 殘差項數（點 x 3 分量 x 6 激發）[mT]
%   nmae = 向量範數版（[MODIFIED 2026-08-18] 與 solve_current/solve_voltage 對齊）：
%          Σ_j Σ_i ‖res_ij‖ / Σ_j Σ_i ‖B_ij‖ · 100，先逐「點×激發」取 3 維殘差長度再相加。
%          舊的分量版 `mean|res|/mean|B|` 已作廢（規則 calibration-transfer-matrix-output 附則三）。
    if ~isfinite(r.l), rms = NaN;  nmae = NaN;  return; end
    Pc  = make_Pc(r.e, Pc_base);
    S   = build_S(r.l, Pc, P0);
    res = S*r.G - B0;
    rms  = sqrt(sum(res(:).^2) / numel(res));
    e_ij = sqrt(sum(reshape(res, 3, []).^2, 1));       % 每點每激發的 ‖res‖
    b_ij = sqrt(sum(reshape(B0,  3, []).^2, 1));
    nmae = 100 * sum(e_ij) / sum(b_ij);
end

% ============================================================================
function Pc = make_Pc(e17, Pc_base)
    E = zeros(3, 6);
    if ~isempty(e17)
        E(:,1) = e17(1:3);   E(:,2) = e17(4:6);
        E(:,3) = e17(7:9);   E(:,4) = e17(10:12);
        E(:,5) = e17(13:15);
        E(1,6) = e17(16);    E(2,6) = e17(17);
        E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);
    end
    Pc = Pc_base + E;
end

% ============================================================================
function S = build_S(l_hat, Pc, P)
    Np   = size(P, 1);
    pbar = P / l_hat;
    S    = zeros(3*Np, 6);
    for k = 1:6
        d  = pbar - Pc(:,k).';
        r3 = sum(d.^2, 2).^1.5;
        S(:,k) = reshape((d ./ r3).', 3*Np, 1);
    end
end

% ============================================================================
function mk(tag, v1, v2, ylab, figdir, R_um, v0)
% 一張圖：橫軸 R [um]（線性），single 藍 / eighteen 紅。
%   v0（選填）：**R=0 的起始值**（使用者指定曲線自 0 起）。給了就在最前面補一個
%     (0, v0) 錨點，橫軸自 0 起；ell 用擬合初值 500、其餘量用 0。
    if nargin < 7, v0 = []; end
    if ~isempty(v0)
        R_um = [0, R_um(:).'];   v1 = [v0, v1(:).'];   v2 = [v0, v2(:).'];
    end
    FS = 36;   LW = 3.0;   MS = 10;
    c1 = [0.05 0.10 0.95];   c2 = [0.85 0.10 0.10];
    fig = figure('Color','w','Position',[100 100 1120 820]);
    ax  = axes(fig);   hold(ax,'on');

    h1 = plot(ax, R_um, v1, '-o', 'Color',c1, 'LineWidth',LW, 'MarkerSize',MS, ...
              'MarkerFaceColor',c1, 'Clipping','off');
    h2 = plot(ax, R_um, v2, '-s', 'Color',c2, 'LineWidth',LW, 'MarkerSize',MS, ...
              'MarkerFaceColor','w', 'Clipping','off');

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.02 .02]);
    ax.Toolbar.Visible = 'off';

    XL = [min(R_um) max(R_um)];   xlim(ax, XL);
    if isempty(v0)
        set(ax,'XTick',[200 300 400]);                   % 3 個等距內部刻度，端點另以 text 標
    else
        set(ax,'XTick',100:100:400);                     % 自 0 起：整數刻度，端點 0/500 以 text 標
    end
    [YL, YT] = axlim_auto(min([v1 v2]), max([v1 v2]), [3 5]);
    ylim(ax, YL);   set(ax,'YTick',YT);

    yoff = YL(1) - 0.022*diff(YL);
    for xv = XL
        text(ax, xv, yoff, sprintf('%g', xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    xlabel(ax, '$\mathbf{R\;(micro\;meter)}$', 'Interpreter','latex', 'FontSize',FS);
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

    out = fullfile(figdir, sprintf('%s_vs_R_conv_maxwell.png', tag));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
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
