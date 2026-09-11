function plot_conv_vs_R(force, MODEL, GEOM, VARIANT, SUF, WHICH)
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
    % [MODIFIED 2026-08-27] 參數化 model/geom/variant + 快取後綴，改讀 *_fullref 快取
    %   （新取樣器 sample_equal_h + 新判準「對全格點 3% 連續 10 級」）。
    if nargin < 2 || isempty(MODEL),   MODEL   = 'long2016_hexapole_halfcut'; end
    if nargin < 3,                     GEOM    = 'tip40um';                   end
    if nargin < 4,                     VARIANT = '';                          end
    if nargin < 5 || isempty(SUF),     SUF     = '_fullref';                  end
    % [ADDED 2026-08-27] WHICH：要產哪幾張圖（預設兩張都產，維持原行為）。
    %   {'rms'} 只產 sqrt(J/N)、{'nmae'} 只產 NMAE。
    if nargin < 6 || isempty(WHICH), WHICH = {'rms','nmae'}; end
    if ischar(WHICH), WHICH = {WHICH}; end
    % R 掃描範圍與逐 R 的收斂設計都來自 full_vs_conv_vs_R_maxwell.mat（見檔頭）

    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    cachef = fullfile(here, 'data', ['conv_vs_R_maxwell' SUF '.mat']);

    %% ---- 計算（快取）------------------------------------------------------
    if exist(cachef,'file') && ~force
        S = load(cachef);   fprintf('由快取載入 %s\n', cachef);
    else
        S = sweep_R(MODEL, GEOM, VARIANT, SUF, l0, here);
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
    % [MODIFIED 2026-08-27 使用者拍板] rms 圖改畫**全格點外推**殘差（S.rms1/rms2）：
    %   用該 R 收斂點的校正參數去預測 R 內全部真實 .fld 格點，J = 殘差平方和。
    %   in-sample 的 sqrt(minJ/(3·N_c·6))（S.rmsc1/rmsc2）仍算、仍存快取，但不出圖
    %   —— 它會因為 eighteen 自由度多而天生偏低，不能用來比「預測能力」。
    if any(strcmp(WHICH,'rms'))
        mk('rms',  S.rms1,  S.rms2,  '$\mathbf{\sqrt{J/N}\;(mT)}$',              figdir, S.R_um, 0, SUF);
    end
    if any(strcmp(WHICH,'nmae'))
        mk('nmae', S.nmae1, S.nmae2, '$\mathbf{NMAE\;(\%)}$',                    figdir, S.R_um, 0, SUF);
    end
end

% ============================================================================
function S = sweep_R(MODEL, GEOM, VARIANT, SUF, l0, here)
% 逐 R：取 full_vs_conv 已定的收斂設計 → 用該設計的減量點校正 → 在真實格點上評估。
    % [MODIFIED 2026-08-27] Tree moved to matlab\Flux\; the hardcoded matlab\Maxwell no
    %   longer exists. Locate relative to this file and add temp_code\scripts, which holds
    %   the new samplers (sample_equal_h / sample_rings / sample_axes4).
    MAIN = fileparts(fileparts(here));                 % ...\ANSYS\main
    CAL  = fullfile(MAIN, 'matlab', 'Flux', 'Maxwell');
    addpath(fullfile(CAL,'function'), fullfile(CAL,'common_path'), fullfile(CAL,'utils'));
    addpath(fullfile(MAIN,'temp_code','scripts'));

    cfg = model_config(MODEL, GEOM);
    if isempty(VARIANT), VARIANT = cfg.default_variant; end
    F   = zeros(6, cfg.N_I);
    for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end

    % 評估集來源：全部 .fld 格點（一次載入，之後按 R 取子集）
    raw = extract_maxwell_data(cfg, 'all', VARIANT);
    ad  = build_actuator_data(raw, cfg);

    % 逐 R 的收斂設計：直接讀 plot_full_vs_conv_vs_R 的快取（判準與 N_c 由那支決定）
    df = fullfile(here, 'data', ['full_vs_conv_vs_R_maxwell' SUF '.mat']);
    assert(exist(df,'file')==2, '找不到 %s —— 先跑 plot_full_vs_conv_vs_R(true)', df);
    D    = load(df);
    R_um = D.R_um(:).';   nR = numel(R_um);
    % [MODIFIED 2026-08-27] tri_c* in the new cache are all zero: the equal-h ladder is a
    %   single index, not an (N_r,N_phi,N_theta) triple. A design is now identified by its
    %   POINT COUNT n_c*, rebuilt from the KIND/SPEC ladder. fb_c* = fell back to full grid.
    % [ADDED 2026-08-28] Two cache formats now exist:
    %   '_fullref' : the retired mixed ehgap ladder -> carries KIND/SPEC + fb_c*, design
    %                rebuilt from them by point count (design_by_npts).
    %   '_axsh'    : sample_axes_shells (six actuator axes on Nr equally spaced shells +
    %                centre, N = 6*Nr+1) -> the point count alone identifies the design,
    %                so N -> Nr is inverted directly. That sweep has no full-grid fallback.
    AXSH = ~isfield(D,'KIND');
    KIND = [];   SPEC = [];
    if ~AXSH, KIND = D.KIND;  SPEC = D.SPEC; end
    NC   = {D.n_c1(:).',  D.n_c2(:).'};
    if isfield(D,'fb_c1')
        FB = {logical(D.fb_c1(:).'), logical(D.fb_c2(:).')};
    else
        % [ADDED 2026-08-28] The axsh cache has no fb_c*; it marks a full-grid fallback
        %   by src=="full grid" (zhi_peng R=500: the axis points fall inside iron, so the
        %   ladder never converges).  Detect it from n_c == npts_f, which holds either way
        %   -- inverting 64769 into an axes-shells Nr would be nonsense.
        FB = {NC{1} == D.npts_f(:).', NC{2} == D.npts_f(:).'};
    end
    o    = struct('model',MODEL, 'geom',GEOM, 'variant',VARIANT, ...
                  'frame','actuator', 'quiet',true);

    [S.Nc1,S.Nc2,S.ell1,S.ell2,S.gI1,S.gI2,S.rms1,S.rms2,S.nmae1,S.nmae2,S.Neval, ...
     S.rmsc1,S.rmsc2] = deal(nan(1,nR));                % rmsc = sqrt(min J / (3·N_c·6))
    S.tri1 = zeros(nR,3);   S.tri2 = zeros(nR,3);
    S.fb1  = FB{1};         S.fb2  = FB{2};    % fell back to full grid -> in-sample

    for a = 1:nR
        R = R_um(a)*1e-6;
        [P0, B0, n0] = cfg.select_ball(ad, R);          % 評估集：該 R 內全部真實格點
        S.Neval(a) = n0;

        for m = 1:2
            nc = NC{m}(a);
            if ~isfinite(nc), continue; end             % full_vs_conv 跳過了這個 R
            if FB{m}(a)
                P = P0;  Bs = B0;                       % 退回全格點：校正集 = 評估集
            elseif AXSH
                % [MODIFIED 2026-08-28] n_c in the cache is the count AFTER iron filtering,
                %   so it is not always 6*Nr+1: at R=500 the outermost shell sits exactly on
                %   the pole tips (R_norm = 500 um) and the points along the pole axes are
                %   dropped (811 -> 808, 3889 -> 3886).  Recover Nr by trying the smallest Nr
                %   whose nominal count reaches nc and checking the kept count.
                P = [];  Bs = [];
                for Nr_ = ceil((nc-1)/6) + (0:2)
                    Pq_ = conv_design_ws(Nr_, R, struct('points_only',true,'R_act',cfg.R_act,'quiet',true));
                    evalc('[Pk_, Bk_] = conv_design_ws([], R, setfield(o,''query'',Pq_));');
                    if size(Pk_,1) == nc,  P = Pk_;  Bs = Bk_;  break;  end
                end
                if isempty(P)
                    fprintf(['  R=%3d um m=%d: no axes-shells Nr gives N=%d after iron' ...
                             ' filtering, skipped' newline], R_um(a), m, nc);
                    continue;
                end
            else
                [P, Bs] = design_by_npts(KIND, SPEC, nc, R, o, cfg);
                if isempty(P)
                    fprintf(['  R=%3d um m=%d: no ladder step with N=%d,' ...
                             ' skipped\n'], R_um(a), m, nc);
                    continue;
                end
            end
            tri = [0 0 0];
            np  = size(P,1);
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
function mk(tag, v1, v2, ylab, figdir, R_um, v0, SUF)
% 一張圖：橫軸 R [um]（線性），single 藍 / eighteen 紅。
%   v0（選填）：**R=0 的起始值**（使用者指定曲線自 0 起）。給了就在最前面補一個
%     (0, v0) 錨點，橫軸自 0 起；ell 用擬合初值 500、其餘量用 0。
    if nargin < 7, v0 = []; end
    if nargin < 8 || isempty(SUF), SUF = ''; end
    if ~isempty(v0)
        R_um = [0, R_um(:).'];   v1 = [v0, v1(:).'];   v2 = [v0, v2(:).'];
    end
    % [MODIFIED 2026-09-02 使用者拍板的新繪圖規則，見 .claude/rules/figure-style.md]
    %   規則1 刻度數字 = 60（原本 FS 同時餵刻度與軸標題，這裡拆成 FS / FSLAB）；
    %   規則2 圖例 = 45 且圖例線段樣本與資料線一併加粗；規則3 框線加粗；規則6 畫布等邊。
    %   軸標題 FSLAB 使用者未指定 -> 沿用本腳本原值 36，不自行編造。
    %   LW / MS 依刻度 36->60（x1.67）同量級放大，與 axsh_gain / axsh_kfro 取同一組值。
    FS = 60;   FSLAB = 36;   FSLEG = 45;   LW = 5.0;   MS = 12;   LWBOX = 5.0;
    c1 = [0.05 0.10 0.95];   c2 = [0.85 0.10 0.10];
    % 畫布等邊且以英吋給（不用像素：超過螢幕時 MATLAB 會靜默縮小視窗、比例跑掉）。
    CANV = 14.5;                                   % 畫布邊長 [in]，與兩支 axsh 同值
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes(fig);   hold(ax,'on');

    % [MODIFIED 2026-08-27] R=20 is NaN (its full grid holds only 2 points) and that NaN
    %   broke the line between the R=0 anchor and R=40. Drop non-finite entries per
    %   series so the curve is continuous; the surviving points keep their x positions.
    R_um = R_um(:).';   v1 = v1(:).';   v2 = v2(:).';
    k1 = isfinite(v1);   k2 = isfinite(v2);
    h1 = plot(ax, R_um(k1), v1(k1), '-o', 'Color',c1, 'LineWidth',LW, 'MarkerSize',MS, ...
              'MarkerFaceColor',c1, 'Clipping','off');
    h2 = plot(ax, R_um(k2), v2(k2), '-s', 'Color',c2, 'LineWidth',LW, 'MarkerSize',MS, ...
              'MarkerFaceColor','w', 'Clipping','off');

    box(ax,'on');  grid(ax,'off');
    % [MODIFIED 2026-09-02 使用者指定] TickDir='out'：刻度線朝外。
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
           'TickLength',[.02 .02],'TickDir','out');
    ax.Toolbar.Visible = 'off';

    XL = [min(R_um) max(R_um)];   xlim(ax, XL);
    if isempty(v0)
        set(ax,'XTick',[200 300 400]);                   % 3 個等距內部刻度，端點另以 text 標
    else
        set(ax,'XTick',100:100:400);                     % 自 0 起：整數刻度，端點 0/500 以 text 標
    end
    % [MODIFIED 2026-09-02 規則 4] 縱軸**起點與終點都不標**。原本用 axlim_auto，當 lo=0 時
    %   它會產生從 0 起算的刻度組 -> 0 被標出來（實測兩張圖的 y 都標了 0）。改用與
    %   axsh_kfro 相同的 ylim_from_zero：把 [0,T] 平分 n+1 段、刻度 = (1:n)*s，首末不進 YTick。
    [YL, YT] = ylim_from_zero(max([v1 v2]));
    ylim(ax, YL);   set(ax,'YTick',YT);
    fprintf('  縱軸 ylim=[%g %g]  YTick=%s  (%d 根，填充 %.0f%%)\n', ...
            YL(1), YL(2), mat2str(YT), numel(YT), 100*max([v1 v2])/YL(2));

    yoff = YL(1) - 0.022*diff(YL);
    for xv = XL
        text(ax, xv, yoff, sprintf('%g', xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    xlabel(ax, '$\mathbf{R\;(micro\;meter)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    ylabel(ax, ylab, 'Interpreter','latex', 'FontSize',FSLAB);

    lg = legend(ax, [h1 h2], {'Single parameter', 'Eighteen parameters'}, ...
                'Interpreter','tex', 'Location','northoutside', 'NumColumns',2);
    lg.FontSize = FSLEG;  lg.FontWeight = 'bold';
    lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;
    % [FIXED 2026-09-02] ItemTokenSize（圖例 icon 區寬）預設固定 30 pt、不隨 FontSize 放大，
    %   字級 45 pt 下 marker 會貼到圖例框線上。與 axsh_gain / axsh_kfro 取同值。
    lg.ItemTokenSize = [55 25];
    hold(ax,'off');

    drawnow;
    axp = get(ax,'Position');   lgp = get(lg,'Position');
    lgw = lgp(3);   lgh = lgp(4);
    GAPN = 0.022;   newTop = 1 - lgh - GAPN - 0.006;
    axp(4) = newTop - axp(2);   set(ax,'Position',axp);
    if lgw < 0.70*axp(3)
        set(lg, 'Position', [axp(1) + (axp(3)-lgw)/2, newTop + GAPN, lgw, lgh]);
    else
        % [FIXED 2026-09-02] 給圖例自然寬度並夾住右緣。原本把寬度設成 axp(3)，但 MATLAB
        %   不會把 legend 縮到比內容窄，只會對稱置中 -> 圖例比座標框寬時兩側一起溢出畫布。
        lgx = max(0.015, min(axp(1), 1 - 0.015 - lgw));
        set(lg, 'Position', [lgx, newTop + GAPN, lgw, lgh]);
    end

    % [MODIFIED 2026-08-28] '_axsh' is stripped like '_fullref': the sampler is which
    %   ladder produced N_c, not a different figure -- the new run overwrites the old one.
    out = fullfile(figdir, sprintf('%s_vs_R_conv_maxwell%s.png', ...
                                   tag, strrep(strrep(SUF,'_fullref',''),'_axsh','')));
    % [MODIFIED 2026-09-02 規則 6] print 取代 exportgraphics：後者會裁掉畫布四周白邊，
    %   正方形畫布匯出後就不是正方形（實測長寬比 0.96 / 1.11）。
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
    print(fig, out, '-dpng', '-r200');             % 14.5 in x 200 dpi = 2900 px 見方
    fprintf('wrote %s\n', out);
end

% ============================================================================
function [yl, tk] = ylim_from_zero(maxv)
% [REPLACED 2026-09-02] 取代原本的 axlim_auto —— 那支在 lo=0 時會產生「從 0 起算」的
%   刻度組，於是縱軸起點 0 被標出來，違反規則 4（縱軸起點與終點都不標）。
%   本函式與 axsh_kfro.m 的同名 local 完全相同（兩處刻意保持一致）：
%   自 0 起，把 [0, T] 平分 n+1 段，刻度 = (1:n)*s，起點與終點都不進 tk。
%   [MODIFIED 2026-09-02 使用者指定] n 由 3 改 **4**（與 axsh_gain 的 gain 那組一致）。
%   ⚠ 4 是偶數、與規則 5「奇數個 tick」相衝突 —— 這是使用者對本組圖的明確指示，
%     跟水平軸的整數刻度例外並列，日後規則正式定案時要一起寫進去。
%   s 取 0.1 的倍數，在 [smin, 1.15*smin] 窗口內挑最漂亮者
%   （整數 > 0.5 的倍數 > 0.2 的倍數 > 其餘，同分取最小），確保刻度是整數或一位小數。
    n = 4;
    smin = 1.08*maxv/(n+1);
    k0 = max(1, ceil(smin/0.1 - 1e-9));
    k1 = max(k0, ceil(1.15*smin/0.1));
    best = k0;   bs = -1;
    for k = k0:k1
        if     mod(k,10) == 0, sc = 3;      % 整數
        elseif mod(k, 5) == 0, sc = 2;      % 0.5 的倍數
        elseif mod(k, 2) == 0, sc = 1;      % 0.2 的倍數
        else,                  sc = 0;
        end
        if sc > bs, bs = sc;   best = k; end
    end
    pick = best*0.1;
    tk = round((1:n)*pick*10)/10;   yl = [0, (n+1)*pick];
end

% ============================================================================
function [P, Bs] = design_by_npts(KIND, SPEC, nc, R, o, cfg)
% [ADDED 2026-08-27] Recover a converged design from the mixed ladder by POINT COUNT.
%   The new criterion sweeps a ladder mixing four samplers (ax / rg / gf / eh); a design
%   is no longer an (N_r,N_phi,N_theta) triple, and the cache keeps only the count, so we
%   rebuild each rung and take the first whose point count matches.
    P = [];  Bs = [];
    for q = 1:numel(KIND)
        switch KIND{q}
            case 'eh',  Pq = sample_equal_h(R, SPEC{q}, struct('quiet',true));
            case 'rg',  Pq = sample_rings(  R, SPEC{q}, struct('quiet',true));
            case 'ax',  Pq = sample_axes4(R, cfg.R_act, ...
                                 struct('quiet',true,'naxes',SPEC{q},'plus',true));
            % [REMOVED 2026-09-02] case 'gf'（等測度球格）已隨等測度取樣法一併廢除。
            otherwise,  continue;
        end
        if isempty(Pq) || size(Pq,1) ~= nc, continue; end
        try,  evalc('[Pk,Bk] = conv_design_ws(Pq, R, o);');  catch, continue;  end
        if size(Pk,1) == nc,  P = Pk;  Bs = Bk;  return;  end
    end
end
