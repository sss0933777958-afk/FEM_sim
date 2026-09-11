function plot_sigma_hist(USE_BIAS, R_FIT, R_EVAL, NSAMP, force, MODEL, GEOM, VARIANT, XR, BINW, NTK)
% plot_sigma_hist -- 三個奇異值 σ₁ σ₂ σ₃ 在整個工作空間的分布疊圖（Long Fei 半切六極）
% =========================================================================
%   使用者拍板 2026-08-18：
%     ① **校正**用 R <= R_FIT（150 µm）—— 與 Section4_C 既有極座標圖同一把尺：
%        取該 R 的**收斂設計 N_c**（等測度網格降取樣）擬合 → (ℓ̂, e, ᴮĤ_I)。
%     ② **評估**用 R <= R_EVAL（500 µm）**整個工作空間的真實 .fld 格點**：
%        逐點組 S(p)·ᴮĤ_I（3×6）做一次 SVD → 三個奇異值 σ₁ ≥ σ₂ ≥ σ₃ [mT/A]。
%        ⚠ 每個位置**一次** SVD（六個激發是矩陣的六個欄，不是六個樣本）→
%          N_p 個點就是 N_p 組 (σ₁,σ₂,σ₃)。
%     ③ 三組 σ 疊在同一張直方圖（各自正規化成百分比）；點夠多（R=500 有 65 353 點）
%        分布看起來是連續的。
%
%   物理：σ₁/σ₂/σ₃ 是「單位電流向量在該點能產生的場」的三個主軸增益（mT/A）。
%   σ₁ = 最強方向、σ₃ = 最弱方向；κ = σ₃/σ₁ 是等向性、𝒞 = ∏σ 是致動體積。
%
%   風格：①粗體框圖 + 疊圖直方圖慣例（nb=180、三組**共用 edges**、FaceAlpha、
%   百分比縱軸自 0 起留 8%、x 端點只標數字不畫 tick）；圖例左欄 σ_k、右欄 mean
%   （照 figure-style「左欄 = 資料系列、右欄 = 該系列的統計值」，**不畫 mean 虛線**）。
%
%   快取 data/sigma_hist_R<eval>_maxwell_<tag>.mat（重讀 .fld 要數十秒）。
%   輸出 → figures/paper_fig/Section4_C/sigma_hist_R<eval>_maxwell[_eighteen].png
% =========================================================================
    clc;
    if nargin < 1 || isempty(USE_BIAS), USE_BIAS = false; end   % false = single(ℓ̂ only)
    if nargin < 2 || isempty(R_FIT),    R_FIT    = 150;   end   % 校正半徑 [µm]
    if nargin < 3 || isempty(R_EVAL),   R_EVAL   = 500;   end   % 評估半徑 [µm]（整個工作空間）
    % [ADDED 2026-08-18 使用者拍板] NSAMP = 評估點數。
    %   ⚠ SVD 不需要 FEM 資料（S(p)·ᴮĤ_I 是**解析式**），所以評估位置可以任意加密 ——
    %     這正是極座標圖看起來連續的原因（它在 161×181 的極座標網格上連續評估）。
    %   [MODIFIED 2026-08-21 使用者拍板] **預設改回真實格點**（NSAMP = []）。
    %   NSAMP = []（預設）→ 用 R<=R_EVAL 的**真實 .fld 格點**（65 353 點）
    %   NSAMP = 數字      → 球內**均勻亂數取樣** NSAMP 點（等體積；曾用 2e6 讓每 bin
    %                       樣本數 ×30、雜訊由 ~3% 降到 ~0.5%，分布更平滑）
    if nargin < 4,                      NSAMP    = [];    end
    if nargin < 5 || isempty(force),    force    = false; end
    % [ADDED 2026-08-21] model 參數化（原本寫死 long2016）→ 志鵬平面六極也能出同一張圖。
    %   非 long2016 時檔名（圖 + 快取）自動加 _<model>，不覆蓋既有的龍飛版。
    %   ⚠ VARIANT 必須明給：同一個 model 有多版場時（zhi_peng: maxwell / maxwell_split）
    %     不傳會退回 cfg.default_variant，**靜默**用到舊場。
    if nargin < 6 || isempty(MODEL),   MODEL   = 'long2016_hexapole_halfcut'; end
    if nargin < 7,                     GEOM    = 'tip40um';  end
    if nargin < 8,                     VARIANT = '';         end
    if strcmp(MODEL,'long2016_hexapole_halfcut') && isempty(GEOM), GEOM = 'tip40um'; end
    % [ADDED 2026-08-21 使用者拍板] 水平軸範圍。σ₁ 的尾巴極長（志鵬 max 281、龍飛 max 58），
    %   舊做法 xr = [floor(min), ceil(max)] 會讓右邊 20% 幾乎全空、且端點是 282 / 59 這種
    %   不乾淨的數（違反 figure-style「水平軸端點必須是乾淨的數」）。改成明給終點：
    %     龍飛 → 55（截掉 0.0x%）   志鵬 → 250（截掉 0.07%）
    %   被截掉的比例會印在 console。傳 XR 可覆寫；傳 [0 0] 可退回舊的貼齊極值行為。
    if nargin < 9 || isempty(XR)
        switch MODEL
            case 'long2016_hexapole_halfcut', XR = [4 55];
            case 'zhi_peng',                  XR = [5 250];
            otherwise,                        XR = [];      % 未知 model → 貼齊極值
        end
    end
    if isequal(XR, [0 0]), XR = []; end
    % [ADDED 2026-08-21] BINW = 直接指定 bin 寬 [mT/A]；[] = 沿用 nb=180（鋪滿各自值域）。
    %   用途：讓不同 model / 不同 R_EVAL 的圖**共用同一把尺**，長條高度才可互比
    %   （面積 = 100%% x bin 寬）。例：long2016 R150 的 bin 寬 0.022042，志鵬 R150
    %   傳同一個值，兩張圖就能直接比高度。
    %   ⚠ 代價：值域寬的那張 bin 數暴增（志鵬 33.08/0.022 = 1501 根），每根樣本數
    %     掉到個位數 -> 統計雜訊大、形狀變毛。
    % [MODIFIED 2026-09-02 使用者要求] long2016 的 R_EVAL=150 預設改為 **明給 0.02**，
    %   讓本圖與 gain_cbrt_hist_A_single_R150 **共用同一把尺**。兩張都是 mT/A、
    %   都在 R=150 評估，bin 寬不同時長條高度不可互比（面積 = 100%% x bin 寬）：
    %     舊值 nb=180 -> 3.99113/180 = 0.022173；gain 圖 fd_bin(nbar=30) -> 0.020000（差 10.9%%）。
    %   改成 0.02 後本圖 180 -> 200 根（每根約 27 點），gain 圖不變。
    %   ⚠ 要改這個值必須同步改 plot_gain_iso_hist.m 的 BW_SOLO_G{1}，否則兩張又跨掉。
    %   其他 model / 其他 R_EVAL 維持 []（志鵬值域 33 mT/A，套 0.02 會變 1650 根）。
    % [MODIFIED 2026-09-02 使用者指定] **同一個模型的圖要用相同 bin 寬**。
    %   本圖（sigma_kk）與該模型的 gain_cbrt_hist（C^(1/3)）都是 mT/A，故逐模型釘成
    %   與 plot_gain_iso_hist 的 BW_SOLO_G / fd_bin 結果相同的值：
    %     long2016      -> 0.02   (gain A    fd_bin(nbar=30) 也是 0.02)
    %     zhi_peng split-> 0.1    (gain B    fd_bin 給 0.1)
    %     zhi_peng gap  -> 0.05   (gain Bgap fd_bin 給 0.05)
    %   ⚠ sigma 的值域比 C^(1/3) 寬很多（gap：23.6 vs 1.69）-> 根數會暴增
    %     （180 -> 333 / 472）、每根點數降到 16 / 11，長條會更毛。這是對齊 bin 寬的代價。
    %   ⚠ kappa 圖不在此列：它是**無因次**的 sigma3/sigma1，mT/A 的 bin 寬套不上去
    %     （0.1 套到值域 0.21 只剩 2 根）；kappa 三張另行共用 0.02。
    if nargin < 10
        vt = regexprep(VARIANT, '^maxwell_?', '');
        BINW = [];
        if R_EVAL == 150
            switch MODEL
                case 'long2016_hexapole_halfcut', BINW = 0.02;
                case 'zhi_peng'
                    switch vt
                        case 'split', BINW = 0.1;
                        case 'gap',   BINW = 0.05;
                    end
            end
        end
    end

    % [ADDED 2026-09-02 使用者指定] NTK = 兩軸各要幾根內部刻度（水平與縱軸同值）。
    %   原本 axis_odd / ylim_odd 都寫死 3。使用者對不同圖指定不同根數
    %   （zhi_peng split = 3、zhi_peng gap = 5），故參數化。預設 3 → 已定案的
    %   long2016 版重跑輸出不變。
    if nargin < 11 || isempty(NTK), NTK = 3; end

    tag  = 'single';   if USE_BIAS, tag = 'eighteen'; end
    msfx = '';  if ~strcmp(MODEL,'long2016_hexapole_halfcut'), msfx = ['_' MODEL]; end
    % [ADDED 2026-08-28] variant 後綴：同一 model 有多版場（zhi_peng 的 maxwell_split /
    %   maxwell_gap）時，只帶 model 名會互相覆蓋。照 short-names 剝掉整棵樹都一樣的
    %   'maxwell' -> '_split' / '_gap'（與 plot_svd_polar 同慣例）。
    CAL0 = fullfile('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main','matlab','Flux','Maxwell');
    addpath(fullfile(CAL0,'function'), fullfile(CAL0,'utils'), fullfile(CAL0,'common_path'));
    cfg0 = model_config(MODEL, GEOM);
    if ~isempty(VARIANT) && ~strcmpi(VARIANT, cfg0.default_variant)
        msfx = [msfx '_' regexprep(VARIANT,'^maxwell_?','')];
    end

    here   = fileparts(fileparts(mfilename('fullpath')));       % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section4_C');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    if isempty(NSAMP), nstr = 'grid'; else, nstr = sprintf('N%dk', round(NSAMP/1e3)); end
    cf = fullfile(here, 'data', sprintf('sigma_hist_R%d_maxwell%s_%s_%s.mat', R_EVAL, msfx, tag, nstr));
    if exist(cf,'file') && ~force
        S = load(cf);   fprintf('由快取載入 %s\n', cf);
    else
        S = compute(MODEL, GEOM, VARIANT, R_FIT, R_EVAL, USE_BIAS, here, NSAMP);
        save(cf, '-struct', 'S');   fprintf('已存 %s\n', cf);
    end

    fprintf('校正 R<=%d µm（N_c=%d，設計 %d,%d,%d）｜ℓ̂=%.1f µm  ĝ_I=%.4f mT/A  NMAE=%.2f%%\n', ...
            S.R_fit, S.Nc, S.tri, S.l_hat*1e6, S.gI, S.NMAE);
    fprintf('評估 R<=%d µm：%d 個點（%s）→ %d 組 (σ1,σ2,σ3)\n', S.R_eval, S.npts, S.smode, S.npts);
    for k = 1:3
        v = S.sig(:,k);
        fprintf('  σ%d [mT/A]：min %7.3f  mean %7.3f  max %7.3f  CV %5.2f%%\n', ...
                k, min(v), mean(v), max(v), std(v)/mean(v)*100);
    end
    fprintf('  σ1/σ3 平均比 = %.2f｜σ1 與 σ2 分布重疊率 %.1f%%\n', ...
            mean(S.sig(:,1))/mean(S.sig(:,3)), ...
            100*mean(S.sig(:,2) > min(S.sig(:,1))));

    render(S, figdir, tag, msfx, XR, BINW, NTK);
end

% ============================================================================
function S = compute(MODEL, GEOM, VARIANT, R_FIT, R_EVAL, USE_BIAS, here, NSAMP)
    % [MODIFIED 2026-08-28] Tree moved under matlab/Flux/; the old hardcoded path is gone,
    %   so addpath silently failed and an APDL copy of model_config/solve_* could shadow
    %   the Maxwell ones.
    MAIN = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main';
    CAL  = fullfile(MAIN,'matlab','Flux','Maxwell');
    addpath(fullfile(CAL,'function'), fullfile(CAL,'utils'), fullfile(CAL,'common_path'));

    cfg = model_config(MODEL, GEOM);
    if isempty(VARIANT), VARIANT = cfg.default_variant; end
    raw = extract_maxwell_data(cfg, 'all', VARIANT);
    ad  = build_actuator_data(raw, cfg);

    F = zeros(6, cfg.N_I);
    for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end

    % ---- ① 校正：R<=R_FIT 的收斂設計 N_c（與 plot_svd_polar 同一把尺）----
    %   ki_gate / ki_req 只對 long2016 開：六極不等強的設計（zhi_peng）K̄_I 的物理結構
    %   條件不適用 → 收斂點只由 ℓ̂ 與 ĝ_I 決定（使用者拍板 2026-08-15）。
    is_l2016 = strcmp(MODEL,'long2016_hexapole_halfcut');
    % ---- 讀收斂設計並就地校正（axes-shells 取樣器）------------------------
    %   [MODIFIED 2026-08-28] main.m 的 calib_*_convN*.mat 出自**舊的等測度取樣器**，與
    %   現行判準不一致（sample_axes_shells：六根致動軸 x Nr 層等距球殼 + 中心、N=6*Nr+1；
    %   穩態 = 連續 20 級步進 < 0.01%，收斂 = 連續 10 級落在穩態值 +/-0.2%）。
    %   改成讀 axsh 掃描的 N_c、在此重建設計就地擬合 —— 與 plot_svd_polar 同一把尺。
    %   快取優先序：在 R_FIT 本身跑的 '..._axsh[_zhi]_R<R_FIT>.mat' 勝過 R-sweep
    %   （sweep 格是 40:20:500，R_FIT=150 不在格上，向鄰居借 N_c 是任意的）。
    zsuf_ = '';   if strcmp(MODEL,'zhi_peng'), zsuf_ = '_zhi'; end
    %   [ADDED 2026-08-28] 再往前一層：**variant 專屬**的 R_FIT 快取（split / gap 的
    %   ĝ_I 差 ~34%，收斂點不保證相同，不可共用一顆）。優先序：
    %     ..._axsh<zsuf>_R<R>_<vtag>.mat  >  ..._axsh<zsuf>_R<R>.mat  >  R-sweep
    vtag_ = regexprep(VARIANT,'^maxwell_?','');
    dfv_ = fullfile(here,'data',sprintf('full_vs_conv_vs_R_maxwell_axsh%s_R%d_%s.mat',zsuf_,R_FIT,vtag_));
    dfx_ = fullfile(here,'data',sprintf('full_vs_conv_vs_R_maxwell_axsh%s_R%d.mat',zsuf_,R_FIT));
    df_  = fullfile(here,'data',['full_vs_conv_vs_R_maxwell_axsh' zsuf_ '.mat']);
    if     exist(dfv_,'file')==2, df_ = dfv_;
    elseif exist(dfx_,'file')==2, df_ = dfx_;   end
    assert(exist(df_,'file')==2, 'missing %s', df_);
    D_  = load(df_);
    NCv = D_.n_c1;   if USE_BIAS, NCv = D_.n_c2; end
    FBv = NCv(:).' == D_.npts_f(:).';        % 退回全格點的 R（校正集 = 評估集）要跳過
    Rv_ = D_.R_um(:).';   ok_ = isfinite(NCv(:).') & ~FBv;
    if ~any(ok_), ok_ = isfinite(NCv(:).'); end
    cd_ = find(ok_);   [~,jj_] = min(abs(Rv_(cd_) - R_FIT));   ix_ = cd_(jj_);
    Nc  = NCv(ix_);
    o_  = struct('model',MODEL,'geom',GEOM,'variant',VARIANT,'frame','actuator','quiet',true);
    Pd_ = [];   Bd_ = [];
    for Nr_ = ceil((Nc-1)/6) + (0:2)         % 濾鐵可能吃掉幾點 -> N 不一定剛好是 6*Nr+1
        Pq_ = conv_design_ws(Nr_, R_FIT*1e-6, struct('points_only',true,'R_act',cfg.R_act,'quiet',true));
        evalc('[Pk_,Bk_] = conv_design_ws([], R_FIT*1e-6, setfield(o_,''query'',Pq_));');
        if size(Pk_,1) == Nc,  Pd_ = Pk_;  Bd_ = Bk_;  break;  end
    end
    assert(~isempty(Pd_), 'no axes-shells Nr gives N=%d at R=%g um', Nc, R_FIT);
    [Pf_,Bf_]       = cfg.select_ball(ad, R_FIT*1e-6);   % out-of-sample NMAE 的評估集
    [e, l_hat]      = fitting(Pd_, Bd_, cfg.Pc_base, 0.5e-3, USE_BIAS);
    [KI, gI, ~, rm] = solve_current(l_hat, e, cfg.Pc_base, Pd_, Bd_, F, [], Pf_, Bf_);
    tri = [0 0 0];                           % 階梯是單一索引，不再是 (Nr,Nphi,Ntheta) 三元組
    fprintf(['  校正設計：axes-shells N_c=%d (Nr=%g) @R=%d um，判準取自 R=%d' newline], ...
            Nc, (Nc-1)/6, R_FIT, Rv_(ix_));

    % ---- ② 評估：取樣 → 逐點一次 SVD（批次 pagesvd）----
    if isempty(NSAMP)
        [P, ~, np] = cfg.select_ball(ad, R_EVAL*1e-6);      % 真實 .fld 格點
        smode = 'real .fld grid points';
    else
        % 球內均勻取樣（等體積）：r = R·u^(1/3)、方向在球面上均勻
        rng(0);   np = NSAMP;
        u  = rand(np,1);   r = (R_EVAL*1e-6) * u.^(1/3);
        z  = 2*rand(np,1) - 1;   ph = 2*pi*rand(np,1);   st = sqrt(1 - z.^2);
        P  = [r.*st.*cos(ph), r.*st.*sin(ph), r.*z];
        smode = 'uniform analytic samples';
    end
    Pc   = make_Pc(e, cfg.Pc_base);
    Hhat = gI * KI;                          % ᴮĤ_I [mT/A]
    fprintf('  評估：%d 個點（%s）→ 批次 SVD\n', np, smode);

    sig = zeros(np, 3);
    CH  = 200000;                            % 分塊避免 3×6×N 陣列吃光記憶體
    pb  = P / l_hat;
    for a = 1:CH:np
        b  = min(a+CH-1, np);   n = b - a + 1;
        A  = zeros(3, 6, n);
        for k = 1:6
            d = pb(a:b,:) - Pc(:,k).';                   % n×3（無因次）
            A(:,k,:) = permute(d ./ (vecnorm(d,2,2).^3), [2 3 1]);
        end
        sv = pagesvd(pagemtimes(A, Hhat), 'vector');     % 3×1×n
        sig(a:b,:) = permute(sv, [3 1 2]);
    end

    S = struct('model',MODEL, 'geom',GEOM, 'variant',VARIANT, 'R_fit',R_FIT, 'R_eval',R_EVAL, ...
               'USE_BIAS',USE_BIAS, 'Nc',Nc, 'tri',tri, 'npts',np, 'smode',smode, ...
               'l_hat',l_hat, 'e',e, 'gI',gI, 'KI',KI, 'NMAE',rm.NMAE, 'sig',sig);
end

% ============================================================================
function render(S, figdir, tag, msfx, XR, BINW, NTK)
% 三組 σ 疊圖（各自正規化）。配色用 house 三色：深藍 / 紅 / 紫。
    % 直方圖（離散長條，照 figure-style「分布圖一律用離散長條、不要連續曲線 / KDE」）。
    %   [2026-08-18] 曾短暫改成連續曲線（plot 連 bin 中心），使用者拍板**改回長條**。
    %   [2026-08-21] 預設評估點改回**真實 .fld 格點**（65 353 點）→ 每 bin 幾百個樣本，
    %   長條帶可見的統計雜訊；要平滑包絡就傳 NSAMP（如 2e6）改用亂數加密取樣。
    % [MODIFIED 2026-09-02 使用者指定] 套上 figure-style 的九條規則：
    %   規則1 刻度數字 60（原本 FS=28 同時餵刻度與軸標題 → 拆成 FS / FSLAB）；
    %   規則3 框線加粗 LWBOX=5.0；規則6 畫布等邊（英吋正方形 + print）。
    %   軸標題 FSLAB 規則未指定值 → 沿用本腳本原值 36。
    %   規則2 圖例 45、規則7 圖例框與座標框同粗。
    FS = 60;   FSLAB = 36;   FSLEG = 45;   LWBOX = 5.0;   CANV = 14.5;   ALPH = 0.60;   nb = 180;
    COL = { [0.05 0.10 0.95], [0.85 0.10 0.10], [0.482 0.322 0.671] };
    v   = {S.sig(:,1), S.sig(:,2), S.sig(:,3)};

    allv = S.sig(:);
    if nargin < 6 || isempty(BINW)
        edg = linspace(min(allv), max(allv), nb+1);   % 三組共用 edges（固定 180 根）
    else
        edg = min(allv) : BINW : (min(allv) + ceil((max(allv)-min(allv))/BINW)*BINW);
        nb  = numel(edg) - 1;                          % 明給 bin 寬（跨圖共用同一把尺）
    end
    fprintf(['  bin 寬 %.6f mT/A（%d 根；全域 %.4f ~ %.4f）' char(10)], ...
            edg(2)-edg(1), nb, min(allv), max(allv));
    for kk = 1:3
        vk = S.sig(:,kk);   hc = histcounts(vk, edg);
        fprintf(['    sigma%d 全寬 %.4f -> %.0f 根，非空 bin 平均 %.1f 點（雜訊 ~%.0f%%）' char(10)], ...
                kk, max(vk)-min(vk), (max(vk)-min(vk))/(edg(2)-edg(1)), ...
                mean(hc(hc>0)), 100/sqrt(max(mean(hc(hc>0)),eps)));
    end
    ctr  = (edg(1:end-1) + edg(2:end))/2;

    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes(fig);   hold(ax,'on');
    h = gobjects(1,3);   pk = [];
    for k = 1:3
        p = histcounts(v{k}, edg) / numel(v{k}) * 100;
        % [MODIFIED 2026-08-21 user] no edge line on the bars (EdgeColor none):
        %   65k grid points spread over 180 narrow bars -> a 0.3pt black edge on every
        %   bar merges into a dark haze that hides the fill colour. Fill + FaceAlpha only.
        h(k) = bar(ax, ctr, p, 1, 'FaceColor',COL{k}, 'FaceAlpha',ALPH, ...
                   'EdgeColor','none');
        pk = [pk p];  %#ok<AGROW>
    end

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
           'TickLength',[.015 .015],'TickDir','out');
    % [MODIFIED 2026-08-18 使用者拍板] 橫軸兩端**貼齊有資料的範圍**（不再從 0 起、
    %   也不補到 60）：外推到整數 → [floor(min), ceil(max)]；內部刻度取整數等距奇數個。
    % [MODIFIED 2026-08-21 使用者拍板] 明給 XR 時用它（終點是乾淨的數、不被 σ₁ 的長尾
    %   拉到 282 / 59 那種值）；XR 空著才退回舊的「貼齊極值」。
    % [MODIFIED 2026-09-02] 規則 5（奇數個等距 tick、與兩端間距也相等）+ 規則 9
    %   （刻度最多一位小數）。舊的 [floor,ceil] = [13,18] 只能配 14/15/16/17（**偶數** 4 根），
    %   跨 5 個單位要奇數根且兩端間距相等只剩 s=0.5（9 根，60 pt 下必撞在一起）。
    %   axis_odd 改取 [12,18] + 13/14/15/16/17（5 根、整數、兩端間距 = 1）。
    if isempty(XR)
        [xr, xt0] = axis_odd(min(allv), max(allv), NTK);
    else
        xr = XR;
        for k = 1:3                                   % 落在視野外的比例（誠實回報）
            f = 100*mean(v{k} < xr(1) | v{k} > xr(2));
            if f > 0, fprintf('  ⚠ σ%d 有 %.3f%% 的點落在視野 [%g, %g] 之外\n', k, f, xr); end
        end
    end
    xlim(ax, xr);
    if exist('xt0','var'), xt = xt0; else, xt = xticks_in(xr); end   % 等距、奇數個、乾淨的數
    if isempty(xt)                                    % 保底：舊做法
        sx = nice_step((xr(2)-xr(1))/6);
        xt = (ceil((xr(1)+sx/2)/sx) : floor((xr(2)-sx/2)/sx)) * sx;
        if mod(numel(xt),2) == 0 && numel(xt) > 1, xt = xt(1:end-1); end
    end
    set(ax,'XTick',xt);
    [yr, yt] = ylim_odd(max(pk), NTK);            % [MODIFIED 2026-09-02] 規則 4+5
    ylim(ax, yr);   set(ax,'YTick',yt);   ytop = yr(2);

    for xv = xr                                        % x 端點只標數字、不畫 tick
        text(ax, xv, -0.022*ytop, sprintf('%g', xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    xlabel(ax, '$\mathbf{\sigma_{kk}\;(mT/A)}$',  'Interpreter','latex', 'FontSize',FSLAB);
    ylabel(ax, '$\mathbf{Percentage\;(\%)}$',    'Interpreter','latex', 'FontSize',FSLAB);

    % ---- 圖例：座標框**內右上角**、縱向三列（使用者拍板 2026-08-18）----
    %   [RESTORED 2026-09-02] 上一輪誤將本圖的圖例拿掉（拿錯張），使用者更正後還原。
    %   三個系列的圖例是必要的（圖上無法區分 sigma11/22/33）；單一系列的
    %   gain/iso solo 圖才是圖例冗餘、已於同日拿掉。
    lb = cell(1,3);
    for k = 1:3
        % ⚠ 不要用 sprintf 帶 '\sigma'：sprintf 會把 \s 當跳脫序列吃掉 -> 標籤變空白
        lb{k} = ['\sigma_{' num2str(k) num2str(k) '}'];
    end
    lg = legend(ax, h, lb, 'Interpreter','tex', 'Location','northeast', 'NumColumns',1);
    lg.FontSize = FSLEG;   lg.FontWeight = 'bold';                  % [規則2]
    lg.Box = 'on';      lg.EdgeColor = 'k';   lg.LineWidth = LWBOX; % [規則7] 圖例框與座標框同粗
    lg.Color = 'w';
    lg.ItemTokenSize = [55 25];   % 45 pt 字級下預設 30 pt 的色塊會貼到框線
    ax.Toolbar.Visible = 'off';   hold(ax,'off');

    sfx = '';  if strcmp(tag,'eighteen'), sfx = '_eighteen'; end
    out = fullfile(figdir, sprintf('sigma_hist_R%d_maxwell%s%s.png', S.R_eval, msfx, sfx));
    % [MODIFIED 2026-09-02] 規則 6：exportgraphics 會裁掉畫布四周白邊，正方形畫布
    %   匯出後就不是正方形 → 改用 print + 明確的 PaperPosition。
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
    print(fig, out, '-dpng', '-r200');                % 14.5 in x 200 dpi = 2900 px 見方
    fprintf('wrote %s\n', out);
    close(fig);
end

% ============================================================================
function [xr, xt] = axis_odd(lo, hi, NTK)
% [ADDED 2026-09-02] 水平軸：規則 5（奇數個等距 tick、且與起點/終點的間距也相等）
%   + 規則 9（刻度數字最多一位小數）。
%   形式：xr = [x0, x0+(n+1)*s]、刻度 = x0 + (1:n)*s → 兩端留白 ≡ 刻度間距。
%   排序鍵：① 視野最小（浪費最少）② 刻度是整數優先 ③ 根數接近 3 ④ 左右留白對稱。
    best = [inf inf inf inf];   xr = [floor(lo) ceil(hi)];   xt = [];
    for k = -2:3
        for c = [1 1.5 2 2.5 5]   % [MODIFIED 2026-09-02] 固定三根後需要 1.5 這級步長才湊得出合理視野
            s = c*10^k;
            if abs(s*10 - round(s*10)) > 1e-9, continue; end       % s 本身也要一位小數
            for n = NTK               % [MODIFIED 2026-09-02 使用者指定] 根數由呼叫端給（預設 3）
                span = (n+1)*s;
                if span < (hi-lo) - 1e-9, continue; end
                for m = (floor(hi/s) - n) + (-3:1)
                    x0 = m*s;   x1 = x0 + span;
                    if x0 > lo + 1e-9 || x1 < hi - 1e-9, continue; end
                    t = x0 + (1:n)*s;
                    if any(abs(t*10 - round(t*10)) > 1e-9), continue; end
                    sc = [span, double(any(abs(t-round(t)) > 1e-9)), abs(n-NTK), ...
                          abs((lo-x0) - (x1-hi))];
                    if lexlt(sc, best), best = sc;   xr = [x0 x1];   xt = t; end
                end
            end
        end
    end
end

% ============================================================================
function [lim, tk] = ylim_odd(maxv, NTK)
% [ADDED 2026-09-02] 縱軸：規則 4（起點與終點都不標）+ 規則 5（等距、且與兩端間距也
%   相等）+ 使用者指定的「三根 tick」。上緣 = 4*s、刻度 = (1:3)*s，0 與上緣不進 tk。
%
% [MODIFIED 2026-09-02 使用者回報「上面留白太多」]
%   舊版的步長只能取 [1 1.5 2 2.5 3 4 5 ...] 這組 nice 值。三根下上緣被鎖成 4*s，
%   步長一跳上緣就跳 → sigma 圖峰值 4.2%% 需 s>=1.05，舊清單只能給 1.5 → 上緣 6、
%   填充率僅 70%%，上方空一大塊。
%   現在允許**兩位有效數字**的步長（1.1 / 1.7 / 2.3 …），取滿足上緣的最小者；
%   但若 nice 步長只多浪費 <= 20%%，就還是用 nice 的（讓 2/4/6、2.5/5/7.5 這種
%   好讀的刻度不會被 1.7/3.4/5.1 取代）。
%   實測：sigma 5.5%% → 1.5/3/4.5 上緣 6 改成 1.1/2.2/3.3 上緣 4.4（填充 70%%→95%%）；
%          gain A 6.3%% 仍是 2/4/6（tight=1.7、nice=2，2/1.7=1.18 <= 1.20）。
%   縱軸不受規則 9（一位小數）約束 —— 那條只綁水平軸。
    n = NTK;                                      % [MODIFIED 2026-09-02] 由呼叫端給
    smin = 1.02*maxv/(n+1);                       % 上緣至少 1.02*maxv（長條不貼框）
    k     = floor(log10(max(smin, realmin)));
    tight = ceil(smin/10^(k-1) - 1e-9) * 10^(k-1);  % 兩位有效數字的最小合格步長
    nice  = inf;
    for kk = k-1:k+1
        for c = [1 1.5 2 2.5 3 4 5 6 8]
            s = c*10^kk;
            if s >= smin - 1e-12 && s < nice, nice = s; end
        end
    end
    s = tight;
    if nice <= 1.20*tight, s = nice; end          % nice 只貴 <=20%% 就用 nice
    tk  = (1:n)*s;
    lim = [0, (n+1)*s];
end

% ============================================================================
function xt = xticks_in(xr)
% 給定 [x0 x1] 內挑內部刻度：等距、乾淨的數、3~5 根。
%   端點本身不放 tick —— 端點值由呼叫端以 text 標。與 plot_gain_iso_hist 同一支。
%
%   [MODIFIED 2026-08-21] **優先序修正**：照 figure-style 慣例 #5，
%       整數 > 奇數個 > 接近 5 根
%   舊版把「奇數個」排在「整數」之前，[13,18] 會選出 13.5/14.5/15.5/16.5/17.5
%   （5 根、奇數，但全是小數）而捨棄 14/15/16/17（4 根、偶數，但全整數）——
%   規則明寫「為了湊整數而讓刻度變成 4 個（偶數）是可以接受的」，故改為整數優先。
%   既有圖不受影響：[4,55] → {10,20,30,40,50}（整數且 5 根，兩種排序都選它）；
%                   [5,250] → {25,75,125,175,225}（整數相位有 4 根與 5 根兩解，取 5 根）。
    span = xr(2) - xr(1);
    cand = [1 2 2.5 5 10];
    xt = [];   best = [inf inf inf];
    for k = (floor(log10(span))-2) : (floor(log10(span))+1)
        for c = cand
            s = c*10^k;
            for o = [0, s/2]
                t = (ceil((xr(1)-o)/s + 1e-9) : floor((xr(2)-o)/s - 1e-9))*s + o;
                t = t(t > xr(1)+1e-9 & t < xr(2)-1e-9);
                n = numel(t);
                if n < 3 || n > 5, continue; end
                if any(abs(t*10 - round(t*10)) > 1e-9), continue; end   % 乾淨的數（整數或 0.N）
                notint = any(abs(t - round(t)) > 1e-9);                 % 是否含小數
                sc = [double(notint), double(mod(n,2)==0), abs(n-5)];   % 三層排序鍵
                if lexlt(sc, best), best = sc;  xt = t; end
            end
        end
    end
end

% ============================================================================
function tf = lexlt(a, b)
% 字典序比較 a < b（三層排序鍵用）
    tf = false;
    for i = 1:numel(a)
        if a(i) < b(i), tf = true;  return; end
        if a(i) > b(i), return; end
    end
end

% ============================================================================
function s = nice_step(x)
% 取「>= x 的最小 nice 步長」（figure-style：等距、疏密適中、整數優先）
    n = [1 2 5 10];
    e = floor(log10(x));   c = x/10^e;
    s = n(find(n >= c - 1e-12, 1)) * 10^e;
end

% ============================================================================
function [xr, xt] = xlim_pick(lo, hi) %#ok<DEFNU>
% 橫軸：貼著資料選 nice 步長（不強迫從 0 起），取 3~7 個內部刻度；端點另以 text 標數字。
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
    s = rng_/5;   xr = [lo-0.1*rng_, hi+0.1*rng_];   xt = lo + (1:4)*s;
end

% ============================================================================
function [lim, tk] = ylim_from_zero(maxv)
% 縱軸自 0 起、上緣只留 8% 裕度；刻度等距、奇數個、且**至少 3 根**（figure-style
%   「tick 不可太擠 / 疏密適中」）。
% [MODIFIED 2026-08-21] 舊版只用 maxv/4 挑一個步長，遇到 n 被「取奇數」由 2 削成 1
%   就只剩**一根**刻度 —— R_eval=150 的 long2016 σ 圖正是如此（peak 4.2% → s=2 →
%   n=floor(4.54/2)=2 → 削成 1）。改成從該步長起往小一級試，取第一個能給 >=3 根的。
%   既有圖不受影響：R500 long2016（peak~6 → 2/4/6）與 zhi_peng（peak~8.6 → 2.5/5/7.5）
%   本來 n 就是 3，第一輪就命中、結果不變。
    cand = [1 2 2.5 3 4 5 10];
    top  = 1.08*maxv;
    x    = maxv/4;   k = floor(log10(x));
    s0   = cand(find(cand*10^k >= x, 1)) * 10^k;          % 舊版的選擇
    S    = sort([cand*10^(k-1), cand*10^k], 'descend');
    S    = S(S <= s0 + 1e-12);                            % 只往小試，不往大
    lim  = [0 top];
    for s = S
        n = floor((top - 1e-12)/s);
        if mod(n,2) == 0, n = n - 1; end
        if n >= 3, tk = (1:n)*s;  return; end
    end
    n = floor((top - 1e-12)/s0);                          % 保底：退回舊行為
    if mod(n,2) == 0, n = n - 1; end
    tk = (1:max(n,1))*s0;
end

% ============================================================================
function Pc = make_Pc(e17, Pc_base)
% 電荷格 Pc = Pc_base + E(e)（含 e6z 約束；與 solve_current 內的版本一致）
    if isempty(e17) || all(e17(:) == 0), Pc = Pc_base;  return; end
    E = zeros(3,6);
    E(:,1) = e17(1:3);    E(:,2) = e17(4:6);
    E(:,3) = e17(7:9);    E(:,4) = e17(10:12);
    E(:,5) = e17(13:15);
    E(1,6) = e17(16);     E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);
    Pc = Pc_base + E;
end
