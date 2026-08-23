function plot_sigma_hist(USE_BIAS, R_FIT, R_EVAL, NSAMP, force, MODEL, GEOM, VARIANT, XR, BINW)
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
    if nargin < 10, BINW = []; end

    tag  = 'single';   if USE_BIAS, tag = 'eighteen'; end
    msfx = '';  if ~strcmp(MODEL,'long2016_hexapole_halfcut'), msfx = ['_' MODEL]; end

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

    render(S, figdir, tag, msfx, XR, BINW);
end

% ============================================================================
function S = compute(MODEL, GEOM, VARIANT, R_FIT, R_EVAL, USE_BIAS, here, NSAMP)
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
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
    % ---- 讀 main.m 產的收斂設計校正結果（不再自己重跑階梯 + 校正）--------
    %   [MODIFIED 2026-08-23 使用者拍板] 校正與收斂判準已搬回 main.m
    %   （conv_design_ws / conv_design_sensor 只負責決定內插點位置與取場），
    %   繪圖端改成**接收 main 產完的結果** -> 圖與結果 PDF 保證出自同一次校正。
    %   ⚠ 該組合必須先跑過 main.m（GRID_NRPT='auto'）；找不到就報錯，不猜。
    %   ⚠ 同一組合可能有多顆 convN 檔（舊實驗留下的，例如 long2016 R150
    %     eighteen 就有 convN6/convN80/convN88）-> 只認 conv_auto==true 那顆。
    md_ = fullfile(CAL, 'data', MODEL, '.mat');
    tg_ = 'single';   if USE_BIAS, tg_ = 'eighteen'; end
    dd_ = dir(fullfile(md_, sprintf('calib_current_%s_convN*_R%03d_%s.mat', ...
                                    VARIANT, round(R_FIT), tg_)));
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
           numel(dd_), MODEL, R_FIT, USE_BIAS);
    cal = load(fullfile(md_, dd_(1).name));
    e = cal.e;   l_hat = cal.l_hat;   KI = cal.KI_bar;   gI = cal.gI_hat;
    tri = cal.GRID_NRPT;   Nc = cal.npts;   rm = struct('NMAE', cal.NMAE);

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
function render(S, figdir, tag, msfx, XR, BINW)
% 三組 σ 疊圖（各自正規化）。配色用 house 三色：深藍 / 紅 / 紫。
    % 直方圖（離散長條，照 figure-style「分布圖一律用離散長條、不要連續曲線 / KDE」）。
    %   [2026-08-18] 曾短暫改成連續曲線（plot 連 bin 中心），使用者拍板**改回長條**。
    %   [2026-08-21] 預設評估點改回**真實 .fld 格點**（65 353 點）→ 每 bin 幾百個樣本，
    %   長條帶可見的統計雜訊；要平滑包絡就傳 NSAMP（如 2e6）改用亂數加密取樣。
    FS = 28;   ALPH = 0.60;   nb = 180;
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

    fig = figure('Color','w','Position',[100 100 1180 860]);
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
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5, ...
           'TickLength',[.015 .015],'TickDir','out');
    % [MODIFIED 2026-08-18 使用者拍板] 橫軸兩端**貼齊有資料的範圍**（不再從 0 起、
    %   也不補到 60）：外推到整數 → [floor(min), ceil(max)]；內部刻度取整數等距奇數個。
    % [MODIFIED 2026-08-21 使用者拍板] 明給 XR 時用它（終點是乾淨的數、不被 σ₁ 的長尾
    %   拉到 282 / 59 那種值）；XR 空著才退回舊的「貼齊極值」。
    if isempty(XR)
        xr = [floor(min(allv)), ceil(max(allv))];
    else
        xr = XR;
        for k = 1:3                                   % 落在視野外的比例（誠實回報）
            f = 100*mean(v{k} < xr(1) | v{k} > xr(2));
            if f > 0, fprintf('  ⚠ σ%d 有 %.3f%% 的點落在視野 [%g, %g] 之外\n', k, f, xr); end
        end
    end
    xlim(ax, xr);
    xt = xticks_in(xr);                               % 等距、奇數個、乾淨的數
    if isempty(xt)                                    % 保底：舊做法
        sx = nice_step((xr(2)-xr(1))/6);
        xt = (ceil((xr(1)+sx/2)/sx) : floor((xr(2)-sx/2)/sx)) * sx;
        if mod(numel(xt),2) == 0 && numel(xt) > 1, xt = xt(1:end-1); end
    end
    set(ax,'XTick',xt);
    [yr, yt] = ylim_from_zero(max(pk));
    ylim(ax, yr);   set(ax,'YTick',yt);   ytop = yr(2);

    for xv = xr                                        % x 端點只標數字、不畫 tick
        text(ax, xv, -0.022*ytop, sprintf('%g', xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    xlabel(ax, '$\mathbf{\sigma_{kk}\;(mT/A)}$',  'Interpreter','latex', 'FontSize',36);
    ylabel(ax, '$\mathbf{Percentage\;(\%)}$',    'Interpreter','latex', 'FontSize',36);

    % ---- 圖例：座標框**內右上角**、縱向三列（使用者拍板 2026-08-18）----------
    %   （先前的「框外 + 切齊框寬 + 手動均勻排版」已作廢：右上角本來就是空的，
    %     放進去省版面、也不必再調欄位分布。）
    %   [MODIFIED 2026-08-21 使用者拍板] **圖例不再印 mean**，只列 sigma_kk 三個系列；
    %     三個平均值改由主函式印在 console（min / mean / max / CV 那幾行）。
    lb = cell(1,3);
    for k = 1:3
        % ⚠ 不要用 sprintf 帶 '\sigma'：sprintf 會把 \s 當跳脫序列吃掉 -> 標籤變空白
        lb{k} = ['\sigma_{' num2str(k) num2str(k) '}'];
    end
    lg = legend(ax, h, lb, 'Interpreter','tex', 'Location','northeast', 'NumColumns',1);
    lg.FontSize = 24;   lg.FontWeight = 'bold';
    lg.Box = 'on';      lg.EdgeColor = 'k';   lg.LineWidth = 2.5;
    lg.Color = 'w';
    ax.Toolbar.Visible = 'off';   hold(ax,'off');

    sfx = '';  if strcmp(tag,'eighteen'), sfx = '_eighteen'; end
    out = fullfile(figdir, sprintf('sigma_hist_R%d_maxwell%s%s.png', S.R_eval, msfx, sfx));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
    close(fig);
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
