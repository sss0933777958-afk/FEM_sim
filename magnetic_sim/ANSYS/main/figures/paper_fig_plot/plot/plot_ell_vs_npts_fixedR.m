function plot_ell_vs_npts_fixedR(RFIX, MODELS, NN, ORDER, force, QTY, MREP)
% plot_ell_vs_npts_fixedR -- paper 圖：固定取樣半徑 R，只增加取樣點數 n，看參數是否穩定
% =========================================================================
%   水平軸 = 取樣點數 n（log）、縱軸 = QTY：
%     QTY='ell' （預設）→ 有效長度 l_hat [um]      顏色 **藍**
%     QTY='gain'        → 電流增益 g_I_hat [mT/A]  顏色 **紅**
%   MODELS = 'both'（預設，上下兩 panel）| 'single' | 'eighteen'
%
%   [MODIFIED 2026-08-12] 使用者拍板兩件事：
%     ① **一張 l_hat、一張 gain**（不是 random/radial 兩張 l_hat）；兩張都上下分
%        single / eighteen 兩個子圖。
%     ② 曲線 = **MREP=10 組隨機抽樣的逐點中位數**（不是單一 seed）——單組帶巨大抽樣運氣
%        （同一 R 的 10 組 N_min 差 5~17 倍），中位數才是可引用的統計量。
%     ③ 視覺編碼：**顏色 = 物理量**（l_hat 藍 / gain 紅）、**marker = 模型**
%        （single 圓 o / eighteen 方 s）。
%
%   ⚠ 與 plot_ell_vs_npts 的差別（兩者橫軸都是點數，意義完全不同）：
%     plot_ell_vs_npts        n 隨 R 一起變（n ~ R^3）→ 點數與區域大小綁在一起，分不開
%     plot_ell_vs_npts_fixedR **R 固定**、只有 n 變 → 自變數乾淨，直接回答「幾點才夠」
%
%   取樣（使用者拍板 2026-08-11）：
%     ① 取樣池 = R <= RFIX 球內的 **10 um 內插格**格點（RFIX=150 時 14,147 個）。
%        「10 um 內插」= 把 Maxwell 0.02mm 匯出格內插到 10 um 均勻格（同 plot_ell_gain_vs_R
%        的 interp_grid_sample）。⚠ 內插不增加資訊，論文圖說須標明場經內插取樣。
%     ② `rng(0)` + randperm 把池**隨機打亂一次**（固定種子 → 可重現）。
%     ③ 第 k 個點數 n_k 的點集 = 該排列的**前 n_k 個** —— **巢狀累積**：
%        n 變大時舊點全留、只新增，不重抽、不重複取樣。
%        ⇒ 曲線的變化純粹來自「加點」，不含「換了一組點」的雜訊。
%     ④ n_k = 20 → 池大小，log 間距 NN 個（預設 40，每步約 +18%）。每個 n 只擬合一次。
%
%   已知結果（RFIX=150, single）：n >= 1588 之後 l_hat = 875.12 um、峰對峰 0.257%
%   （±0.128%），而點數還能再增加 8.9 倍。n < 60 完全不可用（20~55 在 864.6~898.7 亂跳）。
%
%   風格 = 選項①粗體框 @ paper scale；**不放軸標題**（同 ell_vs_npts_maxwell.png）。
%   ⚠ R <= RFIX 的資訊只在檔名，不在圖上（legend 照標準樣式只放兩個模型名）。
%   輸出 → figures/paper_fig/Section2_E/ell_vs_npts_fixedR<R>_maxwell10_<models>.png
% =========================================================================
    clc;
    if nargin < 1 || isempty(RFIX),   RFIX   = 150;    end   % 固定取樣半徑 [um]
    if nargin < 2 || isempty(MODELS), MODELS = 'both'; end   % 'both' | 'single' | 'eighteen'
    if nargin < 3 || isempty(NN),     NN     = 40;     end   % 點數取樣個數（log 間距）
    % [ADDED 2026-08-11] 取點順序（使用者拍板改問法）：
    %   'random'（原本）＝把池隨機打亂、取前 n → n 點**散布在整個 R 球內**，
    %                    問的是「這個範圍要取多密」。⚠ 結果依抽樣運氣而變。
    %   'radial'（新）  ＝把池**依半徑由小到大排序**、取前 n → 前 n 點 = 半徑 <= r(n) 的
    %                    **全部**點，問的是「要取到多外面才夠」。**確定性、無抽樣抖動**。
    %   ⚠ radial 下 l_hat(n) 走的是「l_hat 對半徑」的曲線（n=20 -> r~17um、
    %     n=pool -> r=RFIX），不是「同一個球內加密」。兩者問的是不同問題，別混用。
    if nargin < 4 || isempty(ORDER),  ORDER  = 'random'; end
    if nargin < 5 || isempty(force),  force  = false;  end
    % [ADDED 2026-08-12] QTY 與抽樣組數。radial 是確定性序 → 強制 MREP=1（中位數無意義）。
    if nargin < 6 || isempty(QTY),    QTY    = 'ell';   end   % 'ell' | 'gain'
    if nargin < 7 || isempty(MREP),   MREP   = 10;      end   % 隨機序的抽樣組數（取逐點中位數）
    if strcmpi(ORDER,'radial'), MREP = 1; end
    switch lower(QTY)
        case 'ell',  ylab = '$\mathbf{\hat{\ell}\;(micro\;meter)}$';       CQ = [0.05 0.10 0.95];
        case 'gain', ylab = '$\mathbf{{}^{B}\hat{g}_{I}\;(mT/A)}$';        CQ = [0.85 0.10 0.10];
        otherwise,   error('QTY 必為 ''ell'' | ''gain''');
    end

    here   = fileparts(fileparts(mfilename('fullpath')));      % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    switch lower(MODELS)
        case 'single',   bias = false;          lbl = {'Single parameter'};
        case 'eighteen', bias = true;           lbl = {'Eighteen parameters'};
        case 'both',     bias = [false true];   lbl = {'Single parameter','Eighteen parameters'};
        otherwise, error('MODELS 必為 ''both'' | ''single'' | ''eighteen''');
    end

    nl = cell(1,numel(bias));   el = cell(1,numel(bias));
    for i = 1:numel(bias)
        [nl{i}, el{i}] = sweep_median(RFIX, bias(i), NN, ORDER, MREP, QTY, force, here);
    end

    % [MODIFIED 2026-08-12] 穩定判準線改用**本 R 實際補算的 N_min**（原本寫死 n=1588）。
    %   演算法與 plot_nmin_vs_R 的隨機序逐條相同（見 local nmin_canonical）：
    %   10 um 內插格池 → rng(d-1) 打亂、巢狀累積取前 n（NN=45 log 間距）→ single 擬合
    %   → 平坦段起點（其後 KTAIL=7 點相對變化率 <= TOLPC=0.20%）→ MREP=10 組取中位數。
    %   ⚠ radial 序不適用（少取點 = 縮小範圍，N_min 與 R 無關且不是 R 的答案）→ 不畫線。
    if strcmpi(ORDER,'radial')
        NMARK = [];
    else
        NMARK = nmin_canonical(RFIX, here, force);
    end

    % ---- 畫圖（選項①粗體框；無軸標題）--------------------------------------
    %   MODELS='both' → **上下兩 panel、各自 y 尺度**（使用者拍板 2026-08-11：疊在同一
    %   尺度看不出 eighteen 的結構 —— 它的擺幅只有 single 的 1/14）。依 figure-style
    %   「兩 panel 的 y tick 數量必須一致」，先取共同 tick 數再各自算範圍。
    FS = 36;  LWBOX = 4.0;
    % [MODIFIED 2026-08-12] 視覺編碼改為（使用者拍板）：
    %   **顏色 = 物理量**（l_hat 藍、gain 紅），**marker = 模型**（single 圓 o、eighteen 方 s）。
    %   一張圖只畫一個 QTY → 上下兩 panel 同色（CQ），靠圓/方區分 single / eighteen。
    %   （原本是「顏色 = 模型」：single 藍 / eighteen 紅，與「gain 用紅」衝突。）
    cols  = [CQ; CQ];            % 同一張圖只有一個物理量 → 兩條同色
    mks   = {'-o','-s'};         % single = 圓、eighteen = 方
    if strcmpi(MODELS,'eighteen'), cols = cols(2,:);  mks = mks(2); end
    XL = [min(cellfun(@(v)v(1),nl)) max(cellfun(@(v)v(end),nl))];
    XT = log_ticks(XL);
    two = numel(bias) == 2;

    if two
        fig = figure('Color','w','Position',[100 40 1120 1060]);
        tl  = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
        NY  = pick_common_n(el{1}, el{2}, [3 5]);     % 兩 panel 共用的 tick 數
        axs = gobjects(1,2);
        for i = 1:2
            axs(i) = nexttile(tl);   hold(axs(i),'on');
            h = plot(axs(i), nl{i}, el{i}, mks{i}, 'Color',cols(i,:), 'LineWidth',3.0, ...
                     'MarkerSize',6, 'MarkerFaceColor',cols(i,:), 'Clipping','off');
            if i == 1                                  % 圖例掛上 panel，用假線代表下 panel
                hd = plot(axs(1), nan, nan, mks{2}, 'Color',cols(2,:), 'LineWidth',3.0, ...
                          'MarkerSize',6, 'MarkerFaceColor',cols(2,:));
                hleg = [h hd];
            end
            style_panel(axs(i), FS, LWBOX);
            xlim(axs(i), XL);  set(axs(i),'XScale','log','XTick',XT,'XMinorTick','off');
            [YL, YT] = axlim_auto(min(el{i}), max(el{i}), NY);
            ylim(axs(i), YL);  set(axs(i),'YTick',YT);
            if ~isempty(NMARK)                         % 穩定判準線（兩 panel 同一 n）
                xline(axs(i), NMARK, '--', 'Color',[0.25 0.25 0.25], ...
                      'LineWidth',2.5, 'HandleVisibility','off');
            end
            if i == 1
                axs(i).XTickLabel = {};                % 上 panel x 數字隱藏（共用軸）
            else
                yoff = YL(1) - 0.030*diff(YL);         % 端點只標數字、不畫 tick
                for xv = XL
                    text(axs(i), xv, yoff, sprintf('%g',xv), 'HorizontalAlignment','center', ...
                         'VerticalAlignment','top','FontSize',FS,'FontWeight','bold','Clipping','off');
                end
            end
            hold(axs(i),'off');
        end
        % [MODIFIED 2026-08-11] 使用者要求加回軸標題。兩 panel 同為 l_hat → 用 layout 級
        %   共用 ylabel（各 panel 各標一次會重複）；xlabel 也掛 layout（共用 x 軸）。
        xlabel(tl, '$\mathbf{Number\;of\;points}$',            'Interpreter','latex', 'FontSize',FS);
        ylabel(tl, ylab,    'Interpreter','latex', 'FontSize',FS);
        axleg = axs(1);
    else
        fig = figure('Color','w','Position',[100 100 1120 820]);
        axleg = axes(fig);   hold(axleg,'on');
        hleg = plot(axleg, nl{1}, el{1}, mks{1}, 'Color',cols(1,:), 'LineWidth',3.0, ...
                    'MarkerSize',6, 'MarkerFaceColor',cols(1,:), 'Clipping','off');
        style_panel(axleg, FS, LWBOX);
        xlim(axleg, XL);  set(axleg,'XScale','log','XTick',XT,'XMinorTick','off');
        [YL, YT] = axlim_auto(min(el{1}), max(el{1}), [3 5]);
        ylim(axleg, YL);  set(axleg,'YTick',YT);
        if ~isempty(NMARK)
            xline(axleg, NMARK, '--', 'Color',[0.25 0.25 0.25], ...
                  'LineWidth',2.5, 'HandleVisibility','off');
        end
        yoff = YL(1) - 0.022*diff(YL);
        for xv = XL
            text(axleg, xv, yoff, sprintf('%g',xv), 'HorizontalAlignment','center', ...
                 'VerticalAlignment','top','FontSize',FS,'FontWeight','bold','Clipping','off');
        end
        xlabel(axleg, '$\mathbf{Number\;of\;points}$',         'Interpreter','latex', 'FontSize',FS);
        ylabel(axleg, ylab, 'Interpreter','latex', 'FontSize',FS);
        hold(axleg,'off');
    end

    % 圖例：照 figure-style.md「圖例標準樣式」（Interpreter 必用 tex）
    lg = legend(axleg, hleg, lbl, 'Interpreter','tex', 'Location','northoutside', ...
                'NumColumns', numel(bias));
    lg.FontSize = 24;  lg.FontWeight = 'bold';
    lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;
    drawnow;
    if two
        % ⚠ tiledlayout 內**不可**手動設 axes Position（MATLAB 會警告並忽略）
        %   → 用 layout 自己的 north tile，legend 會橫跨整個 layout 寬度、自動對齊。
        lg.Layout.Tile = 'north';
    else
        axp = get(axleg,'Position');  lgh = get(lg,'Position');  lgh = lgh(4);
        GAPN = 0.022;  newTop = 1 - lgh - GAPN - 0.006;
        axp(4) = newTop - axp(2);  set(axleg,'Position',axp);
        set(lg, 'Position', [axp(1), newTop + GAPN, axp(3), lgh]);
    end

    ostr = ''; if strcmpi(ORDER,'radial'), ostr = '_radial'; end
    out = fullfile(figdir, sprintf('%s_vs_npts_fixedR%d_maxwell10_%s%s.png', lower(QTY), RFIX, lower(MODELS), ostr));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('\nwrote %s\n', out);
end

% ============================================================================
function [nlist, val] = sweep_median(RFIX, USE_BIAS, NN, ORDER, MREP, QTY, force, here)
% [MODIFIED 2026-08-12] 一個模型的 n-sweep，**MREP 組隨機抽樣 → 回逐點中位數**。
%   一次算 l_hat 與 g_I_hat 兩個量、存同一顆快取（QTY 只決定回傳哪一個）
%   → 「l_hat 一張 + gain 一張」只需跑一次擬合。
%   每組 d：rng(d-1) 打亂池 → 巢狀累積取前 n（只增不換）→ fitting 得 l_hat、
%           再 solve_current 得 g_I_hat（F = identity，Maxwell map 就是 identity）。
%   ORDER='radial' 時 MREP 已被強制為 1（確定性序，無抽樣運氣可平均）。
    HGRID = 10e-6;   N0 = 20;
    mstr   = 'single'; if USE_BIAS, mstr = 'eighteen'; end
    ostr   = ''; if strcmpi(ORDER,'radial'), ostr = '_radial'; end
    cachef = fullfile(here, 'data', ...
             sprintf('npts_fixedR%d_maxwell10_%s_m%d%s.mat', RFIX, mstr, MREP, ostr));

    if exist(cachef,'file') && ~force
        S = load(cachef);   nlist = S.nlist;   ELL = S.ell;   GAIN = S.gain;
        fprintf('loaded cache %s（%d 組）\n', cachef, size(ELL,1));
    else
        solver_path();
        cfg = model_config('long2016_hexapole_halfcut','tip40um');
        raw = extract_maxwell_data(cfg, 'all', cfg.default_variant);
        ad  = build_actuator_data(raw, cfg);
        ad  = interp_grid_sample(ad, cfg, RFIX*1e-6, HGRID);   % → 10 um 格、只到 RFIX

        Pool = ad.Pa;   Bpool = ad.Ba;   npool = size(Pool,1);
        F = zeros(6, cfg.N_I);                                  % coil→pole（Maxwell = identity）
        for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end
        nlist = unique(round(logspace(log10(N0), log10(npool), NN)));
        fprintf('  [%s/%s] 池 %d 點（R<=%d um, %g um 格）；%d 組 x %d 個 n（%d..%d）\n', ...
                mstr, lower(ORDER), npool, RFIX, HGRID*1e6, MREP, numel(nlist), nlist(1), nlist(end));

        ELL = nan(MREP, numel(nlist));   GAIN = nan(MREP, numel(nlist));
        for d = 1:MREP
            if strcmpi(ORDER,'radial')
                [~, ord] = sort(ad.r2, 'ascend');               % 由內往外（確定性）
            else
                rng(d-1);   ord = randperm(npool);              % seed = d-1（同 N_min 研究）
            end
            for k = 1:numel(nlist)
                idx = ord(1:nlist(k));                          % 巢狀累積
                P   = Pool(idx, :);
                Bstack = zeros(3*nlist(k), cfg.N_I);
                for j = 1:cfg.N_I
                    Bstack(:,j) = reshape(Bpool(idx,:,j).', [], 1);
                end
                [e, l_hat] = fitting(P, Bstack, ad.Pc_base, 0.5e-3, USE_BIAS);
                [~, gI]    = solve_current(l_hat, e, ad.Pc_base, P, Bstack, F);
                ELL(d,k)  = l_hat*1e6;                          % [um]
                GAIN(d,k) = gI;                                 % [mT/A]
            end
            fprintf('    seed %d done\n', d-1);
        end
        ell = ELL;  gain = GAIN;   %#ok<NASGU>
        save(cachef, 'nlist', 'ell', 'gain', 'npool', 'RFIX', 'USE_BIAS', ...
                     'HGRID', 'ORDER', 'MREP', 'NN');
        fprintf('saved cache %s\n', cachef);
    end

    switch lower(QTY)                                           % 逐點（跨組）中位數
        case 'ell',  val = median(ELL,  1, 'omitnan');
        case 'gain', val = median(GAIN, 1, 'omitnan');
    end
end

% ============================================================================
function Nmin = nmin_canonical(RFIX, here, force)
% [ADDED 2026-08-12] 補算本 R 的 N_min —— **演算法與 plot_nmin_vs_R 的隨機序逐條相同**。
%   之所以要補：N_min 研究的網格是 R = 100:20:500，**沒有 150**（相鄰的 140→1761、
%   160→1256.5，差 1.4 倍，不能內插代用）。
%
%   演算法（= plot_nmin_vs_R 的 random_draws + flat_one + median，參數取定案值）：
%     ① 池 = R<=RFIX 球內的 10 um 內插格點（single 用；濾鐵後）
%     ② 每組 d：rng(d-1) 打亂 → **巢狀累積**取前 n（只增不換、不重抽）
%     ③ n = log 間距 NN=45 個，20 → min(pool, NCAP=15000)
%     ④ 平坦段起點 = 第一個 k，使**其後 KTAIL=7 點**的相對變化率都 <= TOLPC=0.20%
%     ⑤ MREP=10 組（seed 0..9）取**中位數**
%   ⚠ 判準參數取「定案值」而非 plot_nmin_vs_R 的函式預設（NN=40/K=5/M=5/tol=0.26）——
%     已驗證 TOL=0.20/K=7/M=10 能重現既有快取的發表值（最大 2807、佔池中位數 0.51%）。
%   ⚠ 池的列舉順序是 ndgrid 的 (z,y,x) 字典序 → 在 Rmax=RFIX 建格與「建到 500 再篩 r<=RFIX」
%     得到**同一組、同順序**的點，故同 seed 的 randperm 抽到相同子集、與研究一致。
    NN = 45;  NCAP = 15000;  N0 = 20;  HGRID = 10e-6;
    MREP = 10;  TOLPC = 0.20;  KTAIL = 7;
    cachef = fullfile(here, 'data', sprintf('nmin_fixedR%d_maxwell10_single_random.mat', RFIX));

    if exist(cachef,'file') && ~force
        S = load(cachef);   nlist = S.nlist;   E = S.ell;   npool = S.npool;
        fprintf('loaded cache %s（%d 組）\n', cachef, size(E,1));
    else
        solver_path();
        cfg = model_config('long2016_hexapole_halfcut','tip40um');
        raw = extract_maxwell_data(cfg, 'all', cfg.default_variant);
        ad  = build_actuator_data(raw, cfg);
        ad  = interp_grid_sample(ad, cfg, RFIX*1e-6, HGRID);
        npool = size(ad.Pa,1);
        nlist = unique(round(logspace(log10(N0), log10(min(npool,NCAP)), NN)));
        E = nan(MREP, numel(nlist));
        fprintf('  [N_min] R<=%d um 池 %d 點；%d 組 x %d 個 n\n', RFIX, npool, MREP, numel(nlist));
        for d = 1:MREP
            rng(d-1);   ord = randperm(npool);                 % seed = d-1（同 random_draws）
            for k = 1:numel(nlist)
                idx = ord(1:nlist(k));                         % 巢狀累積：只增不換
                P   = ad.Pa(idx,:);
                Bs  = zeros(3*nlist(k), cfg.N_I);
                for j = 1:cfg.N_I
                    Bs(:,j) = reshape(ad.Ba(idx,:,j).', [], 1);
                end
                [~, l_hat] = fitting(P, Bs, ad.Pc_base, 0.5e-3, false);   % single
                E(d,k) = l_hat*1e6;
            end
            fprintf('    seed %d done\n', d-1);
        end
        ell = E;   %#ok<NASGU>
        save(cachef, 'nlist', 'ell', 'npool', 'RFIX', 'NN', 'NCAP', 'HGRID', ...
                     'MREP', 'TOLPC', 'KTAIL');
        fprintf('saved cache %s\n', cachef);
    end

    % 每組的平坦段起點 → 中位數
    nd = nan(1, size(E,1));
    for d = 1:size(E,1)
        ell_d = E(d,:);
        dch   = [NaN, abs(diff(ell_d))./ell_d(2:end)*100];
        for a = 1:numel(nlist)-KTAIL
            if all(dch(a+1 : a+KTAIL) <= TOLPC), nd(d) = nlist(a);  break; end
        end
    end
    Nmin = median(nd, 'omitnan');
    fprintf('  [N_min] R=%d um：10 組 = %s\n', RFIX, mat2str(nd));
    fprintf('  [N_min] 中位數 = %g（佔池 %d 的 %.2f%%）\n', Nmin, npool, Nmin/npool*100);
end

% ============================================================================
function style_panel(ax, FS, LWBOX)
    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX,'TickLength',[.015 .015]);
    ax.Toolbar.Visible = 'off';
end

% ============================================================================
function n = pick_common_n(v1, v2, nlist)
% 兩 panel 共用的 y tick 數（figure-style：上下 panel tick 數量必須一致）：
%   取「較差那個 panel 的填充率最高」者。填充率 = 資料跨度 / 軸總跨距。
    best = nlist(1);   bestFill = -inf;
    for n_ = nlist
        L1 = axlim_auto(min(v1), max(v1), n_);
        L2 = axlim_auto(min(v2), max(v2), n_);
        f  = min( (max(v1)-min(v1))/diff(L1), (max(v2)-min(v2))/diff(L2) );
        if f > bestFill, bestFill = f;  best = n_; end
    end
    n = best;
end

% ============================================================================
function XT = log_ticks(XL)
% log 軸內部刻度：候選 c×10^k（c ∈ {1,2,3,5}，同 c ⇒ log 空間等距），
%   取「奇數個」且「離兩端點最遠」的那組 —— 否則端點數字會與最近的刻度標籤相撞
%   （本圖 14147 與 10^4 只差 0.15 decade，直接疊在一起）。
    L = log10(XL);
    best = [];   bestClr = -inf;
    for c = [1 2 3 5]
        k  = ceil(log10(XL(1)/c)) : floor(log10(XL(2)/c));
        t  = c * 10.^k;
        t  = t(t > XL(1) & t < XL(2));
        if isempty(t) || mod(numel(t),2) == 0, continue; end
        clr = min([log10(t(1)) - L(1), L(2) - log10(t(end))]);
        if clr > bestClr, bestClr = clr;  best = t; end
    end
    if isempty(best)
        k = ceil(L(1)):floor(L(2));   best = 10.^k;
        best = best(best > XL(1) & best < XL(2));
        if mod(numel(best),2) == 0 && ~isempty(best), best = best(1:end-1); end
    end
    XT = best;
end

% ============================================================================
function solver_path()
% 掛 Maxwell 分支、移除 APDL 分支（兩分支有同名函式，只 addpath 會被遮蔽）。
    APDL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    CAL  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    warning('off','MATLAB:rmpath:DirNotFound');
    rmpath(fullfile(APDL,'function'));  rmpath(fullfile(APDL,'common_path'));
    warning('on','MATLAB:rmpath:DirNotFound');
    addpath(fullfile(CAL,'function'));  addpath(fullfile(CAL,'common_path'));
end

% ============================================================================
function adq = interp_grid_sample(ad, cfg, Rmax, h)
% 把場內插到 h 間距的均勻格、只留球內（與 plot_ell_gain_vs_R 的同名 local 同演算法）。
%   ⚠ 內插不增加資訊，是重新取樣。
    g = -Rmax:h:Rmax;
    [X, Y, Z] = ndgrid(g, g, g);
    Pq = [X(:), Y(:), Z(:)];
    rq = vecnorm(Pq, 2, 2);
    in = rq <= Rmax;   Pq = Pq(in,:);   rq = rq(in);

    Pm  = (ad.R_act.' * Pq.').';                          % actuator → measure
    air = filter_iron_nodes(Pm(:,1), Pm(:,2), Pm(:,3) + cfg.SPH_OFST, cfg);
    Pq  = Pq(air,:);   rq = rq(air);

    src = ad.r2 < (Rmax + 60e-6)^2;                       % 源點只取球外一圈
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
                ctr = round(mid/s)*s;
                t   = ctr + (-(n-1)/2 : (n-1)/2)*s;
                L   = [t(1)-s, t(end)+s];
                clr = 0.15*s;
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
