function plot_ell_vs_npts_fixedR(RFIX, MODELS, NN, ORDER, force)
% plot_ell_vs_npts_fixedR -- paper 圖：固定取樣半徑 R，只增加取樣點數 n，看 l_hat 是否穩定
% =========================================================================
%   水平軸 = 取樣點數 n（log）、縱軸 = l_hat [um]。
%   MODELS = 'both'（預設，single + eighteen 疊圖）| 'single' | 'eighteen'
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
    if nargin < 4 || isempty(ORDER),  ORDER  = 'radial'; end
    if nargin < 5 || isempty(force),  force  = false;  end
    % 穩定判準線：n = 1588 是 **random 序**下 single 進 ±0.26% 的門檻；radial 序不適用 → 不畫。
    NMARK = 1588;
    if strcmpi(ORDER,'radial'), NMARK = []; end

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
        [nl{i}, el{i}] = sweep_one(RFIX, bias(i), NN, ORDER, force, here);
    end

    % ---- 畫圖（選項①粗體框；無軸標題）--------------------------------------
    %   MODELS='both' → **上下兩 panel、各自 y 尺度**（使用者拍板 2026-08-11：疊在同一
    %   尺度看不出 eighteen 的結構 —— 它的擺幅只有 single 的 1/14）。依 figure-style
    %   「兩 panel 的 y tick 數量必須一致」，先取共同 tick 數再各自算範圍。
    FS = 36;  LWBOX = 4.0;
    cols = [0.05 0.10 0.95; 0.85 0.10 0.10];   % single 亮藍 / eighteen 紅
    mks  = {'-o','-s'};
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
        ylabel(tl, '$\mathbf{\hat{\ell}\;(micro\;meter)}$',    'Interpreter','latex', 'FontSize',FS);
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
        ylabel(axleg, '$\mathbf{\hat{\ell}\;(micro\;meter)}$', 'Interpreter','latex', 'FontSize',FS);
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
    out = fullfile(figdir, sprintf('ell_vs_npts_fixedR%d_maxwell10_%s%s.png', RFIX, lower(MODELS), ostr));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('\nwrote %s\n', out);
end

% ============================================================================
function [nlist, ell_n] = sweep_one(RFIX, USE_BIAS, NN, ORDER, force, here)
% 一個模型的 n-sweep（有快取就讀）。回 nlist 與 ell_n [um]。
    HGRID = 10e-6;   N0 = 20;
    mstr   = 'single'; if USE_BIAS, mstr = 'eighteen'; end
    ostr   = ''; if strcmpi(ORDER,'radial'), ostr = '_radial'; end   % random 序不加後綴（沿用舊快取）
    cachef = fullfile(here, 'data', sprintf('ell_vs_npts_fixedR%d_maxwell10_%s%s.mat', RFIX, mstr, ostr));
    if exist(cachef,'file') && ~force
        S = load(cachef);   nlist = S.nlist;   ell_n = S.ell_n;
        fprintf('loaded cache %s\n', cachef);   return;
    end

    solver_path();
    cfg = model_config('long2016_hexapole_halfcut','tip40um');
    raw = extract_maxwell_data(cfg, 'all', cfg.default_variant);
    ad  = build_actuator_data(raw, cfg);
    ad  = interp_grid_sample(ad, cfg, RFIX*1e-6, HGRID);      % → 10 um 格、只到 RFIX

    Pool  = ad.Pa;   Bpool = ad.Ba;   npool = size(Pool,1);
    if strcmpi(ORDER,'radial')
        [~, ord] = sort(ad.r2, 'ascend');                     % 由內往外（確定性、無抽樣運氣）
    else
        rng(0);   ord = randperm(npool);                      % 打亂一次，之後不再抽
    end
    nlist = unique(round(logspace(log10(N0), log10(npool), NN)));
    fprintf('  [%s/%s] 取樣池 %d 點（R<=%d um, %g um 格）；測 %d 個點數 %d..%d\n', ...
            mstr, lower(ORDER), npool, RFIX, HGRID*1e6, numel(nlist), nlist(1), nlist(end));

    ell_n = nan(size(nlist));   reff = nan(size(nlist));
    fprintf('\n      n     ell_hat[um]   r_eff[um]\n');
    for k = 1:numel(nlist)
        idx = ord(1:nlist(k));                                % **前 n 個** → 巢狀累積
        P   = Pool(idx, :);
        Bstack = zeros(3*nlist(k), cfg.N_I);
        for j = 1:cfg.N_I
            Bstack(:,j) = reshape(Bpool(idx,:,j).', [], 1);   % 同 select_ball 的堆疊
        end
        [~, l_hat] = fitting(P, Bstack, ad.Pc_base, 0.5e-3, USE_BIAS);
        ell_n(k) = l_hat*1e6;
        reff(k) = sqrt(max(ad.r2(idx)))*1e6;                  % 該 n 的有效半徑 [um]
        fprintf('  %7d   %10.3f      %6.1f\n', nlist(k), ell_n(k), reff(k));
    end
    save(cachef, 'nlist', 'ell_n', 'reff', 'npool', 'RFIX', 'USE_BIAS', 'HGRID', 'ORDER');
    fprintf('saved cache %s\n', cachef);
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
