function plot_ell_vs_npts(QTY, SRC, SAMPLING, XSCALE, RMIN)
% plot_ell_vs_npts -- paper 圖：校正參數隨「取樣點數 N」的變化（參數穩定性）
% =========================================================================
%   水平軸 = 取樣點數 N（無單位）、縱軸 = QTY：
%     QTY='ell'（預設）→ 有效長度 l_hat [um]
%     QTY='gain'       → 電流增益 g_I_hat [mT/A]
%   [ADDED 2026-08-10] 兩個量共用本腳本（同一組 sweep 快取、同一套風格），
%   照 `modify-existing-files.md`「首選參數化、別開新檔」。
%   single(1-param, 無 e) 與 eighteen(18-param bias) **兩條疊在同一張**（使用者拍板
%   2026-08-10）——兩者幾乎重合，正是「l_hat 隨 N 的行為與模型無關」的證據。
%
%   ⚠ N 不是自變數本身：N 由取樣半徑 R 決定（N ~ R^3）。本圖是把 ell_gain_vs_R 的
%   同一組 sweep 換成以 N 當水平軸呈現，回答「要多少節點 l_hat 才穩」。
%
%   **本圖只讀 plot_ell_gain_vs_R 產生的 sweep 快取、不重算**。
%
%   [MODIFIED 2026-08-10] SAMPLING 預設改 'raw'（使用者拍板）——即直接用 Maxwell
%   **0.02mm 匯出格點**（.fld header: Grid Size 0.02mm、框 +-0.6mm、61^3=226,981 點），
%   不走 'interp'（先內插到 10um 細格）。理由：interp 把 20um 格再細分 8 倍(2^3)，
%   小 R 的取樣不足被內插平滑掉 —— R<=200 的 l_hat 峰對峰 raw 71.9um vs interp 5.0um
%   (single)，振盪在 raw 才看得到。要看內插版傳 SAMPLING='interp'。
%   ⚠ 匯出格點本身也是 tet 解的內插；真正的資訊上限是球內 tet 數
%     （Sphere1mm 0.1mm 網格 134,067 tets => 平均 tet 邊長 ~57um；R=150um 球內約 450 tets）。
%
%   風格 = 選項①粗體框 @ paper scale，但**不放軸標題**（使用者 2026-08-10 拍板：
%   刻度不需要標題）——只留刻度數字。其餘照規則：box on / grid off / FS 36 粗體 /
%   刻度奇數等距 / 曲線首末點貼齊左右框邊 / 水平軸起點與終點以 text 補數字。
%
%   輸出 → figures/paper_fig/Section2_E/<qty>_vs_npts_<src>.png（覆蓋迭代）。
% =========================================================================
    clc;
    if nargin < 1 || isempty(QTY),      QTY     = 'ell';     end   % 'ell' | 'gain'
    if nargin < 2 || isempty(SRC),      SRC     = 'maxwell'; end   % 'apdl' | 'maxwell'
    if nargin < 3 || isempty(SAMPLING), SAMPLING = 'raw'; end      % 'raw'=匯出格點原值 | 'interp'=10um 細格
    if nargin < 4 || isempty(XSCALE),   XSCALE  = 'log';     end   % 'log' | 'linear'
    if nargin < 5 || isempty(RMIN),     RMIN    = 40;        end   % 起始取樣半徑 [um]
    TOL = 0.01;                                                    % 穩定區間判準：偏離自身平台 <= 1%
    switch lower(QTY)
        case 'ell',  fld = 'ell_R';  unit = 'um';    fmt = '%.2f';
        case 'gain', fld = 'gI_R';   unit = 'mT/A';  fmt = '%.4f';
        otherwise,   error('QTY 必為 ''ell'' | ''gain''');
    end

    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    sstr = ''; if strcmpi(SRC,'maxwell'),     sstr = '_maxwell'; end
    istr = ''; if strcmpi(SAMPLING,'interp'), istr = '_interp';  end

    % ---- 讀兩顆 sweep 快取（single + eighteen）；ell 與 gain 都讀，供共同穩定區間 ----
    [N1, R1, E1, G1] = load_sweep(here, sstr, istr, '',      RMIN);   % single
    [N2, R2, E2, G2] = load_sweep(here, sstr, istr, '_bias', RMIN);   % eighteen
    assert(isequal(N1, N2), '兩顆快取的 N(R) 不一致 → 取樣點雲不同，不可疊圖');
    if strcmp(fld,'ell_R'), L1 = E1;  L2 = E2;  else, L1 = G1;  L2 = G2;  end

    fprintf(['  single  : R=%g..%g um, N=%d..%d, ' QTY '=' fmt '..' fmt ' %s\n'], ...
            R1(1), R1(end), N1(1), N1(end), L1(1), L1(end), unit);
    fprintf(['  eighteen: R=%g..%g um, N=%d..%d, ' QTY '=' fmt '..' fmt ' %s\n'], ...
            R2(1), R2(end), N2(1), N2(end), L2(1), L2(end), unit);
    k = R1 <= 200;                                   % 小取樣區的振盪幅度（判斷「幾點才穩」的指標）
    fprintf(['  R<=200um 峰對峰：single ' fmt ' / eighteen ' fmt ' %s\n'], ...
            max(L1(k))-min(L1(k)), max(L2(k))-min(L2(k)), unit);

    % ---- 畫圖（單 panel，選項①粗體框）----
    c1 = [0.05 0.10 0.95];   c2 = [0.85 0.10 0.10];    % single 亮藍 / eighteen 紅
    FS = 36;   LW = 3.0;   MS = 8;
    fig = figure('Color','w','Position',[100 100 1120 820]);
    ax  = axes(fig);   hold(ax,'on');

    h1 = plot(ax, N1, L1, '-o', 'Color',c1, 'LineWidth',LW, 'MarkerSize',MS, ...
              'MarkerFaceColor',c1, 'Clipping','off');
    h2 = plot(ax, N2, L2, '-s', 'Color',c2, 'LineWidth',LW, 'MarkerSize',MS, ...
              'MarkerFaceColor',c2, 'Clipping','off');

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.02 .02]);
    ax.Toolbar.Visible = 'off';

    % ---- 水平軸：曲線首末點貼齊左右框邊（figure-style 2026-08-06）----
    XL = [min(N1), max(N1)];   xlim(ax, XL);
    set(ax, 'XScale', XSCALE);
    if strcmpi(XSCALE,'log')
        % 內部刻度取整十次冪（log 空間等距、奇數個）；起訖數字另用 text 補
        e  = ceil(log10(XL(1))) : floor(log10(XL(2)));
        XT = 10.^e;
        XT = XT(XT > XL(1) & XT < XL(2));
        if mod(numel(XT),2) == 0, XT = XT(1:end-1); end      % 刻度數量取奇數
    else
        XT = inner_ticks(XL(1), XL(2), 3);                    % 3 個等距內部刻度
    end
    set(ax,'XTick',XT);

    % ---- 縱軸：奇數等距刻度 + 兩端留白 = 間距（首末不標數字）----
    [YL, YT] = axlim_auto(min([L1 L2]), max([L1 L2]), [3 5]);
    ylim(ax, YL);   set(ax,'YTick',YT);

    % ---- 穩定區間（虛線標出）----------------------------------------------
    % 判準：最長的**連續**視窗，使 **ell 與 gain 的 single/eighteen 共四條曲線**各自落在
    % 該窗內「自己的平均」的 ±TOL 內。**兩張圖共用同一個區間**（使用者拍板 2026-08-10）。
    % [MODIFIED 2026-08-10] 演進：①兩條落在「共同平台」→ 對 gain 失效（single 9.6 vs
    % eighteen 10.8，本來就差 11.7%）②改 per-curve 但各量各算 → ell R100~200 / gain
    % R220~320 不相交，且 TOL 動 0.5pp 結論就翻（±1% 時反而大幅重疊）＝不穩健。
    % ③現行：四條共同求解、一個區間，TOL=1% ⇒ N 494~7,156（R 100~240），四條各自 <0.9%。
    % 界線由資料算出，不是目測；改 TOL 即改嚴格度。
    [i1, i2] = stable_window({E1, E2, G1, G2}, TOL);
    fprintf(['  穩定區間（ell+gain 共同）：N = %d ~ %d（R = %g ~ %g um）、本圖平台 single ' ...
             fmt ' / eighteen ' fmt ' %s（判準 ±%.1f%%）\n'], ...
             N1(i1), N1(i2), R1(i1), R1(i2), mean(L1(i1:i2)), mean(L2(i1:i2)), unit, TOL*100);
    for xv = [N1(i1), N1(i2)]
        xline(ax, xv, '--', 'Color',[0.25 0.25 0.25], 'LineWidth',2.5, 'HandleVisibility','off');
    end

    % ---- 水平軸起點 + 終點補數字（只標數字、不畫 tick mark）----
    yoff = YL(1) - 0.022*diff(YL);
    for xv = XL
        text(ax, xv, yoff, sprintf('%g', xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end

    % [MODIFIED 2026-08-11] 使用者要求加回軸標題（照 figure-style：LaTeX \mathbf、FS 36、單位括號）
    xlabel(ax, '$\mathbf{Number\;of\;points}$', 'Interpreter','latex', 'FontSize',FS);
    if strcmp(fld,'ell_R')
        ylabel(ax, '$\mathbf{\hat{\ell}\;(micro\;meter)}$', 'Interpreter','latex', 'FontSize',FS);
    else
        ylabel(ax, '$\mathbf{{}^{B}\hat{g}_{I}\;(mT/A)}$',  'Interpreter','latex', 'FontSize',FS);
    end

    % ---- 圖例：照 figure-style.md「圖例標準樣式」（2026-08-10 拍板，範本 plot_sensor_B_hist）----
    lg = legend(ax, [h1 h2], {'Single parameter', 'Eighteen parameters'}, ...
                'Interpreter','tex', 'Location','northoutside', 'NumColumns',2);
    lg.FontSize = 24;  lg.FontWeight = 'bold';  lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;

    ax.Toolbar.Visible = 'off';  hold(ax,'off');

    % 對齊：圖例左右緣切齊座標框、置於框正上方固定間距（northoutside 預設不對齊，須手動定位）
    drawnow;
    axp = get(ax,'Position');  lgh = get(lg,'Position');  lgh = lgh(4);
    GAPN = 0.022;  newTop = 1 - lgh - GAPN - 0.006;
    axp(4) = newTop - axp(2);  set(ax,'Position',axp);
    set(lg, 'Position', [axp(1), newTop + GAPN, axp(3), lgh]);

    % [MODIFIED 2026-08-11] 檔名帶取樣模式後綴，否則 raw 與 interp 兩版**互相覆蓋**且從檔名
    %   看不出是哪個（與交接記載的 npts_cost_vs_R 同一個 bug）。無後綴 = raw（不改既有檔名）。
    out = fullfile(figdir, sprintf('%s_vs_npts_%s%s.png', lower(QTY), lower(SRC), istr));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('\nwrote %s\n', out);
end

% ============================================================================
function [N, R, E, G] = load_sweep(here, sstr, istr, bstr, RMIN)
% 讀一顆 ell_gain sweep 快取，回傳 R >= RMIN 的 (N, R[um], ell[um], gI[mT/A])。
%   ell 與 gain **兩個量都回傳**：畫其中一個，但穩定區間要四條曲線一起算。
    f = fullfile(here, 'data', sprintf('ell_gain_sweep%s%s%s.mat', sstr, istr, bstr));
    assert(exist(f,'file')==2, ...
        'sweep 快取不存在：%s\n先跑 plot_ell_gain_vs_R 重建。', f);
    S = load(f);
    m = S.Rum >= RMIN;
    R = S.Rum(m);   N = S.npts_R(m);   E = S.ell_R(m);   G = S.gI_R(m);
    R = R(:).';     N = N(:).';        E = E(:).';       G = G(:).';   % 一律 row
    fprintf('loaded %s\n', f);
end

% ============================================================================
function [i1, i2] = stable_window(C, TOL)
% 最長的連續視窗 [i1,i2]，使 cell C 內**每條曲線各自**落在該窗內自己的平均的 ±TOL 內。
%   per-curve 而非共同平台：各曲線可以有不同的平台值（gain 的兩個模型就差 11.7%），
%   「穩定」問的是各自平不平，不是彼此合不合。
%   窮舉所有視窗（n=24，O(n^2) 可忽略），取最長者；同長度取先出現者。
    n = numel(C{1});   best = 0;   i1 = 1;   i2 = 1;
    for a = 1:n
        for b = a:n
            d = 0;
            for c = 1:numel(C)
                p = mean(C{c}(a:b));
                d = max(d, max(abs(C{c}(a:b)-p))/p);
            end
            if d <= TOL && (b-a+1) > best
                best = b-a+1;   i1 = a;   i2 = b;
            end
        end
    end
end

% ============================================================================
function tk = inner_ticks(lo, hi, n)
% 線性軸的 n 個等距內部刻度（不含首末端點；n 取奇數）。
    cand = [1 2 2.5 3 4 5 10];
    rng_ = hi - lo;
    s = rng_/(n+1);
    k = floor(log10(s));
    [~, i] = min(abs(cand*10^k - s));
    s = cand(i)*10^k;
    ctr = round(((lo+hi)/2)/s)*s;
    tk  = ctr + (-(n-1)/2 : (n-1)/2)*s;
    tk  = tk(tk > lo & tk < hi);
end

% ============================================================================
function [lim, tk] = axlim_auto(lo, hi, nlist)
% 通用軸範圍/刻度：**奇數個等距 tick，且兩端留白 = tick 間距**（figure-style 拍板）。
%   與 plot_npts_cost_vs_R.m 的同名 local 一致（house style：各腳本自帶一份）。
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
                    L_  = [t(1)-s, t(end)+s];                % 兩端留白 = s
                    clr = 0.15*s;                            % 資料離框至少 0.15 格
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
