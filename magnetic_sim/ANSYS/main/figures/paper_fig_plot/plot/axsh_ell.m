function axsh_ell(force, MODEL, GEOM, VARIANT, NMAX, Rum, XCUT, YTOP, NOREF)
%AXSH_ELL  l_hat vs point count along the actuator-axis shell ladder.
% =========================================================================
%   One script, one figure (repo rule). Reads the cache written by
%   axsh_sweep.m and draws nothing else.
%
%       solid  : l_hat from the reduced design   (single blue / eighteen red)
%       dashed : l_hat from the FULL grid inside R, same colour -- the value
%                the ladder converges towards (NOT in the legend; the legend
%                keeps only the two model names, user 2026-08-27).
%
%   Curve starts at x = 0 with the fitting initial guess l0 = 500 um.
%
%   ── Reaching the right frame edge ───────────────────────────────────
%   The ladder is 6*Nr+1 = 7, 13, 19, ... so XCUT = 200 is NOT a rung
%   (Nr=33 -> 199, Nr=34 -> 205). The line is drawn through the first rung
%   PAST the cut and the axes clip it at x = XCUT, so the curve meets the
%   frame with no gap and no invented data point (user, 2026-08-27: 單純補
%   線段，不要線性內插). Markers are drawn separately, only on rungs within
%   the cut, with clipping off so the x = 0 anchor is not sliced in half.
%
%   Style preset (1) bold-framed: FS 36 bold, box on, grid off, LINEAR x,
%   x endpoints 0 / XCUT labelled with text (no tick mark), internal ticks
%   integer, y bottom flush at 500.
%
%   Output -> figures/paper_fig/Section2_E/ell_vs_npts_axsh_R<Rum>[_<model>][_noref].png
% =========================================================================
    clc;
    if nargin < 1 || isempty(force),   force   = false;                       end
    if nargin < 2 || isempty(MODEL),   MODEL   = 'long2016_hexapole_halfcut'; end
    if nargin < 3 || isempty(GEOM),    GEOM    = 'tip40um';                   end
    if nargin < 4,                     VARIANT = '';                          end
    if nargin < 5 || isempty(NMAX),    NMAX    = 500;                         end
    if nargin < 6 || isempty(Rum),     Rum     = 150;                         end
    if nargin < 7 || isempty(XCUT),    XCUT    = 50;                          end
    % [ADDED 2026-08-27 使用者指定] YTOP：手動指定縱軸上緣（框頂＝該值且標數字）。
    %   空 = 自動（底部 500 貼軸、5 個等距刻度、頂端留 0.4 格）。志鵬的 l_hat 最大只到
    %   745，自動值 940 上方空一大截 -> 傳 800 壓下來。長飛最大 874.25，不可壓到 800。
    if nargin < 8, YTOP = []; end
    % [ADDED 2026-08-28 使用者要求] NOREF：不畫全格點基準虛線的版本。
    %   false（預設）＝原圖；true ＝拿掉兩條 dashed，檔名加 _noref（原圖不覆蓋）。
    if nargin < 9 || isempty(NOREF), NOREF = false; end

    % [MODIFIED 2026-08-30] 本檔由 temp_code/scripts/ 搬進 paper_fig_plot/plot/（產論文圖的腳本
    %   一律住 plot/，見 figure-output.md）。路徑改成相對自身推導，不再寫死絕對路徑：
    %   plot/ -> paper_fig_plot/ -> figures/ -> main/
    MAIN = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    addpath(fileparts(mfilename('fullpath')));   % 同夾的 axsh_sweep 等
    S = axsh_sweep(force, MODEL, GEOM, VARIANT, NMAX, Rum);

    figdir = fullfile(MAIN,'figures','paper_fig','Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    msuf = '';
    if ~strcmp(S.MODEL,'long2016_hexapole_halfcut')
        msuf = ['_' regexprep(S.MODEL,'_.*$','')];
    end
    nsuf = '';   if NOREF, nsuf = '_noref'; end
    out = fullfile(figdir, sprintf('ell_vs_npts_axsh_R%d%s%s.png', S.Rum, msuf, nsuf));

    e1 = 100*(S.ell1 - S.ref1(1))/S.ref1(1);
    e2 = 100*(S.ell2 - S.ref2(1))/S.ref2(1);
    fprintf('全格點 %d 點：single %.2f um  eighteen %.2f um\n', S.nref, S.ref1(1), S.ref2(1));
    fprintf('相對誤差 single  ：N=%d %+0.4f%%  ->  N=%d %+0.4f%%（穩態）\n', ...
            S.N(1), e1(1), S.N(end), e1(end));
    fprintf('相對誤差 eighteen：N=%d %+0.4f%%  ->  N=%d %+0.4f%%（穩態）\n', ...
            S.N(1), e2(1), S.N(end), e2(end));

    draw([0 S.N], [500 S.ell1], [500 S.ell2], S.ref1(1), S.ref2(1), ...
         '$\mathbf{\hat{\ell}\;(micro\;meter)}$', XCUT, 500, YTOP, out, NOREF);
end

% ============================================================================
function draw(N, v1, v2, r1, r2, ylab, XMAX, Y0, YTOP, out, NOREF)
    % [MODIFIED 2026-09-01 套繪圖規則 1/2/3] 字級拆成刻度 FS / 軸標題 FSLAB / 圖例 FSLEG；
    %   線寬與 marker 隨刻度字級 36->60（x1.67）同量級放大，框線 2.5->5.0。
    FS = 60;  FSLAB = 36;  FSLEG = 45;  LW = 5.0;  MS = 12;  LWBOX = 5.0;
    c1 = [0.05 0.10 0.95];   c2 = [0.85 0.10 0.10];
    % [MODIFIED 2026-09-01 規則 6] 實體尺寸等邊：畫布用**英吋**給（像素會被 MATLAB 靜默縮小）。
    %   邊長由「圖例放得下」反推，與 axsh_gain / axsh_kfro 一致。
    CANV = 14.5;
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes(fig);   hold(ax,'on');

    % 線：帶到第一個超出 XMAX 的設計，讓座標框把線段裁在 x = XMAX（不補假點）
    sel  = N <= XMAX;
    selL = sel;   k = find(~sel, 1);   if ~isempty(k), selL(k) = true; end
    plot(ax, N(selL), v1(selL), '-', 'Color',c1, 'LineWidth',LW, 'Clipping','on');
    plot(ax, N(selL), v2(selL), '-', 'Color',c2, 'LineWidth',LW, 'Clipping','on');
    % marker：只在框內的真實設計上，關掉裁切讓 x=0 的 marker 不被切一半
    plot(ax, N(sel), v1(sel), 'o', 'Color',c1, 'MarkerSize',MS, ...
         'MarkerFaceColor',c1, 'LineStyle','none', 'Clipping','off');
    plot(ax, N(sel), v2(sel), 's', 'Color',c2, 'MarkerSize',MS, 'LineWidth',LW, ...
         'MarkerFaceColor','w', 'LineStyle','none', 'Clipping','off');
    % 全格點基準（不進圖例）。
    % ⚠ 兩條只差 0.29 um（874.25 vs 873.96，相對 0.03%），在 500~940 的縱軸上是**同一條線**
    %   —— 後畫的會完全蓋住先畫的。使用者要看到藍色，所以**紅先藍後**（藍在上）。
    if ~NOREF
        yline(ax, r2, '--', 'Color',c2, 'LineWidth',LW, 'Alpha',1, 'HandleVisibility','off');
        yline(ax, r1, '--', 'Color',c1, 'LineWidth',LW, 'Alpha',1, 'HandleVisibility','off');
    end
    % 圖例用的代理 handle（線 + marker 一起顯示）
    p1 = plot(ax, NaN, NaN, '-o', 'Color',c1, 'LineWidth',LW, 'MarkerSize',MS, 'MarkerFaceColor',c1);
    p2 = plot(ax, NaN, NaN, '-s', 'Color',c2, 'LineWidth',LW, 'MarkerSize',MS, 'MarkerFaceColor','w');

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX,'TickLength',[.02 .02]);
    ax.Toolbar.Visible = 'off';

    xlim(ax, [0 XMAX]);   set(ax, 'XTick', xticks_inner(XMAX));
    if isempty(YTOP)
        vmax = max([v1(selL) v2(selL)]);
        if ~NOREF, vmax = max([vmax r1 r2]); end     % 虛線拿掉後不再撐高縱軸
        [yl, tk] = ylim_flush(Y0, vmax);
    else
        [yl, tk] = ylim_top(Y0, YTOP);
    end
    ylim(ax, yl);   set(ax,'YTick', tk);

    yr = ylim(ax);   yoff = yr(1) - 0.022*diff(yr);
    for xv = [0 XMAX]
        text(ax, xv, yoff, sprintf('%g',xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    xlabel(ax, '$\mathbf{Number\;of\;points}$', 'Interpreter','latex', 'FontSize',FSLAB);
    ylabel(ax, ylab, 'Interpreter','latex', 'FontSize',FSLAB);

    lg = legend(ax, [p1 p2], {'Single parameter', 'Eighteen parameters'}, ...
                'Interpreter','tex', 'Location','northoutside', 'NumColumns',2);
    lg.FontSize = FSLEG;   lg.FontWeight = 'bold';
    lg.Box = 'on';      lg.EdgeColor = 'k';   lg.LineWidth = LWBOX;   % [規則7] 圖例框與座標框同粗
    drawnow;
    axp = get(ax,'Position');   lgp = get(lg,'Position');
    lgw = lgp(3);   lgh = lgp(4);   GAPN = 0.022;
    newTop = 1 - lgh - GAPN - 0.006;
    axp(4) = newTop - axp(2);   set(ax,'Position',axp);
    if lgw < 0.70*axp(3)
        set(lg, 'Position', [axp(1) + (axp(3)-lgw)/2, newTop + GAPN, lgw, lgh]);
    else
        set(lg, 'Position', [axp(1), newTop + GAPN, axp(3), lgh]);
    end

    % [MODIFIED 2026-09-01 規則 6] exportgraphics 會裁掉畫布四周白邊 -> 正方形匯出後不等邊。
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
    print(fig, out, '-dpng', '-r200');            % 14.5 in x 200 dpi = 2900 px 見方
    fprintf('wrote %s\n', out);
end

% ============================================================================
function tk = xticks_inner(X)
% [規則 5 的水平軸例外，使用者拍板 2026-09-01]「水平軸維持四根 tick」
%   規則 5 要求奇數個，但 xlim=[0,50] 之下奇數且整數的解只有 s=5（9 根），
%   刻度字級 60 會撞成一團被 MATLAB 轉 45 度；奇數且疏的 s=12.5 是非整數。
%   => 水平軸改走**整數刻度、允許偶數根**；仍守「等間距，且與起點 0 / 終點 X
%      的間距也等於刻度間距」，即 s = X/(n+1)。X=50 -> s=10 -> 10/20/30/40。
    best = [];   bd = inf;
    for k = 0:4
        for c = [1 2 2.5 5]
            s = c*10^k;
            if abs(s - round(s)) > 1e-9, continue; end       % s 必須是整數
            if abs(X/s - round(X/s)) > 1e-9, continue; end   % 端點間距 = s
            n = round(X/s) - 1;
            if n < 3 || n > 5, continue; end
            if abs(n-4) < bd, bd = abs(n-4);   best = (1:n)*s; end
        end
    end
    if isempty(best), best = round(linspace(0,X,5));  best = best(2:end-1); end
    tk = best;
end

% ============================================================================
function [yl, tk] = ylim_top(y0, ytop)
% [規則 4/5] 框底 y0 與框頂 ytop **都不標數字**；中間奇數根等距刻度，
%   且與兩端的間距也相等 => [y0,ytop] 平分 n+1 段、刻度 = y0 + (1:n)*s。
%   [MODIFIED 2026-09-01 使用者指定] **固定四根**（與 ylim_flush 一致的例外）：
%   500~800、n=4 -> s=(800-500)/5=60 -> 560/620/680/740，框範圍完全不動、無空白代價。
    for n = 4
        s  = (ytop - y0)/(n+1);
        tk = y0 + (1:n)*s;
        if any(abs(tk*10 - round(tk*10)) > 1e-6), continue; end   % 需為整數或一位小數
        yl = [y0, ytop];   return
    end
    n = 5;   s = (ytop - y0)/(n+1);   tk = y0 + (1:n)*s;   yl = [y0, ytop];
end

% ============================================================================
function [yl, tk] = ylim_flush(y0, ymax)
% [規則 4/5] 框底 y0 與框頂 T **都不標數字**；內部奇數根等距，與兩端間距亦相等
%   => T = y0 + (n+1)*s、刻度 = y0 + (1:n)*s。s 取 nice 值，須讓 T 蓋過 ymax
%   並留約 8% 裕度；候選中取**框頂最低**者（避免上方留一大截空白）。
    %   ⚠ 裕度不可設太大：0.08 會讓 (n,s)=(3,100)、T=900 的解差 4 um 被剔掉，
    %     退而選 (5,100)、T=1100 -> 資料只到 874，上方空掉 226（填充率 66%）。
    %   [MODIFIED 2026-09-01 使用者指定] **本函式的縱軸固定四根**（規則 5 的另一個例外，
    %     與水平軸同樣是「整數刻度優先於奇數根」）：n=4 -> s=100 -> 600/700/800/900、
    %     框頂 1000。奇數解 n=3 雖然框頂 900 較緊（填充 94%），但使用者要四根。
    R = max(ymax - y0, eps);   need = ymax + 0.02*R;
    bestT = inf;   yl = [];   tk = [];
    for n = 4
        for k = -3:5
            for c = [1 2 2.5 3 4 5 10]
                s = c*10^k;
                T = y0 + (n+1)*s;
                if T < need, continue; end
                t = y0 + (1:n)*s;
                if any(abs(t*10 - round(t*10)) > 1e-6), break; end
                if T < bestT, bestT = T;   tk = t;   yl = [y0, T]; end
                break
            end
        end
    end
    if isempty(tk)
        n = 5;   s = (need - y0)/(n+1);   tk = y0 + (1:n)*s;   yl = [y0, y0+(n+1)*s];
    end
end
