function axsh_gain(force, MODEL, GEOM, VARIANT, NMAX, Rum, XCUT, NOREF)
%AXSH_GAIN  g_I_hat vs point count along the actuator-axis shell ladder.
% =========================================================================
%   One script, one figure (repo rule). Reads the cache written by
%   axsh_sweep.m and draws nothing else.
%
%       solid  : g_I_hat from the reduced design  (single blue / eighteen red)
%       dashed : g_I_hat from the FULL grid inside R, same colour (NOT in the
%                legend; the legend keeps only the two model names).
%
%   Curve starts at x = 0 with value 0 (g_I has no initial guess).
%   Right edge: the line is carried through the first rung past XCUT and the
%   axes clip it there, so it meets the frame with no gap and no invented
%   point. See axsh_ell.m for the full note.
%
%   Output -> figures/paper_fig/Section2_E/gain_vs_npts_axsh_R<Rum>[_<model>][_noref].png
% =========================================================================
    clc;
    if nargin < 1 || isempty(force),   force   = false;                       end
    if nargin < 2 || isempty(MODEL),   MODEL   = 'long2016_hexapole_halfcut'; end
    if nargin < 3 || isempty(GEOM),    GEOM    = 'tip40um';                   end
    if nargin < 4,                     VARIANT = '';                          end
    if nargin < 5 || isempty(NMAX),    NMAX    = 500;                         end
    if nargin < 6 || isempty(Rum),     Rum     = 150;                         end
    if nargin < 7 || isempty(XCUT),    XCUT    = 50;                          end
    % [ADDED 2026-08-28 使用者要求] NOREF：不畫全格點基準虛線的版本。
    %   false（預設）＝原圖；true ＝拿掉兩條 dashed，檔名加 _noref（原圖不覆蓋）。
    if nargin < 8 || isempty(NOREF),   NOREF   = false;                       end

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
    out = fullfile(figdir, sprintf('gain_vs_npts_axsh_R%d%s%s.png', S.Rum, msuf, nsuf));

    e1 = 100*(S.gI1 - S.ref1(2))/S.ref1(2);
    e2 = 100*(S.gI2 - S.ref2(2))/S.ref2(2);
    fprintf('全格點 %d 點：single %.4f mT/A  eighteen %.4f mT/A\n', S.nref, S.ref1(2), S.ref2(2));
    fprintf('相對誤差 single  ：N=%d %+0.4f%%  ->  N=%d %+0.4f%%（穩態）\n', ...
            S.N(1), e1(1), S.N(end), e1(end));
    fprintf('相對誤差 eighteen：N=%d %+0.4f%%  ->  N=%d %+0.4f%%（穩態）\n', ...
            S.N(1), e2(1), S.N(end), e2(end));

    draw([0 S.N], [0 S.gI1], [0 S.gI2], S.ref1(2), S.ref2(2), ...
         '$\mathbf{{}^{B}\hat{g}_{I}\;(mT/A)}$', XCUT, 0, out, NOREF);
end

% ============================================================================
function draw(N, v1, v2, r1, r2, ylab, XMAX, Y0, out, NOREF)
    % [MODIFIED 2026-09-01 使用者指定] 原本 FS 同時餵刻度數字與軸標題、圖例寫死 24。
    %   拆成三個：刻度 FS / 軸標題 FSLAB / 圖例 FSLEG，才能分別調。
    %   使用者要的值：刻度 50、圖例 45。軸標題維持原本的 36（未指定，不動）。
    % [MODIFIED 2026-09-01 使用者拍板的新繪圖規則，見 .claude/rules/figure-style.md]
    %   規則1 刻度數字 = 60；規則2 圖例 = 45 且圖例線段樣本與資料線一併加粗；
    %   規則3 框線加粗；規則6 畫布等邊。LW/MS 依刻度 36->60 (x1.67) 的同量級放大。
    FS = 60;  FSLAB = 36;  FSLEG = 45;  LW = 5.0;  MS = 12;  LWBOX = 5.0;
    c1 = [0.05 0.10 0.95];   c2 = [0.85 0.10 0.10];
    % [MODIFIED 2026-09-01] 畫布改**等邊**且以英吋給。邊長 14.5 in 是由圖例反推：
    %   45 pt 的圖例自然寬約 12.7 in，原本 11.7 in 的畫布放不下會被裁掉。
    CANV = 14.5;                                   % 畫布邊長 [in]（正方形）
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes(fig);   hold(ax,'on');

    sel  = N <= XMAX;
    selL = sel;   k = find(~sel, 1);   if ~isempty(k), selL(k) = true; end
    plot(ax, N(selL), v1(selL), '-', 'Color',c1, 'LineWidth',LW, 'Clipping','on');
    plot(ax, N(selL), v2(selL), '-', 'Color',c2, 'LineWidth',LW, 'Clipping','on');
    plot(ax, N(sel), v1(sel), 'o', 'Color',c1, 'MarkerSize',MS, ...
         'MarkerFaceColor',c1, 'LineStyle','none', 'Clipping','off');
    plot(ax, N(sel), v2(sel), 's', 'Color',c2, 'MarkerSize',MS, 'LineWidth',LW, ...
         'MarkerFaceColor','w', 'LineStyle','none', 'Clipping','off');
    if ~NOREF
        yline(ax, r1, '--', 'Color',c1, 'LineWidth',LW, 'Alpha',1, 'HandleVisibility','off');
        yline(ax, r2, '--', 'Color',c2, 'LineWidth',LW, 'Alpha',1, 'HandleVisibility','off');
    end
    p1 = plot(ax, NaN, NaN, '-o', 'Color',c1, 'LineWidth',LW, 'MarkerSize',MS, 'MarkerFaceColor',c1);
    p2 = plot(ax, NaN, NaN, '-s', 'Color',c2, 'LineWidth',LW, 'MarkerSize',MS, 'MarkerFaceColor','w');

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX,'TickLength',[.02 .02]);
    ax.Toolbar.Visible = 'off';

    xlim(ax, [0 XMAX]);   set(ax, 'XTick', xticks_inner(XMAX));
    vmax = max([v1(selL) v2(selL)]);
    if ~NOREF, vmax = max([vmax r1 r2]); end          % 虛線拿掉後不再撐高縱軸
    [yl, tk] = ylim_flush(Y0, vmax);
    ylim(ax, yl);   set(ax,'YTick', tk);
    fprintf('  縱軸 ylim=[%g %g]  YTick=%s  (%d 根，填充 %.0f%%)\n', ...
            yl(1), yl(2), mat2str(tk), numel(tk), 100*vmax/yl(2));

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
    % [FIXED 2026-09-02 使用者回報「圖例的藍點太靠近圖例框」] 同 axsh_kfro：ItemTokenSize
    %   預設 30 pt 不隨 FontSize 放大，字級 45 pt 下 marker 會貼到框線。marker 在 token
    %   中心，token 加寬 W -> marker 右移 W/2。兩支腳本用同一個值以保持八張圖一致。
    lg.ItemTokenSize = [55 25];
    drawnow;
    axp = get(ax,'Position');   lgp = get(lg,'Position');
    lgw = lgp(3);   lgh = lgp(4);   GAPN = 0.022;
    newTop = 1 - lgh - GAPN - 0.006;
    axp(4) = newTop - axp(2);   set(ax,'Position',axp);
    if lgw < 0.70*axp(3)
        set(lg, 'Position', [axp(1) + (axp(3)-lgw)/2, newTop + GAPN, lgw, lgh]);
    else
        % [FIXED 2026-09-02] 同 axsh_kfro：給自然寬度並夾住右緣，避免圖例比座標框寬時
        %   被 MATLAB 對稱置中而溢出畫布。
        lgx = max(0.015, min(axp(1), 1 - 0.015 - lgw));
        set(lg, 'Position', [lgx, newTop + GAPN, lgw, lgh]);
    end

    % print 而非 exportgraphics：後者會裁掉畫布四周白邊，正方形畫布匯出後就不是正方形。
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
    print(fig, out, '-dpng', '-r200');             % 14.5 in x 200 dpi = 2900 px 見方
    fprintf('wrote %s\n', out);
end

% ============================================================================
function tk = xticks_inner(X)
% [規則 5 的水平軸例外，使用者拍板 2026-09-01]
%   規則 5 要求「奇數個 tick」，但 xlim=[0,50] 之下唯一的奇數整數解是 s=5（9 根），
%   在刻度字級 60 會撞成一團被 MATLAB 轉 45 度（實測）；奇數且疏的解 s=12.5 會得到
%   12.5 / 25 / 37.5 這種非整數刻度。使用者選擇**水平軸改走整數刻度、允許偶數根**。
%   仍守的部分：等間距，且與起點 0 / 終點 X 的間距也等於刻度間距（s = X/(n+1)）。
%   => 找「X/s 為整數」的 nice 整數 s，內部根數 n = X/s - 1 落在 3~5，取最接近 4 者。
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
function [yl, tk] = ylim_flush(y0, ymax)
% [規則 4] 縱軸**起點與終點都不標**；內部刻度等間距，且與兩端的間距也相等。
%   => [y0, T] 平分 n+1 段：刻度 = y0 + (1:n)*s、T = y0 + (n+1)*s。
%   s 需讓 T >= 1.08*ymax（上緣留一點裕度），且刻度數字是整數或一位小數。
%   y0 與 T 都不進 tk -> MATLAB 自然不會標它們。
%
% [縱軸 tick 數的例外，使用者拍板 2026-09-01]
%   規則 5 要求「奇數個 tick」，但使用者對本組四張圖（長飛 / 志鵬 × ref / noref）
%   明確指定**縱軸四根**。故 n 固定為 4，不再試 [5 3]。
%   實測：長飛 ymax=10.66 -> s=2.5（2.5/5/7.5/10、上緣 12.5、填充 85%）；
%         志鵬 ymax=37.56 -> s=10 （10/20/30/40、上緣 50、填充 75%）。
    % [MODIFIED 2026-09-01] 步長候選加入 6/7/8（原本 [1 2 2.5 3 4 5 10] 從 5 直接跳到 10）。
    %   志鵬 R150 的 ymax=37.56 需要 s>=6.76，舊清單只能取 s=10 -> 上緣 60、填充率僅 63%，
    %   縱軸上方空掉一大塊。加入後取 s=7（刻度 7/14/21/28/35、上緣 42、填充率 89%）。
    %   長飛 R150（ymax=10.66、smin=1.92）仍取 s=2，輸出不變。
    for n = 4
        smin = (1.08*ymax - y0)/(n+1);
        for k = -3:5
            for c = [1 1.5 2 2.5 3 4 5 6 7 8 10]
                s = c*10^k;
                if s < smin, continue; end
                tk = y0 + (1:n)*s;
                if any(abs(tk*10 - round(tk*10)) > 1e-6), continue; end
                yl = [y0, y0 + (n+1)*s];   return
            end
        end
    end
    s = (ymax - y0)/4;   tk = y0 + (1:4)*s;   yl = [y0, y0 + 5*s];   % fallback，同樣 4 根
end
