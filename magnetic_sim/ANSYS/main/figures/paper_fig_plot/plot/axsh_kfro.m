function axsh_kfro(force, MODEL, GEOM, VARIANT, NMAX, Rum, XCUT, NOLINE)
%AXSH_KFRO  Frobenius error of K_I_bar against the full-grid K_I_bar.
% =========================================================================
%   One script, one figure (repo rule). Reads the cache written by
%   axsh_sweep.m and draws nothing else.
%
%       y = ||K_bar(N) - K_bar(full grid)||_F / ||K_bar(full grid)||_F * 100  [%]
%
%   Note this is NOT the step-to-step Frobenius change that was retired in
%   2026-08-14 -- it is measured against the full-grid answer, so it says how
%   far the reduced design actually is, not merely how much the last rung
%   moved. There is no dashed reference line: the reference IS zero.
%
%   Style preset (1) bold-framed, LINEAR x, y from 0 (the floor is physical).
%
%   [ADDED 2026-08-29 使用者指定] NOLINE：true = **只畫點、不用線連**（marker only），
%   輸出檔名加 '_pts' 後綴、不覆蓋有線版。資料、軸範圍、刻度、圖例全部相同，
%   照 plot_ell_gain_2panel.m 的同名參數辦。
%
%   Output -> figures/paper_fig/Section2_E/kfro_vs_npts_axsh_R<Rum>[_<model>][_pts].png
% =========================================================================
    clc;
    if nargin < 1 || isempty(force),   force   = false;                       end
    if nargin < 2 || isempty(MODEL),   MODEL   = 'long2016_hexapole_halfcut'; end
    if nargin < 3 || isempty(GEOM),    GEOM    = 'tip40um';                   end
    if nargin < 4,                     VARIANT = '';                          end
    if nargin < 5 || isempty(NMAX),    NMAX    = 500;                         end
    if nargin < 6 || isempty(Rum),     Rum     = 150;                         end
    if nargin < 7 || isempty(XCUT),    XCUT    = 50;                          end
    if nargin < 8 || isempty(NOLINE),  NOLINE  = false;                       end

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
    psuf = '';   if NOLINE, psuf = '_pts'; end
    out = fullfile(figdir, sprintf('kfro_vs_npts_axsh_R%d%s%s.png', S.Rum, msuf, psuf));

    fprintf('dK/K  single  : %.4f%% ~ %.4f%%（末端 %.4f%%）\n', ...
            min(S.fro1), max(S.fro1), S.fro1(end));
    fprintf('dK/K  eighteen: %.4f%% ~ %.4f%%（末端 %.4f%%）\n', ...
            min(S.fro2), max(S.fro2), S.fro2(end));

    draw(S.N, S.fro1, S.fro2, ...
         '$\mathbf{\|\bar{K}_{I}-\bar{K}_{I}^{full}\|_{F}/\|\bar{K}_{I}^{full}\|_{F}\;(\%)}$', ...
         XCUT, out, NOLINE);
end

% ============================================================================
function draw(N, v1, v2, ylab, XMAX, out, NOLINE)
    % [MODIFIED 2026-09-01] 依 LAB406 那篇 paper 實測的版式定字級：內文 10 pt、圖說 8 pt、
    %   圖內刻度數字 7.0~7.8 pt、圖內軸標題 10.0~10.2 pt。目標擺放寬度 = 3.0 in
    %   （該篇單欄圖實測 2.54/2.92/2.99 in，並非塞滿 3.5 in 欄寬）。
    %   畫布寬約 10.6 in -> 換算倍率 10.6/3.0 = 3.53：
    %     刻度數字 8 pt  -> 28      軸標題 10 pt -> 35      圖例 8 pt -> 28
    %   原本 FS 同時餵刻度與軸標題，這裡拆成 FS / FSLAB 兩個。
    % [MODIFIED 2026-09-01 使用者拍板的新繪圖規則，見 .claude/rules/figure-style.md]
    %   規則1 刻度數字 = 60；規則2 圖例 = 45 且圖例線段樣本與資料線一併加粗；
    %   規則3 框線加粗；規則6 畫布等邊。LW/MS 依刻度 36->60 (x1.67) 的同量級放大。
    FS = 60;  FSLAB = 35;  FSLEG = 45;  LW = 5.0;  MS = 12;  LWBOX = 5.0;
    c1 = [0.05 0.10 0.95];   c2 = [0.85 0.10 0.10];
    % [MODIFIED 2026-09-01 使用者要求] 畫布改**等邊**，且尺寸以英吋給（不再用像素）。
    %   邊長 12.5 in 是由「圖例放得下」反推的下限：圖例在 FSLEG=37 時自然寬度約
    %   9.3 in，加上 y 軸標題與刻度數字佔掉的左緣約 1.3 in = 10.6 in；畫布若只有
    %   11.4 in，legend 對齊碼會把座標框推到畫布右緣、圖例左右被裁掉（FSLEG=40 時
    %   實測就是這樣爆版）。12.5 in 留了餘裕。
    %   FSLEG=40 時圖例自然寬約 11.3 in，加左緣 1.7 in = 13.0 in -> 12.5 in 會被裁
    %   （實測：圖例右緣、右框線、末端刻度 50 全部切掉）。13.5 in 留 0.5 in 餘裕。
    CANV = 14.5;                                   % 畫布邊長 [in]（正方形）
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes(fig);   hold(ax,'on');

    % 線帶到第一個超出 XMAX 的設計、由座標框裁在 x = XMAX（不補假點）；marker 只畫
    % 框內的真實設計。原本整條 4801 點都畫且 Clipping off -> 尾端跑到框外面去了。
    sel  = N <= XMAX;
    selL = sel;   k = find(~sel, 1);   if ~isempty(k), selL(k) = true; end
    if ~NOLINE                                        % NOLINE 時整條連線都不畫
        plot(ax, N(selL), v1(selL), '-', 'Color',c1, 'LineWidth',LW, 'Clipping','on');
        plot(ax, N(selL), v2(selL), '-', 'Color',c2, 'LineWidth',LW, 'Clipping','on');
    end
    plot(ax, N(sel), v1(sel), 'o', 'Color',c1, 'MarkerSize',MS, ...
         'MarkerFaceColor',c1, 'LineStyle','none', 'Clipping','off');
    plot(ax, N(sel), v2(sel), 's', 'Color',c2, 'MarkerSize',MS, 'LineWidth',LW, ...
         'MarkerFaceColor','w', 'LineStyle','none', 'Clipping','off');
    if NOLINE, LG1 = 'o';  LG2 = 's';  else, LG1 = '-o';  LG2 = '-s';  end
    h1 = plot(ax, NaN, NaN, LG1, 'Color',c1, 'LineWidth',LW, 'MarkerSize',MS, 'MarkerFaceColor',c1);
    h2 = plot(ax, NaN, NaN, LG2, 'Color',c2, 'LineWidth',LW, 'MarkerSize',MS, 'MarkerFaceColor','w');

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX,'TickLength',[.02 .02]);
    ax.Toolbar.Visible = 'off';

    % [MODIFIED 2026-08-28 使用者要求] 水平軸自**第一個設計**（7 點）起，不再從 0 起
    %   —— 這條曲線在 x=0 沒有錨點（不像 ell/gain 有初值），從 0 起只是左邊留一段空白。
    %   照 figure-style「曲線圖首末資料點貼齊框邊」：xlim = [第一點, XMAX]。
    XMIN = N(find(sel,1));
    xlim(ax, [XMIN XMAX]);   set(ax, 'XTick', xticks_inner(XMIN, XMAX));
    % 縱軸自 0 起：內部刻度 (1:n)*s 等距、起點與終點都不標（figure-style 慣例 #7）
    vmax = max([v1(selL) v2(selL)]);
    [yl, tk] = ylim_from_zero(vmax);
    ylim(ax, yl);   set(ax,'YTick', tk);
    fprintf('  縱軸 ylim=[%g %g]  YTick=%s  (%d 根，填充 %.0f%%)\n', ...
            yl(1), yl(2), mat2str(tk), numel(tk), 100*vmax/yl(2));

    yr = ylim(ax);   yoff = yr(1) - 0.022*diff(yr);
    for xv = [XMIN XMAX]
        text(ax, xv, yoff, sprintf('%g',xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    xlabel(ax, '$\mathbf{Number\;of\;points}$', 'Interpreter','latex', 'FontSize',FSLAB);
    ylabel(ax, ylab, 'Interpreter','latex', 'FontSize',FSLAB);

    lg = legend(ax, [h1 h2], {'Single parameter', 'Eighteen parameters'}, ...
                'Interpreter','tex', 'Location','northoutside', 'NumColumns',2);
    lg.FontSize = FSLEG;   lg.FontWeight = 'bold';
    lg.Box = 'on';      lg.EdgeColor = 'k';   lg.LineWidth = LWBOX;   % 圖例框與座標框同粗
    % [FIXED 2026-09-02 使用者回報「圖例的藍點太靠近圖例框」] ItemTokenSize（圖例 icon 區
    %   寬度）預設固定 30 pt，**不會隨 FontSize 放大**。字級提到 45 pt 後 icon 區相對變窄，
    %   marker 被擠到框線邊上：實測藍點距框左內緣僅 18 px（_pts）/ 15 px（有線版），而框線
    %   本身就 14 px 粗。marker 畫在 token 的中心，故 token 加寬 W 會讓 marker 右移 W/2。
    lg.ItemTokenSize = [55 25];
    drawnow;
    axp = get(ax,'Position');   lgp = get(lg,'Position');
    lgw = lgp(3);   lgh = lgp(4);   GAPN = 0.022;
    newTop = 1 - lgh - GAPN - 0.006;
    axp(4) = newTop - axp(2);
    % [MODIFIED 2026-09-01 使用者要求] 圖例框與座標框左右緣切齊。
    %   ⚠ 不能靠「把 legend 的 Position 寬度設成座標框寬」達成 —— MATLAB **不會把
    %     legend 縮到比它的內容還窄**，只會把框撐回自然寬度並對給定矩形置中。實測
    %     FSLEG=35 時圖例框 1970 px、座標框 1815 px，左右各溢出 77/78 px（對稱 =
    %     置中的證據）。字級要維持不變就只剩一條路：**反過來把座標框加寬去配合圖例**。
    %   作法：座標框寬度 = 圖例自然寬度 lgw（右緣不得超出畫布，故取 min 夾住）。
    %   若圖例窄於座標框（小字級時）則沿用舊行為：圖例保持自然寬度、對座標框置中。
    % [FIXED 2026-09-02 使用者回報「水平軸被截到」] 右緣保留由 0.005（≈14 px）加到 0.032
    %   （≈93 px）。端點數字是**以座標框邊緣為中心**畫的，「50」在 FS=60 時半寬實測 77 px
    %   （量 _pts 與 gain 兩張未被裁的圖，框右緣到數字右緣都是 77 px）。原本只留 14 px ->
    %   長飛有線版餘裕 49 px、志鵬 8 px，數字右半被畫布切掉；志鵬連圖例框右緣也一起切掉。
    %   圖例自然寬 ~2480 px 在 2900 px 畫布內本來就放得下，所以不必動畫布或字級，
    %   只要別把座標框推到畫布邊上；圖例若比座標框寬，MATLAB 會自動置中、往左讓出空間。
    %   0.032（93 px）剛好不裁，但數字右緣離畫布只剩 9 px、看起來像貼在邊上；
    %   改 0.040（116 px）讓數字外側還有 ~32 px，與未被裁的 _pts 版（39 px）相當。
    RPAD = 0.040;
    if lgw > axp(3)
        axp(3) = min(lgw, 1 - axp(1) - RPAD);
        set(ax, 'Position', axp);
        % [FIXED 2026-09-02] 圖例用**自然寬度**並夾住右緣，不再把寬度設成 axp(3) 讓
        %   MATLAB 自己置中 —— 那會在圖例比座標框寬時對稱溢出、右邊被畫布切掉
        %   （加大 ItemTokenSize 後志鵬 kfro 立刻重現此問題，右餘裕 0）。
        %   左緣仍以座標框為準；放不下時才整個往左移到剛好放得下。
        lgx = max(0.015, min(axp(1), 1 - 0.015 - lgw));
        set(lg, 'Position', [lgx, newTop + GAPN, lgw, lgh]);
    else
        set(ax, 'Position', axp);
        set(lg, 'Position', [axp(1) + (axp(3)-lgw)/2, newTop + GAPN, lgw, lgh]);
    end

    % print 而非 exportgraphics：後者**會裁掉畫布四周的白邊**，正方形畫布匯出後就不再
    %   是正方形（plot_ell_gain_2panel 的註解記了實測值 0.96 / 1.11）。要保證等邊就必須
    %   用 print 搭配明確的 PaperPosition，把整張畫布原封輸出。
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
    print(fig, out, '-dpng', '-r200');             % 12.5 in x 200 dpi = 2500 px 見方
    fprintf('wrote %s\n', out);
end

% ============================================================================
function tk = xticks_inner(LO, HI)
% 內部整數刻度 3~5 個（優先奇數），不碰兩端（端點另以 text 標數字）。
%   [MODIFIED 2026-08-28] 支援非 0 起點：刻度只取落在 (LO, HI) 內部者。
    if nargin < 2, HI = LO;  LO = 0; end
    % [MODIFIED 2026-08-28] 端點保留帶由 4% 加大到 10%：起點改成 7 之後，s=10 的第一個
    %   刻度「10」離端點「7」只有 7% 圖寬，數字會黏在一起；且那組是 4 個刻度（偶數，
    %   違反 figure-style 的奇數規則）。10% 保留帶會把它剔掉 -> 20/30/40（3 個、奇數）。
    pad = 0.10*(HI-LO);
    best = [];   bscore = inf;
    for k = 0:4
        for c = [1 2 2.5 5 10]
            s = c*10^k;
            if abs(s - round(s)) > 1e-9, continue; end
            t = ceil((LO+pad)/s)*s : s : floor((HI-pad)/s)*s;
            n = numel(t);
            if n < 3 || n > 5, continue; end
            score = 10*(mod(n,2)==0) + abs(n-4);
            if score < bscore, bscore = score;  best = t; end
        end
    end
    if isempty(best), best = round(linspace(LO,HI,5));  best = best(2:end-1); end
    tk = best;
end

% ============================================================================
function [yl, tk] = ylim_from_zero(maxv)
% 自 0 起：把 [0, T] 平分 n+1 段，刻度 = (1:n)*s（n 取奇數 3），起點與終點都不標。
%   [MODIFIED 2026-08-28 使用者：上面留太多、但縱軸維持三根 tick] 步長改照
%   figure-style 慣例 #7 的原始寫法：**s 取 0.1 的倍數**（不再限制在 nice 十進位
%   1/2/2.5/5/10），在 [smin, 1.15*smin] 窗口內挑最漂亮者（整數 > 0.5 的倍數 >
%   0.2 的倍數 > 其餘，同分取最小）。
%   舊版只認 nice 十進位 -> 志鵬 maxv=2.10 時 smin=0.567 被迫跳到 s=1（上緣 4.0、
%   填充率僅 52%）；新版取 s=0.6（刻度 0.6/1.2/1.8、上緣 2.4、填充率 87%）。
    n = 3;
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
