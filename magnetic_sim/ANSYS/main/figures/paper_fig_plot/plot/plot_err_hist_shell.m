function plot_err_hist_shell(USE_BIAS, Rum, Rsplit, SPLIT)
% plot_err_hist_shell -- 「用 N_c 減量校正的模型」對範圍內全部真實格點的殘差直方圖，
%                        依格點半徑分內外兩層疊圖
% =========================================================================
%   ① 校正：取 plot_full_vs_conv_vs_R.m 為該 R 定出的**收斂設計**（讀其快取
%      full_vs_conv_vs_R_maxwell.mat 的 tri_c1/tri_c2；判準＝l_hat 穩定 ∧ g_I 穩定
%      ∧ K_I 符合物理，各持續 10 步），用那 N_c 個等測度網格點擬合 → (l_hat, e, G)。
%   ② 評估：拿 **R 以內全部真實 .fld 格點**（不是校正用的內插查詢點），用 ① 的
%      (l_hat, e, G) 預測，逐「格點 × 激發」算殘差大小 ||b_FEM - b_model|| [mT]。
%      ⚠ G 用校正階段解出來的、不在評估集上重解 → 真正的外推驗證。
%   ③ 分層：依格點半徑 r 切成 r <= Rsplit 與 r > Rsplit 兩組，疊在同一張直方圖。
%      **各組各自正規化成百分比**（兩組樣本數差很多：R=300 內共 14082 點，其中
%      r<=150 只有 1771 點）—— 比較的是分布形狀，不是絕對數量。
%
%   直方圖規則（figure-style「分布 / 疊圖直方圖」）：**固定 bin 寬 BINW=2.8 uT**
%   （2026-08-20 由 nb=180 改來，跨圖共用一把尺）、兩組**共用 edges**、FaceAlpha。
%   [MODIFIED 2026-08-17 使用者拍板] **不再畫 mean 虛線**，圖例也隨之只列系列、不列統計值。
%   mean 仍算出來並印在 console 供追溯。
%   風格①粗體框圖：box on、grid off、刻度 FS28 粗體、軸標題 FS36 LaTeX \mathbf、
%   x 起訖只標數字不畫 tick、y 取 3 個內縮 tick。
%
%   配色：**沿用 err_hist_overlay 那組**（使用者指定）——藍 [0.10 0.35 1.00] / 紅
%   [0.85 0.10 0.10]、FaceAlpha 0.60、細黑邊 0.3；mean 虛線用中性色（黑 / 深綠），
%   中性色壓在兩組長條上都看得清（house 慣例，不隨長條色變）。
%   ⚠ 這裡的藍紅是「半徑分層」，不是 single/eighteen —— 本圖只有一個模型。
%
%   SPLIT=false：不分層，整個 R 以內的殘差畫成**一張藍色直方圖**（給「直接看整體分布」用）。
%
%   用法：plot_err_hist_shell                       → single、R=300、切 150（分層）
%         plot_err_hist_shell(false,300,150,false)  → single、R=300、不分層（全藍）
%         plot_err_hist_shell(true)                 → eighteen
%   輸出 → figures/paper_fig/Section2_E/err_hist_{shell,conv}_maxwell_<model>_R<Rum>.png
% =========================================================================
    clc;
    if nargin < 1 || isempty(USE_BIAS), USE_BIAS = false; end
    if nargin < 2 || isempty(Rum),      Rum      = 300;   end   % 校正 + 評估半徑 [um]
    if nargin < 3 || isempty(Rsplit),   Rsplit   = 150;   end   % 內外分層半徑 [um]
    if nargin < 4 || isempty(SPLIT),    SPLIT    = true;  end
    % [ADDED 2026-08-14] SPLIT 三種模式：
    %   true    → 'overlay'：兩層**各自正規化**疊圖（比較分布形狀；樣本數差 7 倍）
    %   false   → 'single' ：不分層、單一顏色
    %   'stack' → 'stack'  ：**同一個分布**（總高度 = 不分層時的分布），內部依半徑
    %                        分成兩色堆疊（使用者拍板：維持一樣的分布、只是內部分色）
    MODE = 'overlay';
    if islogical(SPLIT) && ~SPLIT,                              MODE = 'single'; end
    if (ischar(SPLIT) || isstring(SPLIT)) && strcmpi(SPLIT,'stack'), MODE = 'stack';  end

    here   = fileparts(fileparts(mfilename('fullpath')));       % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    addpath(fullfile(CAL,'function'), fullfile(CAL,'common_path'), fullfile(CAL,'utils'));
    addpath(fullfile(CAL,'utils','long2016_hexapole_halfcut'));

    %% ---- ① 取該 R 的收斂設計 → 減量校正 -----------------------------------
    cfg = model_config('long2016_hexapole_halfcut','tip40um');
    MODEL_ = 'long2016_hexapole_halfcut';   VARIANT_ = cfg.default_variant;

    % ---- 讀 main.m 產的收斂設計校正結果（不再自己重跑階梯 + 校正）--------
    %   [MODIFIED 2026-08-23 使用者拍板] 校正與收斂判準已搬回 main.m
    %   （conv_design_ws / conv_design_sensor 只負責決定內插點位置與取場），
    %   繪圖端改成**接收 main 產完的結果** -> 圖與結果 PDF 保證出自同一次校正。
    %   ⚠ 該組合必須先跑過 main.m（GRID_NRPT='auto'）；找不到就報錯，不猜。
    %   ⚠ 同一組合可能有多顆 convN 檔（舊實驗留下的，例如 long2016 R150
    %     eighteen 就有 convN6/convN80/convN88）-> 只認 conv_auto==true 那顆。
    md_ = fullfile(CAL, 'data', MODEL_, '.mat');
    tg_ = 'single';   if USE_BIAS, tg_ = 'eighteen'; end
    dd_ = dir(fullfile(md_, sprintf('calib_current_%s_convN*_R%03d_%s.mat', ...
                                    VARIANT_, round(Rum), tg_)));
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
           numel(dd_), MODEL_, Rum, USE_BIAS);
    cal = load(fullfile(md_, dd_(1).name));
    e = cal.e;   l_hat = cal.l_hat;   G = cal.G;
    tri = cal.GRID_NRPT;   Nc = cal.npts;   gI = cal.gI_hat;
    fprintf('校正：R<=%d um、N_c=%d (%d,%d,%d)  l_hat=%.1f um  g_I=%.4f mT/A\n', ...
            Rum, Nc, tri(1), tri(2), tri(3), l_hat*1e6, gI);

    %% ---- ② 評估集：R 以內全部真實格點 -------------------------------------
    raw = extract_maxwell_data(cfg, 'all', cfg.default_variant);
    ad  = build_actuator_data(raw, cfg);
    [P0, B0, n0] = cfg.select_ball(ad, Rum*1e-6);

    Pc   = make_Pc(e, cfg.Pc_base);
    S    = build_S(l_hat, Pc, P0);
    res  = S*G - B0;                                  % 3n0 x 6
    err  = zeros(n0, 6);
    for j = 1:6
        rj = reshape(res(:,j), 3, n0);
        err(:,j) = sqrt(sum(rj.^2, 1)).';             % 每格點 x 激發 的殘差大小 [mT]
    end
    RMSPE = 100 * sqrt(sum(res(:).^2) / sum(B0(:).^2));

    %% ---- ③ 依半徑分內外兩層 ------------------------------------------------
    rr  = vecnorm(P0, 2, 2);                          % 旋轉不變，actuator frame 直接量
    in  = rr <= Rsplit*1e-6;
    e1  = reshape(err(in,:),  [], 1);                 % 內層（含 6 個激發）
    e2  = reshape(err(~in,:), [], 1);                 % 外層
    mu1 = mean(e1);   mu2 = mean(e2);
    fprintf('評估格點 %d（內層 %d / 外層 %d）；RMSPE=%.3f%%\n', n0, sum(in), sum(~in), RMSPE);
    fprintf('  內層 mean=%.4f mT  max=%.4f｜外層 mean=%.4f mT  max=%.4f\n', ...
            mu1, max(e1), mu2, max(e2));

    %% ---- 直方圖 ------------------------------------------------------------
    % 配色：**依「是否分層」決定**（使用者拍板 2026-08-14）
    %   分層圖（兩個半徑層疊圖）：固定 深藍 [0.05 0.10 0.95] / 紅 [0.85 0.10 0.10]，
    %     不隨取樣半徑變 —— 藍紅的色相差最大，兩層交疊處才分得清。
    %   不分層（單一系列）：R=150 一律深藍（兩個模型同色，方便並排比較）；
    %     其他 R（300…）**依模型換色**，避免 single / eighteen 兩張圖撞色：
    %     single 用紫 #7B52AB（沿用 plot_err_hist.m 的 pick_colors）、
    %     eighteen 用深青 #008C99（使用者拍板 2026-08-14：紫色已用過，要不同色）。
    % [MODIFIED 2026-08-20 使用者拍板] 分層圖配色對調：**內層 R<=150 一律紅色**
    %   （跟 R150 單獨圖一致），外層改用藍色。
    cOUT = [0.05 0.10 0.95];                          % 外層：深藍（全局統一的藍）
    % [MODIFIED 2026-08-20 使用者拍板] R150 單獨圖改用**紅色**（= overlay 的 eighteen 色），
    %   因為它跟 overlay 的紅色是同一組資料。分層圖仍用藍(內)/紅(外)。
    % ⚠ [FIXED 2026-08-20] 紅色只給「**eighteen** 的 R150 單獨圖」跟分層圖的內層：
    %   前者跟 overlay 的紅色是同一組資料。single(fix) 的 R150 圖**不是**那組資料，
    %   仍用深藍，否則兩張不同模型的圖會擞成同色。
    if ~strcmp(MODE,'single') || (Rum == 150 && USE_BIAS)
        cIN = [0.85 0.10 0.10];                       % 內層 / eighteen R150 單獨圖：紅
    elseif Rum == 150
        cIN = [0.05 0.10 0.95];                       % single(fix) R150：深藍
    elseif USE_BIAS
        cIN = [0.000 0.549 0.600];
    else
        cIN = [0.482 0.322 0.671];
    end
    cM1  = [0.00 0.00 0.00];    cM2  = [0.00 0.60 0.00];   % mean 虛線一律中性色：黑 / 深綠
    ALPH = 0.60;
    % [MODIFIED 2026-08-20 使用者拍板] 原本 `nb = 180` 固定 **bin 數量**、跨距由該圖自己的
    %   max 決定 -> R150 與 R300 兩張的 bin 寬差 17.2 倍（0.613 vs 10.53 uT），並排比較時
    %   面積 = 100% x bin 寬不相等、形狀不可比（figure-style「同一張圖 bin 寬必須相同」
    %   的跨圖版）。改成固定 **bin 寬**、與 plot_err_hist_overlay 同值，整個 err_hist
    %   家族共用一把尺。
    % [MODIFIED 2026-08-20 使用者拍板] bin 寬對齊 `plot_err_hist_overlay.m` 的 2.8 uT：
    %   這張 R150 圖跟 overlay 的**紅色**是同一組資料（eighteen、R<=150、mean 0.0277 mT），
    %   解析度與顏色都要一致，否則同一組資料在兩張圖長不一樣。
    BINW = 2.8e-3;                                    % 共用 bin 寬 [mT]（= 2.8 uT，同 overlay）

    eAll = [e1; e2];                                  % 全部殘差（分層與否都由它定 edges / x 上界）
    maxE = max(eAll);
    edg  = 0 : BINW : (ceil(maxE/BINW)*BINW);         % 共用 edges（跨圖同寬）
    ctr  = (edg(1:end-1) + edg(2:end)) / 2;
    fprintf(['bin 寬 %.2f uT：全域 %d 個 bin（資料 max %.4f mT）' char(10)], BINW*1e3, numel(edg)-1, maxE);

    % [ADDED 2026-08-20 使用者拍板] 黑邊不要「一條一條畫」：座標軸寬約 660 pt，
    %   視野內 N 根長條時每根只有 660/N pt（N=714 -> 0.92 pt），而邊線 0.3 pt
    %   —— 黑線會佔掛一大截、整張糊成一團黑。
    %   ⚙ 提高匯出 DPI 沒用——長條寬與線寬都是「點」，比例不變。
    %   解法（沿用 plot_err_hist_overlay 已定案的做法）：長條**填色不畫邊**，改成
    %   ① stairs 畫一條連續的**頂部輪廓**（這才是直方圖的形狀）、
    %   ② 每 stepV 根才畫一條**全高垂直黑線**。垂直線疑度降下來，但輪廓還在。
    [xr, xt, fout, p995] = xlim_auto(eAll);           % 先算 x 視窗（才知道長條可見寬度）
    AXW_PT = 0.80 * 1100 * 72/96;                     % 座標軸寬約值 [pt]
    nvis   = max(1, ceil(diff(xr)/BINW));             % 視野內長條數
    barw   = AXW_PT / nvis;                           % 每根寬度 [pt]
    % [MODIFIED 2026-08-20 使用者拍板] **黑色線段全部拿掉**：頂部輪廓、垂直線、
    %   堆疊分界線一律不畫，只留純填色長條。（要復原就把 DRAW_OUTLINE 改回 true）
    DRAW_OUTLINE = false;
    stepV  = max(1, round(8.0/barw));                 % 垂直黑線間距目標 ~8 pt（對齊 overlay 那張：
                                                      %   1.95 pt 長條 x stepV=4 = 7.8 pt）；
                                                      %   長條夠寬時 stepV=1，等同每根都畫邊。
    fprintf(['  視野內 %d 根、每根 %.2f pt -> 每 %d 根一條垂直黑線（間距 %.1f pt）' char(10)], ...
            nvis, barw, stepV, barw*stepV);

    FS   = 28;
    fig  = figure('Color','w','Position',[100 100 1100 830]);
    ax   = axes(fig);   hold(ax,'on');

    if strcmp(MODE,'overlay')
        p1 = histcounts(e1, edg) / numel(e1) * 100;   % 各組自我正規化（樣本數差 7 倍）
        p2 = histcounts(e2, edg) / numel(e2) * 100;
        h1 = bar(ax, ctr, p1, 1, 'FaceColor',cIN,  'FaceAlpha',ALPH, 'EdgeColor','none');
        h2 = bar(ax, ctr, p2, 1, 'FaceColor',cOUT, 'FaceAlpha',ALPH, 'EdgeColor','none');
        % [MODIFIED 2026-08-17 使用者拍板] 不再畫 mean 虛線（見檔頭說明）
        if DRAW_OUTLINE, hist_outline(ax, edg, p1, stepV);  hist_outline(ax, edg, p2, stepV); end
        pAll = [p1 p2];
    elseif strcmp(MODE,'stack')
        % **共用分母 numel(eAll)** → 兩段疊起來剛好等於「不分層」時的整體分布，
        %   只是每根長條依格點半徑分成內層（藍）／外層（紅）兩段。
        p1 = histcounts(e1, edg) / numel(eAll) * 100;
        p2 = histcounts(e2, edg) / numel(eAll) * 100;
        hb = bar(ax, ctr, [p1(:) p2(:)], 1, 'stacked', 'EdgeColor','none');
        hb(1).FaceColor = cIN;    hb(1).FaceAlpha = ALPH;
        hb(2).FaceColor = cOUT;   hb(2).FaceAlpha = ALPH;
        h1 = hb(1);   h2 = hb(2);
        muA = mean(eAll);                             % [MODIFIED 2026-08-17] 不再畫 mean 虛線
        pAll = p1 + p2;                               % 堆疊後的總高度
        if DRAW_OUTLINE
            hist_outline(ax, edg, pAll, stepV);       % 輪廓 + 垂直線畫在堆疊總高
            stairs(ax, edg, [p1 p1(end)], 'Color','k', 'LineWidth',0.5, 'HandleVisibility','off');
        end
        fprintf('  堆疊：整體 mean=%.4f mT；內層佔 %.1f%%、外層佔 %.1f%% 的樣本\n', ...
                muA, 100*numel(e1)/numel(eAll), 100*numel(e2)/numel(eAll));
    else
        pA = histcounts(eAll, edg) / numel(eAll) * 100;
        hA = bar(ax, ctr, pA, 1, 'FaceColor',cIN, 'FaceAlpha',ALPH, 'EdgeColor','none');
        muA = mean(eAll);                             % [MODIFIED 2026-08-17] 不再畫 mean 虛線
        pAll = pA;
        if DRAW_OUTLINE, hist_outline(ax, edg, pA, stepV); end
        fprintf('  不分層：mean=%.4f mT  max=%.4f mT（%d 個殘差樣本）\n', muA, maxE, numel(eAll));
    end

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.015 .015],'TickDir','out');
    xlim(ax, xr);   set(ax,'XTick', xt);      % xr/xt 已在上面算好
    if fout > 0
        fprintf('  x 視窗 [0, %g] mT（99.5 百分位 %.4f mT 取整而來）；超出視野的長尾佔 %.2f%%、最大 %.3f mT\n', ...
                xr(2), p995, fout, maxE);
    end
    % [MODIFIED 2026-08-20] 原本 `round(linspace(0,ytop,5),1)` 的四舍五入會**破壞等距**
    %   （ytop=1.3 -> 0.325/0.65/0.975 -> 0.3/0.7/1.0，間隔 0.4 跟 0.3）——正是 figure-style
    %   明文警告的坑。改用跟 plot_sigma_hist / plot_gain_iso_hist 同一支的
    %   ylim_from_zero：自 0 起、上線只留 8% 裕度、nice 等距步長、刻度取奇數個。
    % [MODIFIED 2026-08-20] 內部刻度 3 -> **4**：[0,T] 平分 N+1 段的規則下，N 越大
    %   T 越貼近資料峰值 -> 上方留白變少（R300 由 75% 提到 80%）。
    [yr, yt] = ylim_from_zero(max(pAll), 4);
    % [MODIFIED 2026-08-20 使用者拍板] 終點「只標數字、**不畫 tick mark**」
    %   （與水平軸端點同一慣例）—— YTick 只放內部值，頂端那個用 text 補在軸外。
    % 內部刻度標數字（MATLAB 預設）；終點與起點都不標、也沒 tick。
    ylim(ax, yr);   set(ax,'YTick', yt);   ytop = yr(2);

    % x 起訖：只標數字、不畫 tick mark
    for xv = xr
        text(ax, xv, -0.022*ytop, sprintf('%g', xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    xlabel(ax, '$\mathbf{Residual\;(mT)}$',    'Interpreter','latex', 'FontSize',36);
    ylabel(ax, '$\mathbf{Percentage\;(\%)}$',  'Interpreter','latex', 'FontSize',36);

    % 圖例（使用者拍板 2026-08-19）：放**座標框內右上角**、每列一個系列，
    %   mean 併進該列（不畫 mean 虛線）。原本掛在框外上方並壓低座標軸的做法已移除。
    if strcmp(MODE,'overlay') || strcmp(MODE,'stack')
        % [MODIFIED 2026-08-20 使用者拍板] 圖例拿掉 Mean（數值仍印在 console 供追溯）
        lb = { ['R \leq ' num2str(Rsplit) ' {\mu}m'], ...
               [num2str(Rsplit) ' < R \leq ' num2str(Rum) ' {\mu}m'] };
        hh = [h1 h2];
        % [ADDED 2026-08-21 使用者拍板] **空的層不進圖例**：Rsplit == Rum 時（例 R150 切 150）
        %   外層 0 筆，舊版照樣列一個空的藍色色塊、標籤還寫成荒謬的「150 < R <= 150 um」。
        %   改成依樣本數篩掉空組。
        drop = [isempty(e1), isempty(e2)];
        hh(drop) = [];   lb(drop) = [];
        lg = legend(ax, hh, lb, 'Interpreter','tex', ...
                    'Location','northeast', 'NumColumns',1);
    else
        lg = legend(ax, hA, {['R \leq ' num2str(Rum) ' {\mu}m']}, ...
                    'Interpreter','tex', 'Location','northeast', 'NumColumns',1);
    end
    lg.FontSize = 24;  lg.FontWeight = 'bold';
    lg.Box = 'on';     lg.EdgeColor = 'k';   lg.LineWidth = 2.5;   lg.Color = 'w';
    ax.Toolbar.Visible = 'off';   hold(ax,'off');

    mstr = 'single';  if USE_BIAS, mstr = 'eighteen'; end
    % shell = 各自正規化的疊圖；conv = 單一分布（不分層 或 依半徑分色堆疊）
    tag  = 'conv';    if strcmp(MODE,'overlay'), tag = 'shell'; end
    out  = fullfile(figdir, sprintf('err_hist_%s_maxwell_%s_R%d.png', tag, mstr, Rum));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function [lim, tk] = ylim_from_zero(maxv, N)
% 縱軸（使用者拍板 2026-08-20，定案）：
%   • 自 0 起；**終點（軸上緣）不標數字、也不畫 tick**；起點 0 也不標。
%   • **內部刻度要標數字**，且每個數字必須是**整數或一位小數 0.N**。
%   • N = 內部刻度數（**不含起點與終點**）。
%   • 內部刻度彼此等距，**且與起點 0、終點 T 也等距** -> [0,T] 平分 N+1 段，
%     所以 tk = (1:N)*s、T = (N+1)*s。
%
%   選 s 的方式：s 必須是 **0.1 的倍數**（這樣 (1:N)*s 才保證是整數或 0.N），
%   且 (N+1)*s >= 1.08*maxv。在 [smin, 1.15*smin] 範圍內挑**最漂亮**的：
%   整數 > 0.5 的倍數 > 0.2 的倍數 > 其餘。（只取最小的 s 填充率最高，
%   但會出現 1.1/2.2/3.3 這種難讀的刻度。窗口 1.15 倍 = 只用很小的填充率
%   換一組漂亮數字；**實測 1.35 倍會讓 R150 揉到 68%、上方留白太多**。）
%
%   例：maxv=4.066 -> smin=1.098 -> s=1.2 -> 刻度 1.2/2.4/3.6、T=4.8（不標）、填充 85%。
%       maxv=0.6305 -> smin=0.170 -> s=0.2 -> 刻度 0.2/0.4/0.6、T=0.8（不標）。
    if nargin < 2 || isempty(N), N = 3; end
    smin = 1.08*maxv/(N+1);
    k0   = max(1, ceil(smin/0.1 - 1e-9));            % 以 0.1 為單位的最小步數
    k1   = max(k0, ceil(1.15*smin/0.1));
    best = k0;   bs = -1;
    for k = k0:k1
        if     mod(k,10) == 0, sc = 3;               % 整數
        elseif mod(k,5)  == 0, sc = 2;               % 0.5 的倍數
        elseif mod(k,2)  == 0, sc = 1;               % 0.2 的倍數
        else,                  sc = 0;
        end
        if sc > bs, bs = sc;  best = k; end          % 同分取最小 k（填充率高）
    end
    s   = best*0.1;
    tk  = round((1:N)*s*10)/10;
    lim = [0 (N+1)*s];
end






% ============================================================================
function hist_outline(ax, edg, p, stepV)
% [ADDED 2026-08-20] 直方圖黑線：頂部連續輪廓 + 每 stepV 根一條全高垂直線。
%   做法完全沿用 plot_err_hist_overlay.m（使用者指定的參考圖）。
%   重點是「輪廓連續、垂直線稀疏」——長條很細時垂直線一條一條畫會糊成一團黑，
%   但頂部輪廓是水平的、不管多細都不會互相重疊。
    stairs(ax, edg, [p p(end)], 'Color','k', 'LineWidth',0.5, 'HandleVisibility','off');
    for i = 1:stepV:numel(p)
        if p(i) > 0
            line(ax, [edg(i) edg(i)], [0 p(i)], 'Color','k', 'LineWidth',0.5, ...
                 'HandleVisibility','off');
        end
    end
end

% ============================================================================
function [xr, xt, frac, maxE] = xlim_auto(eAll, PCT)
% 橫軸（殘差 mT）：0 起、等距內縮刻度（端點另以 text 標數字，照直方圖既有做法）。
% [MODIFIED 2026-08-14] 上界改以 **PCT 百分位**（預設 99.5）為準，不再用 max
%   （使用者反映水平軸拉太寬）：長尾極稀疏 —— R=300 的 max 是 1.90 mT，但那是極少數點，
%   用它定上界會讓分布主體擠在左側 1/3。被截到視窗外的比例仍由 frac 回報、console 印出。
    if nargin < 2 || isempty(PCT), PCT = 99.5; end
    maxE = prctile(eAll, PCT);
    cand = [1 2 2.5 5 10];
    k0 = floor(log10(max(maxE, realmin)/5));   s = [];
    for k = k0:(k0+3)
        for c = cand
            if 5*c*10^k >= maxE, s = c*10^k;  break; end
        end
        if ~isempty(s), break; end
    end
    if isempty(s), s = maxE/5; end
    xr = [0, 5*s];   xt = (1:4)*s;
    if xr(2) > 1.25*maxE                       % 右端留白過大 → 收緊
        s2 = nice_step(maxE/3.5);
        n  = max(3, floor(maxE/s2));
        up = (n + 0.5)*s2;
        while up < maxE, n = n + 1;  up = (n + 0.5)*s2; end
        xr = [0, up];   xt = (1:n)*s2;
    end
    frac = 100 * mean(eAll > xr(2));
end

function s = nice_step(x)
    k = floor(log10(x));   m = x/10^k;
    cand = [1 2 2.5 5 10];
    [~, i] = min(abs(cand - m));
    s = cand(i)*10^k;
end

% ============================================================================
function Pc = make_Pc(e17, Pc_base)
% 電荷格 Pc = Pc_base + E(e)（含 e6z 約束；與 solve_current 內的版本一致）
    if isempty(e17) || all(e17(:) == 0), Pc = Pc_base;  return; end
    E = zeros(3, 6);
    E(:,1) = e17(1:3);     E(:,2) = e17(4:6);
    E(:,3) = e17(7:9);     E(:,4) = e17(10:12);
    E(:,5) = e17(13:15);
    E(1,6) = e17(16);      E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);
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
