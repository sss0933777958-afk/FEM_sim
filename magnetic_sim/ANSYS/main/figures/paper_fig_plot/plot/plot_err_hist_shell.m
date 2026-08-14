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
%   直方圖規則（figure-style「分布 / 疊圖直方圖」）：nb=180、兩組**共用 edges**、
%   FaceAlpha、mean 以 xline 虛線標、圖例兩欄（左欄系列、右欄該系列的統計值）。
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
    if nargin < 4 || isempty(SPLIT),    SPLIT    = true;  end   % false = 不分層、單一藍色

    here   = fileparts(fileparts(mfilename('fullpath')));       % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    addpath(fullfile(CAL,'function'), fullfile(CAL,'common_path'), fullfile(CAL,'utils'));
    addpath(fullfile(CAL,'utils','long2016_hexapole_halfcut'));

    %% ---- ① 取該 R 的收斂設計 → 減量校正 -----------------------------------
    cfg = model_config('long2016_hexapole_halfcut','tip40um');
    F   = zeros(6, cfg.N_I);
    for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end

    tri = conv_design(Rum*1e-6, cfg, F, USE_BIAS, here);   % 共用引擎（先查快取、沒有才搜）

    [gx,gy,gz,gB] = sphere_grid_sample(Rum*1e-6, [], ...
                        struct('frame','actuator','NRPT',tri));
    Pc_ = [gx gy gz];   Nc = size(Pc_,1);
    Bc  = zeros(3*Nc, size(gB,3));
    for j = 1:size(gB,3), Bc(:,j) = reshape(gB(:,:,j).', [], 1); end
    [e, l_hat] = fitting(Pc_, Bc, cfg.Pc_base, 0.5e-3, USE_BIAS);
    [~, gI, G] = solve_current(l_hat, e, cfg.Pc_base, Pc_, Bc, F);
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
    cOUT = [0.85 0.10 0.10];
    if SPLIT || Rum == 150
        cIN = [0.05 0.10 0.95];
    elseif USE_BIAS
        cIN = [0.000 0.549 0.600];
    else
        cIN = [0.482 0.322 0.671];
    end
    cM1  = [0.00 0.00 0.00];    cM2  = [0.00 0.60 0.00];   % mean 虛線一律中性色：黑 / 深綠
    ALPH = 0.60;                nb   = 180;                % 依直方圖規則

    eAll = [e1; e2];                                  % 全部殘差（分層與否都由它定 edges / x 上界）
    maxE = max(eAll);
    edg  = linspace(0, maxE, nb+1);                   % （分層時）兩組共用 edges
    ctr  = (edg(1:end-1) + edg(2:end)) / 2;

    FS   = 28;
    fig  = figure('Color','w','Position',[100 100 1100 830]);
    ax   = axes(fig);   hold(ax,'on');

    if SPLIT
        p1 = histcounts(e1, edg) / numel(e1) * 100;   % 各組自我正規化（樣本數差 7 倍）
        p2 = histcounts(e2, edg) / numel(e2) * 100;
        h1 = bar(ax, ctr, p1, 1, 'FaceColor',cIN,  'FaceAlpha',ALPH, 'EdgeColor','k', 'LineWidth',0.3);
        h2 = bar(ax, ctr, p2, 1, 'FaceColor',cOUT, 'FaceAlpha',ALPH, 'EdgeColor','k', 'LineWidth',0.3);
        m1 = xline(ax, mu1, '--', 'Color',cM1, 'LineWidth',2.8);
        m2 = xline(ax, mu2, '--', 'Color',cM2, 'LineWidth',2.8);
        pAll = [p1 p2];
    else
        pA = histcounts(eAll, edg) / numel(eAll) * 100;
        hA = bar(ax, ctr, pA, 1, 'FaceColor',cIN, 'FaceAlpha',ALPH, 'EdgeColor','k', 'LineWidth',0.3);
        muA = mean(eAll);
        mA  = xline(ax, muA, '--', 'Color',cM1, 'LineWidth',2.8);
        pAll = pA;
        fprintf('  不分層：mean=%.4f mT  max=%.4f mT（%d 個殘差樣本）\n', muA, maxE, numel(eAll));
    end

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.015 .015],'TickDir','out');
    [xr, xt, fout] = xlim_auto(eAll);
    xlim(ax, xr);   set(ax,'XTick', xt);
    if fout > 0
        fprintf('  x 視窗 [0, %g] mT（99.5 百分位取整）；超出視野的長尾佔 %.2f%%、最大 %.3f mT\n', ...
                xr(2), fout, maxE);
    end
    ytop = ceil(max(pAll)*10)/10;   ylim(ax,[0 ytop]);
    yd = round(linspace(0, ytop, 5), 1);   set(ax,'YTick', yd(2:4));   % 3 個內縮 tick

    % x 起訖：只標數字、不畫 tick mark
    for xv = xr
        text(ax, xv, -0.022*ytop, sprintf('%g', xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    xlabel(ax, '$\mathbf{Residual\;(mT)}$',    'Interpreter','latex', 'FontSize',36);
    ylabel(ax, '$\mathbf{Percentage\;(\%)}$',  'Interpreter','latex', 'FontSize',36);

    % 圖例：左欄 = 系列、右欄 = 該系列的統計值（house 標準）
    if SPLIT
        % R 大寫；外層標成明確區間 150 < R <= 300（不用開放的 R > 150）
        lg = legend(ax, [h1 h2 m1 m2], ...
            {sprintf('R \\leq %d {\\mu}m', Rsplit), ...
             sprintf('%d < R \\leq %d {\\mu}m', Rsplit, Rum), ...
             sprintf('Inner mean = %.3f mT', mu1),  sprintf('Outer mean = %.3f mT', mu2)}, ...
            'Interpreter','tex', 'Location','northoutside', 'NumColumns',2);
    else
        lg = legend(ax, [hA mA], ...
            {sprintf('R \\leq %d {\\mu}m', Rum), ...
             sprintf('Mean = %.3f mT', mean(eAll))}, ...
            'Interpreter','tex', 'Location','northoutside', 'NumColumns',2);
    end
    lg.FontSize = 24;  lg.FontWeight = 'bold';
    lg.Box = 'on';     lg.EdgeColor = 'k';   lg.LineWidth = 2.5;
    ax.Toolbar.Visible = 'off';   hold(ax,'off');

    drawnow;
    axp = get(ax,'Position');   lgp = get(lg,'Position');   lgh = lgp(4);
    GAPN = 0.022;   newTop = 1 - lgh - GAPN - 0.006;
    axp(4) = newTop - axp(2);   set(ax,'Position',axp);
    set(lg, 'Position', [axp(1), newTop + GAPN, axp(3), lgh]);

    mstr = 'single';  if USE_BIAS, mstr = 'eighteen'; end
    tag  = 'conv';    if SPLIT,    tag  = 'shell';    end     % shell = 分層、conv = 不分層
    out  = fullfile(figdir, sprintf('err_hist_%s_maxwell_%s_R%d.png', tag, mstr, Rum));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function [xr, xt, frac] = xlim_auto(eAll)
% 橫軸（殘差 mT）：0 起、等距內縮刻度（端點另以 text 標數字，照直方圖既有做法）。
%   沿用 plot_err_hist / plot_err_hist_overlay 的 xlim_pick：先取「5 格 nice 步長」，
%   若上界比實際 max 寬超過 25% 再收緊一次（使用者反映 x 範圍太寬）。
%   ⚠ 不用百分位裁切：本分布的尾巴很厚（R=300 時 p90=1.24、p95=1.63 mT），
%   砍尾會丟掉有意義的資料。收緊後 **一格不丟**、只是把空白壓掉。
    maxE = max(eAll);
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
