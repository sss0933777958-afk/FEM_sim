function plot_err_hist_overlay()
% plot_err_hist_overlay -- long2016 半切六極 R=150µm：fix vs bias 絕對殘差「疊圖」
% =========================================================================
%   同一取樣球(R=150µm)、同一絕對殘差定義 |b_FEM − S·ĝ_I·K̄·I| (mT)：
%     藍 = 無 bias(fix 單一 ℓ)、紅 = 有 bias(18-param)。
%   共用 fix 的 x 刻度(xlim[0 1]、tick 0.2~0.8 + 起終點)；共用 edges(nb=180)。
%   percentage 縱軸、圖例標顏色對應、**不標平均線**、無軸標題。
%   輸出 → figures/paper_fig/Section2_E/err_hist_overlay.png(覆蓋迭代)。
% =========================================================================
    clc;
    here   = fileparts(mfilename('fullpath'));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CAL,'function'));  addpath(fullfile(CAL,'common_path'));

    % ---- 前段 pipeline(只做一次)----
    cfg = model_config('long2016_hexapole_halfcut','tip40um');
    raw = extract_ansys_data(cfg, 'all', 'graded');
    ad  = build_actuator_data(raw, cfg);   Pc_base = ad.Pc_base;
    F   = zeros(6, cfg.N_I);  for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end
    Rum = 150;
    [P, Bstack, npts] = cfg.select_ball(ad, Rum*1e-6);

    % ---- 兩模型各擬合一次 → 逐節點×激發 絕對殘差 (mT) ----
    err0 = fit_abs_resid(P, Bstack, Pc_base, F, npts, false);   % single_param (無 bias)
    err1 = fit_abs_resid(P, Bstack, Pc_base, F, npts, true);    % eighteen_param (有 bias)
    cv0 = std(err0)/mean(err0)*100;   cv1 = std(err1)/mean(err1)*100;   % CV = σ/μ [%]（不受 bin 寬影響）
    fprintf('single_param  : mean=%.4f mT  max=%.4f mT  CV=%.1f%%\n', mean(err0), max(err0), cv0);
    fprintf('eighteen_param: mean=%.4f mT  max=%.4f mT  CV=%.1f%%\n', mean(err1), max(err1), cv1);

    % ---- 各自用自己的 edges(不共用;各自然解析度、nb=180),只是畫在同一 x 軸 ----
    nb  = 180;
    edg0 = linspace(0, max(err0), nb+1);  ctr0 = (edg0(1:end-1)+edg0(2:end))/2;   % single 0~0.88
    edg1 = linspace(0, max(err1), nb+1);  ctr1 = (edg1(1:end-1)+edg1(2:end))/2;   % eighteen 0~0.21（自己的細 bin）
    pct0 = histcounts(err0, edg0) / numel(err0) * 100;
    pct1 = histcounts(err1, edg1) / numel(err1) * 100;

    cB = [0.10 0.35 1.00];   cR = [0.85 0.10 0.10];     % 藍=無bias / 紅=有bias
    FS = 28;
    fig = figure('Color','w','Position',[100 100 1100 830]);  ax = axes(fig);  hold(ax,'on');   % [MODIFIED] 加高補償外置兩列圖例
    h0 = bar(ax, ctr0, pct0, 1, 'FaceColor',cB, 'FaceAlpha',0.60, 'EdgeColor','k', 'LineWidth',0.3);   % 亮藍(single) + 黑細邊
    h1 = bar(ax, ctr1, pct1, 1, 'FaceColor',cR, 'FaceAlpha',0.60, 'EdgeColor','none');                 % 紅(eighteen) 無邊；帶透明看得到交錯
    % 紅色：沿輪廓畫連續黑線(silhouette)+ 挑幾根 bar 從上到下畫垂直黑邊
    stairs(ax, edg1, [pct1 pct1(end)], 'Color','k', 'LineWidth',0.5);   % 頂部輪廓(粗度同垂直邊 0.5)
    stepV = 4;
    for i = 1:stepV:numel(pct1)
        if pct1(i) > 0
            line(ax, [edg1(i) edg1(i)], [0 pct1(i)], 'Color','k', 'LineWidth',0.5);   % 全高垂直邊
        end
    end
    % 兩個分布的 mean 虛線(深藍=single、深紅=eighteen)
    mu0 = mean(err0);  mu1 = mean(err1);
    ml0 = xline(ax, mu0, '--', 'Color',[0.00 0.00 0.00], 'LineWidth',2.8);   % single mean 0.250（黑，高對比）
    ml1 = xline(ax, mu1, '--', 'Color',[0.00 0.60 0.00], 'LineWidth',2.8);   % eighteen mean 0.038（深綠，高對比）

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.015 .015],'TickDir','out');
    xlim(ax,[0 1.0]);   set(ax,'XTick',[0.2 0.4 0.6 0.8]);          % 用 fix 刻度
    ymax = max([pct0 pct1]);  ytop = ceil(ymax*10)/10;  ylim(ax,[0 ytop]);
    yd = round(linspace(0, ytop, 5), 1);  set(ax,'YTick', yd(2:4));  % 3 內縮 tick(首末端點不標)
    % 起點 0 + 終點：只標數字、不畫 tick(左右角落補字)
    xr = xlim(ax);
    text(ax, xr(1), -0.022*ytop, sprintf('%g',xr(1)), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS,'FontWeight','bold','Clipping','off');
    text(ax, xr(2), -0.022*ytop, sprintf('%g',xr(2)), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS,'FontWeight','bold','Clipping','off');
    hLegR = bar(ax, NaN, NaN, 'FaceColor',cR, 'FaceAlpha',0.60, 'EdgeColor','k', 'LineWidth',1.2);   % 只給圖例：紅框加黑邊(實際紅 bar 仍無邊)
    % [MODIFIED] 圖例移到圖框外「上方」、**兩列**(NumColumns=2 為欄優先填充 → 條目 1,2 = 左欄兩列、
    %   3,4 = 右欄兩列，故 [h0 hLegR ml0 ml1] 恰為「第一列 Single、第二列 Eighteen」)。
    %   欄距不用白字墊 —— 下面會強制圖例寬度 = 座標軸寬度，兩欄自動平均散開。
    lg = legend([h0 hLegR ml0 ml1], {'Single parameter', 'Eighteen parameters', ...
                sprintf('Single mean = %.3f mT', mu0), sprintf('Eighteen mean = %.3f mT', mu1)}, ...
                'Interpreter','tex', 'Location','northoutside', 'NumColumns',2);
    lg.FontSize = 24;  lg.FontWeight = 'bold';  lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;
    xlabel(ax, '$\mathbf{Residual\;(mT)}$', ...
           'Interpreter','latex', 'FontSize',36);      % 絕對殘差軸標題(標準數學字體、36 粗)
    ylabel(ax, '$\mathbf{Percentage\;(\%)}$', ...
           'Interpreter','latex', 'FontSize',36);      % [ADDED] 百分比縱軸標題(同 err_hist 風格)
    ax.Toolbar.Visible = 'off';  hold(ax,'off');

    % [ADDED] 圖例框:①左右邊界對齊圖框(寬度=ax 寬、x 起點同 ax，兩欄自動平均散開)
    %                ②壓低 ax 上緣，讓圖例底線與圖框上緣「留明確間隙」、不相黏
    drawnow;
    axp = get(ax,'Position');  lgh = get(lg,'Position');  lgh = lgh(4);
    GAPN   = 0.022;                          % 圖例底線 ↔ 圖框上緣的間隙(normalized)
    newTop = 1 - lgh - GAPN - 0.006;         % ax 新上緣(頂端留 0.006 給圖例外框)
    axp(4) = newTop - axp(2);                % 只壓上緣、下緣不動
    set(ax, 'Position', axp);
    set(lg, 'Position', [axp(1), newTop + GAPN, axp(3), lgh]);

    out = fullfile(figdir,'err_hist_overlay.png');
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
end

% ---- 擬合 + 逐節點×激發 絕對殘差 |S·G − Bstack| (mT) ----
function err = fit_abs_resid(P, Bstack, Pc_base, F, npts, USE_BIAS)
    [e, l_hat, ~] = fitting(P, Bstack, Pc_base, 0.5e-3, USE_BIAS);
    [~, ~, G]     = solve_current(l_hat, e, Pc_base, P, Bstack, F);
    Pc = make_Pc(e, Pc_base);
    pbar = P / l_hat;  S = zeros(3*npts, 6);
    for k = 1:6, d = pbar - Pc(:,k).'; r3 = sum(d.^2,2).^1.5; S(:,k) = reshape((d./r3).', 3*npts, 1); end
    resid = S*G - Bstack;
    err = zeros(npts, 6);
    for j = 1:6, rj = reshape(resid(:,j), 3, npts); err(:,j) = sqrt(sum(rj.^2,1)).'; end
    err = err(:);
end

% ---- 電荷格 Pc = Pc_base + E(e)(與 solve_current 內一致；fix e=0 → Pc_base)----
function Pc = make_Pc(e17, Pc_base)
    if isempty(e17) || all(e17(:)==0), Pc = Pc_base; return; end
    E = zeros(3, 6);
    E(:,1) = e17(1:3);     E(:,2) = e17(4:6);
    E(:,3) = e17(7:9);     E(:,4) = e17(10:12);
    E(:,5) = e17(13:15);
    E(1,6) = e17(16);      E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);
    Pc = Pc_base + E;
end
