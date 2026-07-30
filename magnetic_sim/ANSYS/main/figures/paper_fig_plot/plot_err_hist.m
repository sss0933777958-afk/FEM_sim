function plot_err_hist(USE_BIAS)
% plot_err_hist -- long2016 半切六極 R=150µm 校正的「絕對殘差直方圖」
% =========================================================================
%   graded / current、R=150µm 校正球。逐節點(每 node×每激發)絕對殘差 (mT):
%     |b_FEM − Sᵢ·ᴮĝ_I·K̄·I_j| = ||S·G − Bstack||_node   (S 用擬合後電荷位置 Pc)。
%   直方圖:亮藍填色 + 細黑邊、mean 虛線、圖例(sampling range + mean)、percentage 縱軸、無軸標題。
%   USE_BIAS=false → fix(單一 ℓ)→ err_hist.png；true → 18-param bias → err_hist_bias.png。
%   同時印 RMSPE = sqrt(Σε²/Σb²)·100(與 PDF 一致)。
% =========================================================================
    clc;
    if nargin < 1, USE_BIAS = false; end   % false = fix → err_hist.png；true = bias → err_hist_bias.png
    here   = fileparts(mfilename('fullpath'));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CAL,'function'));  addpath(fullfile(CAL,'common_path'));

    % ---- 前段 + R=150 fix 擬合 ----
    cfg = model_config('long2016_hexapole_halfcut','tip40um');
    raw = extract_ansys_data(cfg, 'all', 'graded');
    ad  = build_actuator_data(raw, cfg);   Pc_base = ad.Pc_base;
    F   = zeros(6, cfg.N_I);  for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end
    Rum = 150;
    [P, Bstack, npts] = cfg.select_ball(ad, Rum*1e-6);
    [e, l_hat, ~] = fitting(P, Bstack, Pc_base, 0.5e-3, USE_BIAS);
    [~, ~, G] = solve_current(l_hat, e, Pc_base, P, Bstack, F);

    % ---- 重建 S(用擬合後電荷位置 Pc = Pc_base + E(e)；fix 時 e=0 → Pc=Pc_base)、模型場、殘差 ----
    Pc = make_Pc(e, Pc_base);
    pbar = P / l_hat;  S = zeros(3*npts, 6);
    for k = 1:6, d = pbar - Pc(:,k).'; r3 = sum(d.^2,2).^1.5; S(:,k) = reshape((d./r3).', 3*npts, 1); end
    resid = S*G - Bstack;                         % 3npts×6

    % ---- 逐節點×激發 絕對殘差 |b_FEM − S·ĝ_I·K̄·I| (mT) ----
    %   resid = S*G − Bstack；S·G ≡ S·(ĝ_I·K̄_I·I)（G = ĝ_I·K̄_I·F 的重排），故此殘差即 |b_FEM − 模型|
    err = zeros(npts, 6);
    for j = 1:6
        rj = reshape(resid(:,j), 3, npts);   % 殘差向量 (mT)
        err(:,j) = sqrt(sum(rj.^2,1)).';     % 每節點×激發 殘差大小 (mT)
    end
    err = err(:);
    RMSPE = sqrt(sum(resid(:).^2) / sum(Bstack(:).^2)) * 100;
    mu = mean(err);
    fprintf('npts=%d  RMSPE=%.3f%%  mean|resid|=%.4f mT  max=%.4f mT\n', npts, RMSPE, mu, max(err));

    % ---- 直方圖(percentage 縱軸;亮藍 + 粗黑邊)----
    nb  = 180;  maxE = max(err);                       % 全範圍、密(依直方圖規則密度)
    edg = linspace(0, maxE, nb+1);
    [cnt, edg] = histcounts(err, edg);
    pct = cnt / numel(err) * 100;
    ctr = (edg(1:end-1) + edg(2:end)) / 2;

    FS = 28;
    fig = figure('Color','w','Position',[100 100 1100 720]);  ax = axes(fig);  hold(ax,'on');
    hb = bar(ax, ctr, pct, 1, 'FaceColor',[0.10 0.35 1.00], 'FaceAlpha',0.95, ...
             'EdgeColor','k', 'LineWidth',0.3);                       % 亮藍 + 細黑邊(密)
    ml = xline(ax, mu, '--', 'Color',[0.85 0.10 0.10], 'LineWidth',3.0);  % mean 虛線

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.015 .015],'TickDir','out');
    if ~USE_BIAS
        xlim(ax,[0 1.0]);   set(ax,'XTick',[0.2 0.4 0.6 0.8]);       % fix:橫軸 mT，稍超 range(0.88)、tick 0.2~0.8（起始點不畫 tick）
    else
        xlim(ax,[0 0.25]);  set(ax,'XTick',[0.05 0.10 0.15 0.20]);   % bias:橫軸 mT，稍超 range(0.21)、tick 0.05~0.20（起始點不畫 tick）
    end
    ymax = max(pct);  ytop = ceil(ymax*10)/10;  ylim(ax,[0 ytop]);
    yd = round(linspace(0, ytop, 5), 1);  set(ax,'YTick', yd(2:4));  % 3 內縮 tick(首末端點不標)
    % 起始點 0 + 終點：只標數字、不畫 tick（左右角落補字，text）
    xr = xlim(ax);
    text(ax, xr(1), -0.022*ytop, sprintf('%g',xr(1)), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS,'FontWeight','bold','Clipping','off');   % 起點
    text(ax, xr(2), -0.022*ytop, sprintf('%g',xr(2)), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS,'FontWeight','bold','Clipping','off');   % 終點
    lg = legend([hb ml], {sprintf('Sampling range \\leq %d {\\mu}m', Rum), sprintf('Mean = %.3f mT', mu)}, ...
                'Interpreter','tex', 'Location','northeast');
    lg.FontSize = 24;  lg.FontWeight = 'bold';
    xlabel(ax, '$\mathbf{\| b_{FEM} - S_i\, {}^{B}\hat{g}_{I}\, \bar{K}\, I_{j} \|\;(mT)}$', ...
           'Interpreter','latex', 'FontSize',36);      % 絕對殘差軸標題(標準數學字體、36 粗)
    ax.Toolbar.Visible = 'off';  hold(ax,'off');       % 刻度數字保留

    bstr = ''; if USE_BIAS, bstr = '_bias'; end
    out = fullfile(figdir, sprintf('err_hist%s.png', bstr));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
end

% ---- 電荷格 Pc = Pc_base + E(e)（含 e6z 約束；與 solve_current 內 make_Pc 一致）----
function Pc = make_Pc(e17, Pc_base)
    if isempty(e17) || all(e17(:)==0)   % fix：無偏移
        Pc = Pc_base;  return;
    end
    E = zeros(3, 6);
    E(:,1) = e17(1:3);     E(:,2) = e17(4:6);
    E(:,3) = e17(7:9);     E(:,4) = e17(10:12);
    E(:,5) = e17(13:15);
    E(1,6) = e17(16);      E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);
    Pc = Pc_base + E;
end
