function plot_geomean_hist_hung_ws()
% PLOT_GEOMEAN_HIST_HUNG_WS  讀 σ .mat → 每點奇異值幾何平均 (σ1σ2σ3)^(1/3) 直方圖（R300/R700 各一張）。
%   分工（見 calibration-shared-structure）：SVD 已由 utils/hung_hexapole/svd_sigma_hung_ws.m 算好存 .mat；
%   本腳本只「讀 .mat → 逐點開三次方根 → 畫圖」（elementwise，非重算 SVD）。
%   g = (σ1·σ2·σ3)^(1/3) = 𝒞^(1/3)：等體積等向橢球半徑 = 濃縮成單一「有效增益」（單位 mT/A）。
%   風格：選項① 粗體框、nb=180、percentage 縱軸、mean 虛線、legend 含 mean。
%   輸出：figures/hung_hexapole/current/eighteen/geomean_hist_{R300,R700}.png（覆蓋迭代）。

    here = fileparts(mfilename('fullpath'));                          % plot/hung_hexapole/current
    CAL  = fileparts(fileparts(fileparts(here)));                     % Calibration_using_FEM_modeling
    matf = fullfile(CAL,'data','hung_hexapole','.mat','svd_sigma_hung_ws_R150.mat');
    assert(exist(matf,'file')==2, '缺 %s（先跑 utils/hung_hexapole/svd_sigma_hung_ws.m）', matf);
    L = load(matf);

    figdir = fullfile(CAL,'figures','hung_hexapole','current','eighteen');
    if ~exist(figdir,'dir'), mkdir(figdir); end

    col = [0.40 0.18 0.55];    % 幾何平均（有效增益）單色

    for c = 1:numel(L.variants)
        V  = L.variants{c};   sd = L.(V);
        g  = (sd.sigma1 .* sd.sigma2 .* sd.sigma3).^(1/3);           % 逐點幾何平均 [mT/A]
        fprintf('%s: N=%d, mean (σ1σ2σ3)^{1/3} = %.4g mT/A\n', V, sd.npts, mean(g));
        render_overlayN({g}, '$(\sigma_1\sigma_2\sigma_3)^{1/3}\;[\mathrm{mT/A}]$', col, {V}, ...
                        fullfile(figdir, sprintf('geomean_hist_%s.png', V)), false);
    end
end

% ============================================================================
function render_overlayN(Y, xlab, cols, labels, outpng, xlog)
% N 條疊圖（共用 bins、percentage 縱軸、mean 虛線、legend 含 mean、選項① 粗體框）。單條亦可。
    if nargin < 6, xlog = false; end
    nb = 180;
    allY = cat(1, Y{:});
    if xlog
        edges = logspace(log10(min(allY)), log10(max(allY)), nb+1);
    else
        edges = linspace(min(allY), max(allY), nb+1);
    end
    fig = figure('Color','w','Position',[100 100 820 600]);  ax = axes(fig);  hold(ax,'on');
    h = gobjects(1,numel(Y));  leg = cell(1,numel(Y));
    for k = 1:numel(Y)
        h(k) = histogram(ax, Y{k}, edges, 'Normalization','percentage', ...
                         'FaceColor',cols(k,:), 'FaceAlpha',0.55, 'EdgeColor','w', 'LineWidth',0.3);
        leg{k} = sprintf('%s  (mean = %.4g, CV = %.2f%%)', labels{k}, mean(Y{k}), 100*std(Y{k})/mean(Y{k}));   % mean + CV（tex：% 不被吃）
    end
    for k = 1:numel(Y)
        xline(ax, mean(Y{k}), '--', 'Color',cols(k,:), 'LineWidth',2);
    end
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    box(ax,'on'); grid(ax,'off');
    if xlog, set(ax,'XScale','log');  xlim(ax,[edges(1) edges(end)]); end
    yl = ylim(ax); ylim(ax,[0 yl(2)*1.20]);
    yt = get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    xint = 'tex'; if startsWith(strtrim(xlab),'$'), xint = 'latex'; end
    xlabel(ax, xlab, 'FontWeight','bold','Interpreter',xint);
    ylabel(ax,'Percentage (%)','FontWeight','bold');
    lg = legend(ax, h, leg, 'Location','northeast');   % tex（預設）：% 不被吃
    set(lg,'FontSize',15,'FontWeight','bold','Box','on');
    ax.Toolbar.Visible = 'off';
    exportgraphics(fig, outpng, 'Resolution',600);
    fprintf('saved %s\n', outpng);
    close(fig);
end
