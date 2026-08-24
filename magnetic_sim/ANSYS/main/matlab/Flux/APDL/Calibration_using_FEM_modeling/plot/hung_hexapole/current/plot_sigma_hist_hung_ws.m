function plot_sigma_hist_hung_ws()
% PLOT_SIGMA_HIST_HUNG_WS  讀 utils 產的空間 SVD σ .mat → σ1/σ2/σ3 疊圖直方圖（R300/R700 各一張）。
%   分工（見 calibration-shared-structure）：本腳本只「讀 .mat → 畫圖」，SVD 計算在
%   utils/hung_hexapole/svd_sigma_hung_ws.m（讀校正結果算 σ、存 data/.mat）。
%   資料：data/hung_hexapole/.mat/svd_sigma_hung_ws_R150.mat（每變體 sigma1/2/3，R≤150µm 球內均勻取樣）。
%   σ = svd(T)，T = S(p)·Ĥ_I（3×6），單位 mT/A；σ1≥σ2≥σ3。
%   風格：選項① 粗體框、共用 bins(nb=180)、percentage 縱軸、mean 虛線、legend 含 mean+CV。
%   輸出：figures/hung_hexapole/current/eighteen/sigma_hist_{R300,R700}.png（覆蓋迭代）。

    here = fileparts(mfilename('fullpath'));                          % plot/hung_hexapole/current
    CAL  = fileparts(fileparts(fileparts(here)));                     % Calibration_using_FEM_modeling
    matf = fullfile(CAL,'data','hung_hexapole','.mat','svd_sigma_hung_ws_R150.mat');
    assert(exist(matf,'file')==2, '缺 %s（先跑 utils/hung_hexapole/svd_sigma_hung_ws.m）', matf);
    L = load(matf);

    figdir = fullfile(CAL,'figures','hung_hexapole','current','eighteen');
    if ~exist(figdir,'dir'), mkdir(figdir); end

    cols = [0.80 0.10 0.10;     % σ1 紅
            0.10 0.55 0.20;     % σ2 綠
            0.12 0.30 0.75];    % σ3 藍
    labs = {'\sigma_1','\sigma_2','\sigma_3'};

    for c = 1:numel(L.variants)
        V  = L.variants{c};   sd = L.(V);
        Y  = {sd.sigma1, sd.sigma2, sd.sigma3};
        fprintf('%s: N=%d, mean σ=[%.3f %.3f %.3f] mT/A\n', V, sd.npts, mean(sd.sigma1), mean(sd.sigma2), mean(sd.sigma3));
        render_overlayN(Y, '$\sigma\;[\mathrm{mT/A}]$', cols, labs, ...
                        fullfile(figdir, sprintf('sigma_hist_%s.png', V)), false);
    end
end

% ============================================================================
function render_overlayN(Y, xlab, cols, labels, outpng, xlog)
% N 條疊圖（共用 bins、percentage 縱軸、mean 虛線、legend 含 mean + CV、選項① 粗體框）。
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
        leg{k} = sprintf('%s  (mean = %.4g)', labels{k}, mean(Y{k}));   % 只留 mean（tex；label 如 \sigma_1）
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
    lg = legend(ax, h, leg, 'Location','northeast');   % tex（預設）：\sigma_1 → σ₁、% 不被吃
    set(lg,'FontSize',15,'FontWeight','bold','Box','on');
    ax.Toolbar.Visible = 'off';
    exportgraphics(fig, outpng, 'Resolution',600);
    fprintf('saved %s\n', outpng);
    close(fig);
end
