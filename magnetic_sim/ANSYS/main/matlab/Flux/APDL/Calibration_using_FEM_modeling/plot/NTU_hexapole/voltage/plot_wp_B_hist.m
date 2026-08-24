function plot_wp_B_hist(R_EVAL_um)
% plot_wp_B_hist — WP 球內（R <= R_EVAL_um）逐節點 |B| 的分布：襯套導磁 vs 不導磁 疊圖。
%   資料：coil1 自激發（P1），**真實 FEM 節點原值、不內插**（per plot-real-nodes）。
%   WP = (0,0,9.30) mm（6 tip 形心）；R<=150µm 球內全在空氣（下極板頂 9.05 / 上極板底 9.55）。
%   橫軸 |B| [mT]、縱軸 **百分比 = 該 bin 節點數 / 該 variant 球內總節點數 × 100**（使用者定義）。
%   風格：離散 bars、**兩組共用 edges**、nb=180、FaceAlpha 0.55、mean 用 xline 虛線 + legend 標值。
%
%   出 1 檔 → figures/shared/wp_B_hist.png

    if nargin < 1 || isempty(R_EVAL_um), R_EVAL_um = 150; end

    % ======================== 參數 ========================
    NB   = 180;                 % bin 數（per figure-style「直方圖 bin 間距」）
    WP   = [0, 0, 9.30];        % 工作點 [mm]（= mt_constants 的 c.WP）
    % =====================================================

    here   = fileparts(mfilename('fullpath'));
    vbase  = fileparts(fileparts(here));
    figdir = fullfile(vbase,'figures','shared');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    % [MODIFIED 2026-08-08] 脫離 backup（規則 no-backup-data）→ live 樹。
    CALROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CALROOT,'function'), fullfile(CALROOT,'utils'), fullfile(CALROOT,'common_path'));
    droot = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\NTU_hexapole\data';

    variants = {'full_assembly_sleeve','full_assembly'};
    vlab     = {'sleeve \mu_r = 280', 'sleeve \mu_r = 1 (non-magnetic)'};
    cols     = [0.85 0.10 0.10; 0.20 0.40 0.90];

    % ---- Pass 1：兩 variant 取球內節點 |B| ----
    S = cell(1,2);  allB = [];
    for v = 1:2
        d = import_ansys_data(fullfile(droot,variants{v},'coil1'),'all','coil1');
        X = [d.x d.y d.z]*1e3;  B3 = [d.bx d.by d.bz]*1e3;              % mm / mT
        in = vecnorm(X - WP, 2, 2) <= R_EVAL_um*1e-3;
        S{v}.B = vecnorm(B3(in,:), 2, 2);
        allB = [allB; S{v}.B];                                          %#ok<AGROW>
        fprintf('%-22s: %5d nodes in R<=%gum | |B| mean %.4f  median %.4f  [%.4f %.4f] mT\n', ...
                variants{v}, numel(S{v}.B), R_EVAL_um, mean(S{v}.B), median(S{v}.B), ...
                min(S{v}.B), max(S{v}.B));
    end
    edges = linspace(min(allB), max(allB), NB+1);                       % 共用 edges
    fprintf('共用 edges = [%.4f %.4f] mT | nb = %d | bin width = %.5f mT\n', ...
            edges(1), edges(end), NB, edges(2)-edges(1));
    fprintf('mean 比：sleeve / no-sleeve = %.4f (%+.2f%%)\n', ...
            mean(S{1}.B)/mean(S{2}.B), 100*(mean(S{1}.B)/mean(S{2}.B)-1));

    % ---- Pass 2：疊圖（縱軸 = 數量/該 variant 總數 × 100）----
    % histogram 無「百分比」Normalization → 自己 histcounts 算 pct，用 BinEdges/BinCounts 畫
    f = figure('Color','w','Position',[100 100 1180 850]); hold on;
    h = gobjects(1,2);  ymax = 0;
    for v = 1:2
        cnt = histcounts(S{v}.B, edges);
        pct = 100 * cnt / numel(S{v}.B);
        h(v) = histogram('BinEdges',edges, 'BinCounts',pct, ...
                         'FaceColor',cols(v,:), 'FaceAlpha',0.55, ...
                         'EdgeColor','w', 'LineWidth',0.3);
        ymax = max(ymax, max(pct));
    end
    for v = 1:2
        xline(mean(S{v}.B), '--', 'Color',cols(v,:), 'LineWidth',2.2);
    end
    legend(h, arrayfun(@(v) sprintf('%s   (mean = %.3f mT)', vlab{v}, mean(S{v}.B)), 1:2, 'uni',0), ...
           'Location','northeast','FontSize',14,'FontWeight','bold','Box','on');

    % ---- 風格：粗體框選項① ----
    ax = gca;
    ylim([0 ymax*1.20]);                                                % headroom：legend 不壓 bar
    grid(ax,'off'); box(ax,'on');
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    xlabel(ax,'|B| (mT)','FontWeight','bold');
    ylabel(ax,'percentage (%)','FontWeight','bold');
    xt=get(ax,'XTick'); if numel(xt)>3; set(ax,'XTick',xt(1:2:end)); end
    yt=get(ax,'YTick'); if numel(yt)>3; set(ax,'YTick',yt(1:2:end)); end
    ax.Toolbar.Visible='off';

    out = fullfile(figdir,'wp_B_hist.png');
    exportgraphics(f, out, 'Resolution', 150);  close(f);
    fprintf('saved %s\n', out);
end
