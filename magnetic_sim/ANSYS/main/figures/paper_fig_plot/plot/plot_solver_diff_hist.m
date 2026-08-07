function plot_solver_diff_hist(force, WHICH)
%PLOT_SOLVER_DIFF_HIST  APDL vs Maxwell 在**同一組位置**的場差直方圖（整個工作區、六個激發合併）。
%   plot_solver_diff_hist()            % 有快取就用，畫大小誤差
%   plot_solver_diff_hist(true)        % 強制重算（載兩套場 + 建內插器，約 3–5 分鐘）
%   plot_solver_diff_hist([], WHICH)   % WHICH = 'mag'（預設）| 'vec' | 'ang'
%
%   比較方式：以 **Maxwell 的 20 µm 規則格點**（R≤500 µm、已濾鐵）為查詢位置，
%   把 **APDL graded 場線性內插**到同一批點，六個激發逐點比較：
%     'mag'  大小誤差(相對) x = (‖b_M‖−‖b_A‖)/‖b_A‖×100  [%]（預設；有號）
%     'magmT'大小誤差(絕對) x = | ‖b_M‖ − ‖b_A‖ |         [mT]
%     'vec'  向量差         x = ‖ b_M − b_A ‖             [mT]
%     'ang'  方向誤差       x = ∠(b_A, b_M)               [deg]
%   縱軸一律 count / total_node × 100 [%]。
%
%   ⚠ **這是內插比較**（兩邊格點不重合）；內插會平滑掉高頻差異 → 差值是**下限**。
%   ⚠ coil→pole map：APDL `[1,3,6,5,2,4]`（欄為 coil 序）、Maxwell identity。比較前必須用
%     paper2coil 對齊，否則會拿 P3 去比 P2（實測會出現 ~90° 的假訊號）。
%
%   風格①粗體框圖 + 專案直方圖慣例（nb=180、mean 用 xline 虛線）。
%   輸出 figures/paper_fig/Section2_E/solver_diff_<which>_hist_R500.png

    if nargin < 1 || isempty(force), force = false; end
    if nargin < 2 || isempty(WHICH), WHICH = 'mag';  end
    here   = fileparts(fileparts(mfilename('fullpath')));          % .../paper_fig_plot
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'), mkdir(figdir); end
    CACHE  = fullfile(here, 'data', 'solver_diff_R500.mat');
    FS = 36;  R = 500e-6;

    if ~force && exist(CACHE,'file')
        S = load(CACHE);  fprintf('cache hit: %s\n', CACHE);
    else
        A  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
        MW = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
        sw(MW, A);                                                  % --- Maxwell：查詢點 ---
        cfgM = model_config('long2016_hexapole_halfcut','tip40um');
        rwM  = extract_maxwell_data(cfgM,'all',cfgM.default_variant);
        adM  = build_actuator_data(rwM,cfgM);
        [PM, BM] = cfgM.select_ball(adM, R);   NQ = size(PM,1);
        sw(A, MW);                                                  % --- APDL：內插源 ---
        cfgA = model_config('long2016_hexapole_halfcut','tip40um');
        rwA  = extract_ansys_data(cfgA,'all','graded');
        adA  = build_actuator_data(rwA,cfgA);
        sel  = adA.r2 < (620e-6)^2;                                 % 球外留一圈當內插源
        XA   = adA.Pa(sel,:);
        fprintf('查詢點 %d（Maxwell 格）；APDL 內插源節點 %d\n', NQ, sum(sel));
        Fi   = scatteredInterpolant(XA(:,1),XA(:,2),XA(:,3), adA.Ba(sel,1,1), 'linear','none');
        p2c  = zeros(1,6);  p2c(cfgA.apdl_to_paper_idx) = 1:6;      % paper P → APDL coil 欄
        S.nA = zeros(NQ,6);  S.nM = zeros(NQ,6);  S.vec = zeros(NQ,6);  S.ang = zeros(NQ,6);
        for p = 1:6
            bA = zeros(NQ,3);
            for c = 1:3, Fi.Values = adA.Ba(sel,c,p2c(p));  bA(:,c) = Fi(PM(:,1),PM(:,2),PM(:,3)); end
            bM = reshape(BM(:,p), 3, []).';
            nA = vecnorm(bA,2,2);  nM = vecnorm(bM,2,2);
            S.nA(:,p)  = nA;  S.nM(:,p) = nM;                       % 兩邊的場大小 [mT]（衍生量都由此算）
            S.vec(:,p) = vecnorm(bM - bA, 2, 2);                    % ‖ b_M − b_A ‖      [mT]
            S.ang(:,p) = acosd(min(1,max(-1, sum(bA.*bM,2)./(nA.*nM))));
        end
        S.r = vecnorm(PM,2,2);  S.R = R;  S.NQ = NQ;
        if ~exist(fullfile(here,'data'),'dir'), mkdir(fullfile(here,'data')); end
        save(CACHE,'-struct','S');  fprintf('saved %s\n',CACHE);
    end

    switch lower(WHICH)
        case 'mag',   v = (S.nM(:)./S.nA(:) - 1)*100;  xlab = 'Magnitude\;error\;(\%)';  ufmt = '%.2f %%';
        case 'magmt', v = abs(S.nM(:)-S.nA(:));        xlab = 'Magnitude\;error\;(mT)';  ufmt = '%.2f mT';
        case 'vec',   v = S.vec(:);   xlab = 'Field\;difference\;(mT)';    ufmt = '%.2f mT';
        case 'ang',   v = S.ang(:);   xlab = 'Direction\;error\;(degree)'; ufmt = '%.2f^{\\circ}';
        case 'angpct',v = sind(S.ang(:))*100;                   % 角度 → 橫向分量佔比 |b_\perp|/|b|
                      xlab = 'Direction\;error\;(\%)';          ufmt = '%.2f %%';
        otherwise,    error('WHICH 必為 ''mag'' | ''magmT'' | ''vec'' | ''ang'' | ''angpct''');
    end

    %% ---- 畫圖（風格①；圖例外置上方加框，同 err_hist_overlay）----
    nb = 180;  lo = prctile(v,1);  hi = prctile(v,99);       % [MODIFIED] 收窄顯示範圍（尾巴來自近極尖）
    edges = linspace(lo, hi, nb+1);
    % [MODIFIED] 直接以「百分比」當 BinCounts 畫，避免 probability(0–1) 與 %(0–100) 混用
    cnt = histcounts(v, edges);   pct = cnt/numel(v)*100;
    fig = figure('Color','w','Position',[100 100 1100 830]);   % 加高補償外置圖例
    ax  = axes(fig);  hold(ax,'on');
    col = [0.20 0.40 0.85];
    hb = histogram(ax, 'BinEdges',edges, 'BinCounts',pct, 'FaceColor',col, ...
                   'FaceAlpha',0.55, 'EdgeColor',[0.25 0.25 0.25], 'LineWidth',0.3);
    mu = mean(v);
    ml = xline(ax, mu, '--', 'Color','k', 'LineWidth',2.8);
    hold(ax,'off');
    ytop = ceil(max(pct)*1.05*10)/10;                       % 貼齊資料最大值 + 5% headroom（進位 0.1）
    yt   = round(linspace(0, ytop, 5), 2);  yt = yt(2:4);   % 3 個內縮 tick、首末不標
    [xr, xt] = ticks_x(lo, hi);
    xlim(ax,xr);  ylim(ax,[0 ytop]);
    set(ax,'XTick', xt(2:end-1), 'YTick', yt);   % [MODIFIED] 起訖點不進 XTick（框角不畫刻度線）
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.015 .015],'TickDir','out');
    box(ax,'on');  grid(ax,'off');
    % 起訖點用 text 補數字（有數字、無刻度線；figure-style 慣例 #5）
    for e = 1:2
        text(ax, xr(e), -0.022*ytop, sprintf('%g',xr(e)), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    xlabel(ax, ['$\mathbf{' xlab '}$'], 'Interpreter','latex', 'FontSize',FS);
    ylabel(ax, '$\mathbf{Percentage\;(\%)}$', 'Interpreter','latex', 'FontSize',FS);
    lg = legend([hb ml], {'Sampling range \leq 500 \mum', sprintf(['Mean = ' ufmt], mu)}, ...
                'Interpreter','tex', 'Location','northoutside', 'NumColumns',2);
    lg.FontSize = 24;  lg.FontWeight = 'bold';  lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;
    % [MODIFIED] 先把 axes 上緣壓低，再把圖例放到框外上方 → 不壓到圖框（同 err_hist_overlay 作法）
    drawnow;
    axp = get(ax,'Position');  lgh = get(lg,'Position');  lgh = lgh(4);  GAP = 0.022;
    axp(4) = axp(4) - (lgh + GAP);
    set(ax, 'Position', axp);
    set(lg, 'Position', [axp(1), axp(2)+axp(4)+GAP, axp(3), lgh]);
    ax.Toolbar.Visible = 'off';

    out = fullfile(figdir, sprintf('solver_diff_%s_hist_R500.png', lower(WHICH)));
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('wrote %s\n', out);  close(fig);
    fprintf('  N=%d  mean=%.3f  median=%.3f  p5=%.3f  p95=%.3f  (顯示 %.3f ~ %.3f)\n', ...
            numel(v), mu, median(v), prctile(v,5), prctile(v,95), lo, hi);
end

% ---- 水平軸：等距 nice tick、範圍貼齊資料（不硬拉回 0）、起訖點即為 tick ----
function [xr, tk] = ticks_x(lo, hi)
    nice = [0.05 0.1 0.2 0.25 0.5 1 2 2.5 5 10 20 25 50];
    for k = 1:numel(nice)
        s = nice(k);
        a = floor(lo/s)*s;  b = ceil(hi/s)*s;  n = round((b-a)/s)+1;
        if n>=5 && n<=8                                    % 5~8 個 tick（含兩端）
            xr = [a b];  tk = a:s:b;  return
        end
    end
    s = (hi-lo)/5;  xr = [lo hi];  tk = lo:s:hi;
end

% ---- 掛上 CAL、移除 OTHER（兩包有同名函式）----
function sw(CAL, OTHER)
    warning('off','MATLAB:rmpath:DirNotFound');
    rmpath(fullfile(OTHER,'function'));  rmpath(fullfile(OTHER,'common_path'));
    warning('on','MATLAB:rmpath:DirNotFound');
    addpath(fullfile(CAL,'function'));   addpath(fullfile(CAL,'common_path'));
end
