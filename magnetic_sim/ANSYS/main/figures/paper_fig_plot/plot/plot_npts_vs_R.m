function plot_npts_vs_R(SRC, SCALE, MODEL, XVAR, REF, force)
% plot_npts_vs_R -- 取樣球內節點數 N 隨取樣半徑 R 的變化（累積計數）
% =========================================================================
%   N(R) = 「|p| <= R 的節點數」，p 為 actuator frame 的節點位置（已濾鐵）。
%
%   ⚠ **這張圖不涉及球殼分箱**：N(R) 是累積量，每個 R 都是把全部節點與 R 比一次大小的
%     精確計數，沒有分箱、沒有統計雜訊。曲線唯一的參數是取點間距 DR（預設 2 µm），
%     它只決定曲線解析度，不影響任何一點的數值。
%     （前面分析用的「球殼厚度 W」只在反推**局部密度** rho(r) 時才需要 —— 那要做差分，
%       才有薄->雜訊 / 厚->抹平梯度 的取捨。累積量沒有這個問題。）
%
%   資料源：SRC='maxwell'（預設）讀 WP 細格 .fld（步距 0.02 mm）；'apdl' 讀 graded .dat。
%   前段沿用校正管線：extract -> build_actuator_data（轉 actuator frame + 濾鐵），
%   與 select_ball 取點的判準完全一致，所以 N(R) 就是校正時實際會拿到的點數。
%
%   節點半徑快取在 data/npts_radii_<src>.mat（重跑只需數秒；重讀 .fld 要數十秒）。
%
%   風格①粗體框圖：字級 36、框線 4.0、box on、無 grid、單位括號、
%   刻度奇數個等距、水平軸端點只標數字不畫 tick、首末資料點貼齊框邊。
%
%   輸出 → figures/paper_fig/Section2_E/npts_vs_R_<src>.png（覆蓋迭代）
% =========================================================================
    clc;
    if nargin < 1 || isempty(SRC),   SRC   = 'maxwell'; end
    if nargin < 2 || isempty(SCALE), SCALE = 'lin';     end   % 'lin' | 'ylog'
    if nargin < 3,                   MODEL = [];        end   % [a b]：疊上閉式解 N(R;a,b)
    if nargin < 4 || isempty(XVAR),  XVAR  = 'R';       end   % 'R'（橫軸 R）| 'R3'（橫軸 R³）
    if nargin < 5 || isempty(REF),   REF   = false;     end   % true = 疊過原點參考線 + 算 INL
    if nargin < 6 || isempty(force), force = false;     end

    DR   = 2;                                   % 曲線取點間距 [µm]（只影響解析度）
    RMAX = 500;                                 % 繪圖上限 [µm]

    here   = fileparts(fileparts(mfilename('fullpath')));      % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    rr = load_radii(SRC, here, force);           % 全部節點的 |p| [µm]

    R = (0:DR:RMAX).';
    N = arrayfun(@(x) nnz(rr <= x), R);          % 逐點精確計數（定義無歧義）

    fprintf('%s: 節點總數 %d｜N(150)=%d  N(300)=%d  N(500)=%d\n', ...
            upper(SRC), numel(rr), N(R==150), N(R==300), N(R==500));

    Nm = [];
    if ~isempty(MODEL)
        a = MODEL(1);  b = MODEL(2);
        Nm = model_N(R, a, b);
        q  = R >= 40 & N > 0;                    % 比較區（避開內側點數過少）
        e  = Nm(q)./N(q) - 1;
        fprintf('  模型 a=%.4f µm, b=%.6f：R>=40 µm 相對誤差 %+.1f%% ~ %+.1f%%（RMS %.1f%%）\n', ...
                a, b, 100*min(e), 100*max(e), 100*sqrt(mean(e.^2)));
        for RR = [100 150 200 300 400 500]
            j = find(R==RR,1);
            fprintf('    R=%3d: 實測 %6d  模型 %6.0f  (%+.1f%%)\n', RR, N(j), Nm(j), 100*(Nm(j)/N(j)-1));
        end
    end
    % ---- 過原點參考線 + INL（積分非線性）----
    %   參考線 N = k·R³（強制過原點：N(0)=0 是物理必然，不讓擬合去湊截距）
    %   INL = max|N − k·R³| / FS × 100%，FS = N(R_max)
    ref = [];
    if REF
        x3 = R.^3;   k = (x3.'*N)/(x3.'*x3);
        ref = k*x3;
        dev = N - ref;  [mx, im] = max(abs(dev));
        ref = struct('y', ref, 'k', k, 'INL', 100*mx/N(end), 'Rmax', R(im), 'dev', mx);
        fprintf('  過原點參考線：k = %.4e 點/µm³｜INL = %.2f%% FS（最大偏差 %.0f 點 @ R = %g µm）\n', ...
                k, ref.INL, mx, R(im));
    end
    render(R, N, Nm, MODEL, ref, SRC, DR, SCALE, XVAR, figdir);
end

% ============================================================================
function N = model_N(R, a, b)
% 閉式解 N(R) = (4/3)pi (R/a)^3 * Phi(x)，x = bR/a
%   Phi(x) = (3/x^3)[ ln(1+x) + 2/(1+x) - 1/(2(1+x)^2) - 3/2 ]，Phi(0)=1
%   小 x 用級數避免 0/0 抵消（分子是三個 O(1) 項相消到 O(x^3)）。
%   級數（展開到 x^5）：ln(1+x)+2/(1+x)-1/(2(1+x)^2)-3/2 = x^3/3 - (3/4)x^4 + (6/5)x^5 - ...
%   => Phi(x) = 1 - (9/4)x + (18/5)x^2 - ...
    x   = b*R/a;
    Phi = 1 - 2.25*x + 3.6*x.^2;                  % 小 x 用級數（避免 0/0 抵消）
    big = x > 1e-3;                               % x 夠大才用閉式（雙精度已穩）
    u   = 1 + x(big);
    Phi(big) = (3./x(big).^3) .* ( log(u) + 2./u - 1./(2*u.^2) - 1.5 );
    N = (4/3)*pi*(R/a).^3 .* Phi;
end

% ============================================================================
function rr = load_radii(SRC, here, force)
% 取全部節點在 actuator frame 的 |p|（已濾鐵），快取成小 .mat
    cf = fullfile(here, 'data', sprintf('npts_radii_%s.mat', lower(SRC)));
    if exist(cf,'file') && ~force
        S = load(cf);  rr = S.rr;
        fprintf('由快取載入 %s（%d 節點）\n', cf, numel(rr));
        return
    end
    if strcmpi(SRC,'maxwell')
        ROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
        addpath(fullfile(ROOT,'function'));
        cfg = model_config('long2016_hexapole_halfcut','tip40um');
        raw = extract_maxwell_data(cfg, 'all', cfg.default_variant);
    else
        ROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
        addpath(fullfile(ROOT,'function'), fullfile(ROOT,'common_path'), fullfile(ROOT,'utils'));
        cfg = model_config('long2016_hexapole_halfcut','tip40um');
        raw = extract_ansys_data(cfg, 'all', 'graded');
    end
    ad = build_actuator_data(raw, cfg);          % actuator frame + 濾鐵
    rr = vecnorm(ad.Pa, 2, 2)*1e6;               % m → µm
    save(cf, 'rr');
    fprintf('已存 %s（%d 節點）\n', cf, numel(rr));
end

% ============================================================================
function s = nice_step(x)
% 取「>= x 的最小 nice 步長」（figure-style：等距、疏密適中）
    n = [1 2 2.5 3 4 5 10];
    e = floor(log10(x));   c = x/10^e;
    i = find(n >= c - 1e-12, 1);   s = n(i)*10^e;
end

% ============================================================================
function render(R, N, Nm, MODEL, ref, SRC, DR, SCALE, XVAR, figdir)
    FS = 36;  LWBOX = 4.0;
    fig = figure('Color','w','Position',[80 60 1250 950]);  ax = axes(fig);  hold(ax,'on');

    r3 = strcmpi(XVAR,'R3');
    if r3, xv = (R.^3)/1e7;  xlab = '$\mathbf{R^{3}\;(\times 10^{7}\,\mu m^{3})}$';
    else,  xv = R;           xlab = '$\mathbf{R\;(\mu m)}$';
    end

    ylog = ~strcmpi(SCALE,'lin');
    if ylog
        q = N >= 1;                              % log 軸不能畫 N = 0
        h1 = plot(ax, xv(q), N(q), '-', 'Color',[0.10 0.35 1.00], 'LineWidth',4.0, 'Clipping','off');
        if ~isempty(Nm), h2 = plot(ax, xv(q), Nm(q), '--', 'Color',[0.85 0.10 0.10], 'LineWidth',3.5); end
    else
        h1 = plot(ax, xv, N/1e3, '-', 'Color',[0.10 0.35 1.00], 'LineWidth',4.0, 'Clipping','off');
        if ~isempty(Nm), h2 = plot(ax, xv, Nm/1e3, '--', 'Color',[0.85 0.10 0.10], 'LineWidth',3.5); end
        if ~isempty(ref), h3 = plot(ax, xv, ref.y/1e3, '--', 'Color',[0.20 0.20 0.20], 'LineWidth',3.0); end
    end

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
           'TickLength',[.015 .015],'TickDir','out');

    XL = [xv(1) xv(end)];                        % 首末資料點貼齊左右框邊
    xlim(ax, XL);
    if r3, set(ax,'XTick',[3 6 9]);  else, set(ax,'XTick',[125 250 375]);  end   % 3 個等距內部刻度
    if ylog
        set(ax,'YScale','log');                             % 以 10 為底
        YL = [1e-1 1e5];  ylim(ax, YL);  set(ax,'YTick',10.^(0:4));   % 5 個等距(每格 1 decade)
        set(ax,'YMinorTick','off');                         % 關掉 log 軸預設的密集次刻度
        fprintf('  XL = [%g %g]; YL = [%g %g]（log10，YTick 10^0~10^4）；第一個 N>=1 在 R = %g µm\n', ...
                XL, YL, R(find(q,1)));
    else
        % [MODIFIED 2026-08-11] 刻度數量在 {5,3} 中取「填充率最高」者 —— 原本寫死 3 個,
        %   nice_step 進位後 YL 上緣常遠高於資料(本例 max 511.6 對 YL=800,上方空掉 36%)。
        %   兩者都守 figure-style:奇數個、等距、兩端留白 = 間距(YL = [0,(n+1)·ny])。
        nyc = [5 3];   fill = zeros(size(nyc));   nys = zeros(size(nyc));
        for k = 1:numel(nyc)
            nys(k)  = nice_step(max(N)/1e3/(nyc(k)+1));
            fill(k) = (max(N)/1e3) / ((nyc(k)+1)*nys(k));
        end
        [~, kb] = max(fill);   nt = nyc(kb);   ny = nys(kb);
        YL = [0 (nt+1)*ny];  ylim(ax, YL);  set(ax,'YTick',(1:nt)*ny);
        fprintf('  XL = [%g %g]; YL = [%g %g]×10³（%d 個 YTick,間距 %g,填充率 %.1f%%）\n', ...
                XL, YL, nt, ny, fill(kb)*100);
    end

    % 水平軸端點只標數字、不畫 tick（log 軸要在 log 空間下算偏移，線性公式會跑到負值）
    if ylog, yoff = 10^(log10(YL(1)) - 0.030*(log10(YL(2))-log10(YL(1))));
    else,    yoff = YL(1) - 0.025*diff(YL);
    end
    for j = 1:2
        text(ax, XL(j), yoff, sprintf('%g',XL(j)), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top','FontSize',FS,'FontWeight','bold','Clipping','off');
    end

    xlabel(ax, xlab, 'Interpreter','latex', 'FontSize',FS);
    if ylog, ylabel(ax, '$\mathbf{N}$', 'Interpreter','latex', 'FontSize',FS);
    else,    ylabel(ax, '$\mathbf{N\;(\times 10^{3})}$', 'Interpreter','latex', 'FontSize',FS);
    end
    % 圖例：照 figure-style「圖例標準樣式」（northoutside + 對齊框寬 + 黑粗框 + FS24）
    % ⚠ Interpreter 必須用 'tex'：用 'latex' 會讓 LaTeX 接管字型，FontWeight='bold' 被忽略
    ltxt = {'FEM data'};   hh = h1;   nc = 1;
    if ~isempty(Nm)
        ltxt{end+1} = sprintf('Model: a = %.4f \\mum, b = %.6f', MODEL(1), MODEL(2));
        hh(end+1) = h2;   nc = 2;
    end
    if ~isempty(ref)
        ltxt{end+1} = sprintf('Linear ref. (Linearity = %.2f%%)', 100 - ref.INL);
        hh(end+1) = h3;   nc = 2;
    end
    lg = legend(hh, ltxt, 'Interpreter','tex', 'Location','northoutside', 'NumColumns',nc);
    lg.FontSize = 24;  lg.FontWeight = 'bold';
    lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;
    drawnow;
    axp = get(ax,'Position');  lgp = get(lg,'Position');
    GAPN = 0.022;  newTop = 1 - lgp(4) - GAPN - 0.006;
    axp(4) = newTop - axp(2);  set(ax,'Position',axp);
    if nc >= 2      % 多則圖例：寬度切齊座標框（figure-style 標準）
        set(lg, 'Position', [axp(1), newTop + GAPN, axp(3), lgp(4)]);
    else            % 單則圖例：寬度隨內容、置中；拉滿框寬會留一大片空白
        set(lg, 'Position', [axp(1) + (axp(3)-lgp(3))/2, newTop + GAPN, lgp(3), lgp(4)]);
    end
    ax.Toolbar.Visible = 'off';  hold(ax,'off');

    sfx = '';  if r3, sfx = '3'; end
    out = fullfile(figdir, sprintf('npts_vs_R%s_%s.png', sfx, lower(SRC)));
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('wrote %s  (DR = %g µm)\n', out, DR);
end
