function plot_gain_iso_hung_ws(R_EVAL_um, XLOG)
% PLOT_GAIN_ISO_HUNG_WS  純 hung 的 gain(𝒞)/iso(κ) 疊圖：工作空間半徑 R700 vs R300。
%   USE_BIAS（18-param）校正、全 no-gap。charge model 平滑函數上均勻取樣 𝒞/κ（非 FEM 節點，避免密度偏差）。
%   讀各變體 bias fit（fit_bias_R150um_<variant>.mat 的 ell/gB/Khat/Pc；Pc=含 bias 的電荷 lattice）——
%   **擬合一律 R ≤ 150µm（兩變體同條件）**，本腳本只決定「評估/取樣」半徑。
%   （R500=no_gap baseline 因 .dat 跨 coil 網格不一致暫略，見 project 記憶。）
%
%   [ADDED] R_EVAL_um = 取樣球半徑 [µm]，兩種模式：
%     省略 / []  → 各變體用「自己的工作半徑」（R700→700µm、R300→300µm）＝**外推**出擬合區
%                   （外推倍率不對等：R700 4.7× vs R300 2×，解讀差異時要留意）
%                   → gain_hung_ws.png / iso_hung_ws.png
%     給數字(150) → 兩變體都用該半徑 ＝ R150 版**擬合區內、不外推、同區域可直接比**
%                   → gain_hung_ws_R150.png / iso_hung_ws_R150.png
%   點數只由 ball_grid 的 nr=32 決定（137,376 點），與半徑無關。
%   風格：選項① 粗體框、共用 bins(nb=180)、percentage 縱軸、mean 虛線、legend 含 mean、headroom。
%
%   [ADDED] XLOG（僅作用於 gain 𝒞 圖；κ 是 0~1 無因次、兩變體範圍重疊，維持 linear）：
%     false/省略 → linear bins（固定絕對寬 (max-min)/180）→ gain_hung_ws<sfx>.png
%     true       → log 等距 bins（固定相對寬 (max/min)^(1/180)）+ XScale log → gain_hung_ws<sfx>_log.png
%     為何：R700(𝒞̄884) 與 R300(𝒞̄1.29e4) 相距 14.6×，共用 linear edges 時 bin 寬 79.4
%           → R700 整條只佔 ~4 bin（單根 >70%）。log 把解析度依相對寬度重分配
%           （884 附近 bin 寬 ~15、12910 附近 ~222）→ R700 攤成 ~22 bin。bin 總數不變。

    if nargin < 1, R_EVAL_um = []; end                                      % [ADDED] 預設 = 各自工作半徑（現值行為）
    if nargin < 2, XLOG = false; end                                        % [ADDED] 預設 linear（現值行為）

    HUNG_CAL = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
                'hung_hexapole\Calibration_using_FEM_modeling'];

    % ---- 變體：variant 夾 / legend / 工作半徑 / 顏色 ----
    variants = {'R700', 'R300'};
    labels   = {'R700', 'R300'};
    cols     = [0.12 0.18 0.55;    % R700 navy
                0.80 0.10 0.10];   % R300 red

    % ---- [ADDED] 取樣半徑 + 檔名後綴（依 R_EVAL_um 分支）----
    if isempty(R_EVAL_um)
        Rwork = [700, 300] * 1e-6;                                         % 各自工作半徑取樣球（外推）
        sfx   = '';
    else
        Rwork = repmat(R_EVAL_um*1e-6, 1, numel(variants));                % 兩變體同半徑（擬合區內）
        sfx   = sprintf('_R%d', R_EVAL_um);
    end

    figdir  = fullfile(HUNG_CAL,'current_base','figures','eighteen_param');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    datadir = fullfile(HUNG_CAL,'current_base','data');

    % ---- 逐變體：載 bias fit → 取樣球內均勻取樣 𝒞/κ ----
    Cc = cell(1,numel(variants));  Kk = cell(1,numel(variants));
    for k = 1:numel(variants)
        matf = fullfile(datadir, sprintf('fit_bias_R150um_%s.mat', variants{k}));
        assert(exist(matf,'file')==2, '缺 %s（先跑 USE_BIAS main.m，VARIANT=%s）', matf, variants{k});
        S = load(matf, 'ell','gB','Khat','Pc');
        ell_m = S.ell*1e-6;  Hhat = S.gB*S.Khat;                          % ell 存 µm；Hhat=ĝ·K̄ [mT/A]
        P = ball_grid(Rwork(k), 32);
        [Cc{k}, Kk{k}] = ck_nodes(P, ell_m, Hhat, S.Pc);                  % 用該 fit 的電荷 lattice Pc（含 bias）
        fprintf('%s (R=%dµm, N=%d)  meanC=%.4g  meanK=%.4f\n', labels{k}, round(Rwork(k)*1e6), size(P,1), mean(Cc{k}), mean(Kk{k}));
    end

    % gain：可選 log 等距 bins + log x 軸（檔名加 _log）；iso：一律 linear
    lsfx = ''; if XLOG, lsfx = '_log'; end                                 % [ADDED]
    render_overlayN(Cc, '$\mathcal{C}\;[(\mathrm{mT/A})^{3}]$', cols, labels, fullfile(figdir,['gain_hung_ws' sfx lsfx '.png']), XLOG);
    render_overlayN(Kk, '$\kappa$',                            cols, labels, fullfile(figdir,['iso_hung_ws'  sfx '.png']), false);
end

% ============================================================================
function P = ball_grid(R, nr)
    h = R/nr;  v = -R:h:R;
    [X,Y,Z] = ndgrid(v,v,v);
    in = (X.^2 + Y.^2 + Z.^2) <= R^2;
    P = [X(in), Y(in), Z(in)];
end

function [Cv,Kv] = ck_nodes(P, ell_m, Hhat, Pc)
    n = size(P,1);  Cv = zeros(n,1);  Kv = zeros(n,1);
    for i = 1:n
        p  = P(i,:).';
        Dk = p/ell_m - Pc;                                    % 3×6
        sv = svd((Dk ./ (vecnorm(Dk).^3)) * Hhat);
        Cv(i) = prod(sv);  Kv(i) = sv(3)/sv(1);
    end
end

function render_overlayN(Y, xlab, cols, labels, outpng, xlog)
% N 條疊圖（共用 bins、percentage 縱軸、mean 虛線、legend 含 mean + CV、選項① 粗體框）。
%   xlog=true → log 等距 edges + XScale log（解析度依相對寬度分配，bin 數不變）。
    if nargin < 6, xlog = false; end
    nb = 180;
    allY = cat(1, Y{:});
    if xlog
        edges = logspace(log10(min(allY)), log10(max(allY)), nb+1);         % [ADDED] 固定比：(max/min)^(1/nb)
    else
        edges = linspace(min(allY), max(allY), nb+1);                       % 固定差：(max-min)/nb
    end
    fig = figure('Color','w','Position',[100 100 820 600]);  ax = axes(fig);  hold(ax,'on');
    h = gobjects(1,numel(Y));  leg = cell(1,numel(Y));
    for k = 1:numel(Y)
        h(k) = histogram(ax, Y{k}, edges, 'Normalization','percentage', ...   % [MODIFIED] 縱軸改 percentage
                         'FaceColor',cols(k,:), 'FaceAlpha',0.55, 'EdgeColor','w', 'LineWidth',0.3);
        % [MODIFIED] legend 加 CV = std/mean（相對離散度，跨變體可直接比）
        leg{k} = sprintf('%s  (mean = %.4g,  CV = %.2f%%)', labels{k}, mean(Y{k}), 100*std(Y{k})/mean(Y{k}));
    end
    for k = 1:numel(Y)
        xline(ax, mean(Y{k}), '--', 'Color',cols(k,:), 'LineWidth',2);
    end
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    box(ax,'on'); grid(ax,'off');
    if xlog                                                                 % [ADDED]
        set(ax,'XScale','log');  xlim(ax,[edges(1) edges(end)]);
    end
    yl = ylim(ax); ylim(ax,[0 yl(2)*1.20]);
    yt = get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    xint = 'tex'; if startsWith(strtrim(xlab),'$'), xint = 'latex'; end
    xlabel(ax, xlab, 'FontWeight','bold','Interpreter',xint);
    ylabel(ax,'Percentage (%)','FontWeight','bold');                          % [MODIFIED]
    lg = legend(ax, h, leg, 'Location','northeast');  set(lg,'FontSize',14,'FontWeight','bold','Box','on');
    ax.Toolbar.Visible = 'off';
    exportgraphics(fig, outpng, 'Resolution',600);
    fprintf('saved %s\n', outpng);
    close(fig);
end
