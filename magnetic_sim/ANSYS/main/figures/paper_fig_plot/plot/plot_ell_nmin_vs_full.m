function plot_ell_nmin_vs_full(QTY, RSEL, TOLPC, KWIN, force)
% plot_ell_nmin_vs_full -- paper 圖：減量取樣 vs 全取樣的 l_hat(R)，並算兩者的 RMSPE
% =========================================================================
%   兩條曲線（同一組 10 um 內插格、同一條校正管線、single 模型）：
%     ① **減量取樣**：每個 R 只取 N_min(R) 個隨機點（10 組排列取中位數）
%     ② **全取樣**  ：每個 R 取球內**全部**格點
%   縱軸 l_hat [um]、水平軸取樣半徑 R [um]；另算 RMSPE = rms((①-②)/②)。
%
%   N_min(R) 的定義（使用者拍板 2026-08-12）：對每組排列，找第一個 n 使
%   **其後 KWIN 點**的相對變化率都 <= TOLPC%；10 組取中位數。
%   ⚠ K=10 是必要的：K=5 會在震盪期產生偽陽性（n=42/57/67 這種不可能的值）。
%
%   🔑 **本圖不做任何新擬合**，全部從既有快取取值：
%     ① data/nmin_vs_R_maxwell10_single_random.mat（10 組 × 每個 n 的 l_hat）
%     ② data/ell_gain_sweep_maxwell_interp.mat（每個 R 用完整 10 um 池的 l_hat）
%     兩者是同一套網格與管線，可直接逐點相比。
%
%   風格 = 選項①粗體框 @ paper scale；圖例照 figure-style「圖例標準樣式」。
%   輸出 → figures/paper_fig/Section2_E/ell_nmin_vs_full_maxwell10_single.png
% =========================================================================
    clc;
    % [MODIFIED 2026-08-12] 起點拉回 R=100、視窗改 K=7（使用者拍板）。
    %   K=10 時 R=100 判不出來 —— 原因是**技術性的**：n 網格 log 間距、固定 45 點，池越小
    %   跨度越短，安定區（n>~1500）之後只剩 8.4 個格點，湊不出要檢查的 10 點。
    %   K=7 讓 21/21 個 R 都有解，且 10 組中位數**沒有一個**被偽陽（<200）污染。
    % [ADDED 2026-08-12] QTY='ell'（預設，從快取取值、零成本）| 'gain'（ĝ_I，需在
    %   那 21 個 (R, N_min) 組合上重跑擬合 —— 快取沒存 ĝ_I；結果另存 nmin_gain_*.mat）。
    if nargin < 1 || isempty(QTY),   QTY   = 'ell';      end
    if nargin < 2 || isempty(RSEL),  RSEL  = 100:20:500; end
    if nargin < 3 || isempty(TOLPC), TOLPC = 0.2;        end
    if nargin < 4 || isempty(KWIN),  KWIN  = 7;          end
    if nargin < 5 || isempty(force), force = false;      end
    % [ADDED 2026-08-12] QTY='nmae'：不比參數值，而是比**兩組參數各自把該 R 內全部點擬合得多好**。
    %   兩條曲線都在**同一個評估集（該 R 內全部格點）**上算 NMAE：
    %     ① 減量：用 N_min 個點校正出的參數 → 這是**泛化誤差**（訓練集 ⊂ 評估集）
    %     ② 全取樣：用全部點校正出的參數   → 這是**自身殘差**（訓練集 = 評估集，必為下界）
    %   NMAE 沿用 doc/error_definition 的**平均場正規化**（分母 = 平均場，L2 換成 L1）：
    %     NMAE = <||B_model - B_FEM||> / <||B_FEM||> * 100 [%]，對「全部點 x 6 個激發」取平均。
    switch lower(QTY)
        case 'ell',  ylab = '$\mathbf{\hat{\ell}\;(micro\;meter)}$';   fmt = '%8.3f';
        case 'gain', ylab = '$\mathbf{{}^{B}\hat{g}_{I}\;(mT/A)}$';    fmt = '%8.4f';
        case 'nmae', ylab = '$\mathbf{NMAE\;(\%)}$';                   fmt = '%8.4f';
        % [ADDED 2026-08-12] 只畫一條：|NMAE_全取樣 - NMAE_減量|（單位 = 百分點）。
        %   兩條 NMAE 曲線疊在一起看不出差別 → 用差值圖把「減量的代價」放大呈現。
        case 'nmaediff', ylab = '$\mathbf{|\Delta NMAE|\;(\%)}$';      fmt = '%8.4f';
        % [ADDED 2026-08-12] 每點平均 cost：J/N(R)。J 是專案既有定義（同 plot_npts_cost_vs_R）
        %   J = Σ_{j=1..6} ||S(l_hat,Pc)·g_j − b_j||² [mT²]，涵蓋全部點 x 3 分量 x 6 激發；
        %   除以 N(R) = 該 R 內的**點數**（不是 3N 或 18N），去掉 J 的外延性，兩個 R 才可比。
        case 'cost', ylab = '$\mathbf{Cost/N\;(mT^{2})}$';             fmt = '%10.3e';
        otherwise,   error('QTY 必為 ''ell'' | ''gain'' | ''nmae'' | ''nmaediff'' | ''cost''');
    end
    ONE = strcmpi(QTY,'nmaediff');                     % 單序列模式（無第二條曲線、無圖例）

    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    S = load(fullfile(here,'data','nmin_vs_R_maxwell10_single_random.mat'));   % 減量（10 組）
    F = load(fullfile(here,'data','ell_gain_sweep_maxwell_interp.mat'));       % 全取樣

    nR = numel(RSEL);
    ellA = nan(1,nR);  ellB = nan(1,nR);  Nuse = nan(1,nR);  pool = nan(1,nR);
    fprintf('   R    N_used   pool     佔比%%    ell_減量   ell_全取樣    誤差%%\n');
    for t = 1:nR
        i = find(S.RLIST == RSEL(t));
        assert(~isempty(i), '快取沒有 R=%d', RSEL(t));
        n = S.nall{i};   E = S.ellall{i};   pool(t) = S.pool(i);

        nd = nan(1, size(E,1));                       % 各組的 N_min
        for d = 1:size(E,1)
            nd(d) = flat_one(n, E(d,:), TOLPC, KWIN);
        end
        Nmed = median(nd, 'omitnan');
        [~, k] = min(abs(n - Nmed));                  % 取最接近中位數的 n 格
        Nuse(t) = n(k);
        ellA(t) = median(E(:,k), 'omitnan');          % 減量：該 n 的 10 組中位數
        ellB(t) = F.ell_R(F.Rum == RSEL(t));          % 全取樣
    end

    if strcmpi(QTY,'gain')                             % ĝ_I 需重跑擬合（見檔頭 [ADDED]）
        [ellA, ellB] = gain_curves(here, RSEL, Nuse, S, F, force);
    elseif strcmpi(QTY,'nmae') || ONE || strcmpi(QTY,'cost')
        met = 'nmae';  if strcmpi(QTY,'cost'), met = 'cost'; end
        [ellA, ellB] = nmae_curves(here, RSEL, Nuse, S, F, force, met);
        if ONE                                         % 差值圖：一條 |全取樣 - 減量|
            ellA = abs(ellB - ellA);   ellB = nan(size(ellA));
        end
    end
    if ONE                                             % [ADDED] 單序列：只印差值
        fprintf('\n   R    N_min     pool    佔比%%    |dNMAE| [百分點]\n');
        for t = 1:nR
            fprintf(['  %3d  %7d  %7d  %6.2f   ' fmt '\n'], ...
                    RSEL(t), Nuse(t), pool(t), Nuse(t)/pool(t)*100, ellA(t));
        end
        RMSPE = NaN;
        fprintf('\n  |dNMAE| 範圍 %.4f ~ %.4f 百分點（中位數 %.4f）；點數佔比 %.2f ~ %.2f %%\n', ...
                min(ellA), max(ellA), median(ellA), min(Nuse./pool)*100, max(Nuse./pool)*100);
    else
    for t = 1:nR
        fprintf(['  %3d  %7d  %7d  %6.2f   ' fmt '   ' fmt '   %+7.3f\n'], ...
                RSEL(t), Nuse(t), pool(t), Nuse(t)/pool(t)*100, ellA(t), ellB(t), ...
                (ellA(t)-ellB(t))/ellB(t)*100);
    end
    rel   = (ellA - ellB) ./ ellB;
    RMSPE = sqrt(mean(rel.^2)) * 100;
    fprintf('\n  RMSPE = %.4f %%   最大絕對誤差 %.3f %%   點數佔比 %.2f ~ %.2f %%\n', ...
            RMSPE, max(abs(rel))*100, min(Nuse./pool)*100, max(Nuse./pool)*100);
    end

    % ---- 畫圖（單 panel，選項①粗體框）----
    FS = 36;  LWBOX = 4.0;
    cA = [0.05 0.10 0.95];   cB = [0.85 0.10 0.10];
    fig = figure('Color','w','Position',[100 100 1160 830]);
    ax  = axes(fig);  hold(ax,'on');
    h1 = plot(ax, RSEL, ellA, '-o', 'Color',cA, 'LineWidth',3.0, 'MarkerSize',8, ...
              'MarkerFaceColor',cA, 'Clipping','off');
    if ~ONE
        h2 = plot(ax, RSEL, ellB, '--s','Color',cB, 'LineWidth',3.0, 'MarkerSize',8, ...
                  'MarkerFaceColor',cB, 'Clipping','off');
    end
    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX,'TickLength',[.015 .015]);
    ax.Toolbar.Visible = 'off';

    XL = [RSEL(1) RSEL(end)];   xlim(ax, XL);
    set(ax,'XTick', inner_ticks(XL(1), XL(2), 3));
    if ONE || strcmpi(QTY,'cost')
        % [ADDED] 差值與 cost/N 都是「距離 0 多遠」的量（有真實 0 下界、且資料貼近 0）
        %   → 基準線 0 貼下框、上緣只留 8% 裕度（figure-style「自 0 起的軸」；
        %   此時**不套**「兩端留白 = 間距」，否則下方會空掉一大截）。
        [YL, YT] = axlim_from_zero(max([ellA ellB], [], 'omitnan'), [5 3]);   % ONE 時 ellB 全 NaN
    else
        [YL, YT] = axlim_auto(min([ellA ellB]), max([ellA ellB]), [3 5]);
    end
    ylim(ax, YL);   set(ax,'YTick',YT);

    % [REMOVED 2026-08-12] 原本每點標 N_min 的數字已拿掉（使用者要求）——
    %   點數佔比改由獨立的火柴棒圖 plot_nmin_ratio_stem.m 呈現。
    yoff = YL(1) - 0.022*diff(YL);                    % 端點只標數字、不畫 tick
    for xv = XL
        text(ax, xv, yoff, sprintf('%g',xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top','FontSize',FS,'FontWeight','bold','Clipping','off');
    end
    xlabel(ax, '$\mathbf{Sampling\;range\;(micro\;meter)}$', 'Interpreter','latex', 'FontSize',FS);
    ylabel(ax, ylab, 'Interpreter','latex', 'FontSize',FS);

    % [MODIFIED 2026-08-12] QTY='nmae' 時兩條曲線本身就是誤差 → 圖例不再掛「RMSPE」
    %   （那是「兩條參數曲線的相對差」，對誤差曲線沒有意義）。
    % [MODIFIED 2026-08-12] 單序列（nmaediff）**不放圖例** —— 只有一條線、標籤已在軸標題裡，
    %   掛一個單則圖例只是多一個空框。此時軸也不必為圖例讓出上方空間。
    if ONE
        hold(ax,'off');  drawnow;
    else
        if strcmpi(QTY,'nmae') || strcmpi(QTY,'cost')   % 兩條本身就是誤差量 → 不掛 RMSPE
            lbls = {'Reduced sampling', 'Full sampling'};
        else
            lbls = {sprintf('Reduced sampling (RMSPE = %.2f%%)', RMSPE), 'Full sampling'};
        end
        lg = legend(ax, [h1 h2], lbls, ...
             'Interpreter','tex', 'Location','northoutside', 'NumColumns',2);
        lg.FontSize = 24;  lg.FontWeight = 'bold';
        lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;
        hold(ax,'off');  drawnow;
        axp = get(ax,'Position');  lgh = get(lg,'Position');  lgh = lgh(4);
        GAPN = 0.022;  newTop = 1 - lgh - GAPN - 0.006;
        axp(4) = newTop - axp(2);  set(ax,'Position',axp);
        set(lg, 'Position', [axp(1), newTop + GAPN, axp(3), lgh]);
    end

    out = fullfile(figdir, sprintf('%s_nmin_vs_full_maxwell10_single.png', lower(QTY)));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('\nwrote %s\n', out);
end

% ============================================================================
function [A, B] = nmae_curves(here, RSEL, Nuse, S, F, force, METRIC)
% [ADDED 2026-08-12] 兩組參數在**同一評估集（該 R 內全部格點）**上的場擬合 NMAE。
%   ① 減量 A(t)：用該 R 的 N_min 個點校正 → 預測到全部點（**泛化誤差**）；10 組取中位數
%   ② 全取樣 B(t)：用全部點校正 → 預測到全部點（**自身殘差**，必為 A 的下界）
%   🔑 **不重跑優化器**：兩邊的 l_hat 都已在既有快取（減量在 nmin cache 的 ellall、
%      全取樣在 ell_gain_sweep_maxwell_interp）。這裡只需用**各自的訓練點集**解一次
%      線性最小平方得電荷 G，再乘上「全部點的 kernel」預測 —— 全部是閉式解。
%   ⚠ 減量的 kernel 必須用該組自己的 l_hat 重建（不可借用全取樣的 l_hat）。
%   ⚠ 抽樣序與 N_min 研究一致：格建到 max(RSEL)、`sel=find(r<=R)`、`rng(d-1)` 後
%      `ord = sel(randperm(numel(sel)))` —— 同 seed 抽到同一組點。
    if nargin < 7 || isempty(METRIC), METRIC = 'nmae'; end
    cachef = fullfile(here,'data','nmae_nmin_vs_full_maxwell10_single.mat');
    if exist(cachef,'file') && ~force
        C = load(cachef);
        % [MODIFIED 2026-08-12] 舊版快取只有 A/B（NMAE）；缺 Ac/Bc（cost/N）就重算。
        if isequal(C.RSEL(:).',RSEL(:).') && isequal(C.Nuse(:).',Nuse(:).') ...
           && isfield(C,'Ac') && isfield(C,'Bc')
            fprintf('loaded cache %s\n', cachef);
            [A, B] = pick_metric(C, METRIC);   return;
        end
        fprintf('  [warn] 快取的 RSEL/Nuse 不同或缺 cost 欄位 → 重算\n');
    end

    solver_path();
    cfg = model_config('long2016_hexapole_halfcut','tip40um');
    raw = extract_maxwell_data(cfg, 'all', cfg.default_variant);
    ad  = build_actuator_data(raw, cfg);
    ad  = interp_grid_sample(ad, cfg, max(RSEL)*1e-6, 10e-6);
    rr  = sqrt(ad.r2);   Pc = ad.Pc_base;   N_I = size(ad.Ba,3);

    nR = numel(RSEL);   A = nan(1,nR);   B = nan(1,nR);
    Ac = nan(1,nR);     Bc = nan(1,nR);                 % [ADDED] cost/N（J/N，mT²）
    fprintf('\n   R    N_min     pool    NMAE_減量%%  NMAE_全取樣%%   cost/N_減量   cost/N_全取樣\n');
    for t = 1:nR
        sel = find(rr <= RSEL(t)*1e-6);   Ns = numel(sel);
        Pall = ad.Pa(sel,:);
        Ball = zeros(3*Ns, N_I);
        for j = 1:N_I, Ball(:,j) = reshape(ad.Ba(sel,:,j).', [], 1); end

        % ---- ② 全取樣：l_hat 取自 sweep 快取，G 用全部點解 ----
        lF   = F.ell_R(F.Rum == RSEL(t));
        SaF  = kernel_S(lF*1e-6, Pc, Pall);
        Gf   = (SaF.'*SaF) \ (SaF.'*Ball);
        B(t)  = nmae_val(SaF*Gf, Ball);
        Bc(t) = sum((SaF*Gf - Ball).^2, 'all') / Ns;    % J/N [mT²]

        % ---- ① 減量：每組用自己的 l_hat + 自己的 N_min 點解 G，再預測到全部點 ----
        i = find(S.RLIST == RSEL(t));   nl = S.nall{i};   E = S.ellall{i};
        [~, k] = min(abs(nl - Nuse(t)));                  % 對到 N_min 那個 n 格
        v = nan(1, size(E,1));   vc = nan(1, size(E,1));
        for d = 1:size(E,1)
            rng(d-1);   ord = sel(randperm(Ns));          % 與 N_min 研究同序
            idx = ord(1:Nuse(t));
            Ptr = ad.Pa(idx,:);
            Btr = zeros(3*Nuse(t), N_I);
            for j = 1:N_I, Btr(:,j) = reshape(ad.Ba(idx,:,j).', [], 1); end
            lD  = E(d,k)*1e-6;                            % 該組自己的 l_hat
            Str = kernel_S(lD, Pc, Ptr);                  % 訓練集 kernel
            Gd  = (Str.'*Str) \ (Str.'*Btr);              % 只用訓練集解電荷
            Sad = kernel_S(lD, Pc, Pall);                 % 同一 l_hat 的全點 kernel
            v(d)  = nmae_val(Sad*Gd, Ball);
            vc(d) = sum((Sad*Gd - Ball).^2, 'all') / Ns;  % J/N [mT²]
        end
        A(t)  = median(v,  'omitnan');
        Ac(t) = median(vc, 'omitnan');
        fprintf('  %3d  %7d  %7d   %9.4f   %10.4f   %11.4e   %12.4e\n', ...
                RSEL(t), Nuse(t), Ns, A(t), B(t), Ac(t), Bc(t));
    end
    save(cachef, 'RSEL', 'Nuse', 'A', 'B', 'Ac', 'Bc');
    fprintf('saved cache %s\n', cachef);
    [A, B] = pick_metric(struct('A',A,'B',B,'Ac',Ac,'Bc',Bc), METRIC);
end

% ============================================================================
function [A, B] = pick_metric(C, METRIC)
    if strcmpi(METRIC,'cost'), A = C.Ac;  B = C.Bc;   else, A = C.A;  B = C.B;  end
end

% ============================================================================
function [lim, tk] = axlim_from_zero(maxv, nlist) %#ok<INUSD>
% [ADDED 2026-08-12] 與 plot_nmin_ratio_stem.m 同一份實作（figure-style「自 0 起的軸」）：
%   下緣固定 0（基準線貼下框）、上緣**只留 8% 裕度**、刻度 (1:n)*s 取不超出上緣的最大奇數 n。
%   ⚠ 此軸**不套**「兩端留白 = 間距」（那條只適用兩端都無界的軸）。
    cand = [1 2 2.5 3 4 5 10];
    x = maxv/4;   k = floor(log10(x));
    s = cand(find(cand*10^k >= x, 1)) * 10^k;
    top = 1.08 * maxv;
    n = floor((top - 1e-12)/s);
    if mod(n,2) == 0, n = n - 1; end
    n = max(n, 1);
    lim = [0, top];   tk = (1:n)*s;
end

% ============================================================================
function Sm = kernel_S(l_hat, Pc, P)
% 無因次點電荷 kernel（3Np x 6），與 solve_current / fitting 的 build_S 同式。
    Np   = size(P,1);
    pbar = P / l_hat;
    Sm   = zeros(3*Np, 6);
    for k = 1:6
        dd = pbar - Pc(:,k).';
        r3 = sum(dd.^2, 2).^1.5;
        Sm(:,k) = reshape((dd ./ r3).', 3*Np, 1);
    end
end

% ============================================================================
function v = nmae_val(Bm, Bt)
% NMAE [%]：對「全部點 x 6 個激發」的**向量誤差範數**取平均，除以平均場範數。
%   = <||B_model - B_FEM||> / <||B_FEM||> * 100   （doc/error_definition 的平均場正規化）
    dm = reshape(Bm - Bt, 3, []);      % 3 x (Np*6)
    dt = reshape(Bt,      3, []);
    v  = mean(vecnorm(dm,2,1)) / mean(vecnorm(dt,2,1)) * 100;
end

% ============================================================================
function [gA, gB] = gain_curves(here, RSEL, Nuse, S, F, force)
% ĝ_I 的兩條曲線：減量（每個 R 在 n=Nuse 抽 M 組、取中位數）與全取樣（既有 sweep）。
%   ⚠ 這是唯一需要重跑擬合的部分 —— 快取只存了 l_hat，沒存 ĝ_I。
%   點數都很小（222~2807），21 個 R × M 組，成本遠低於原本的 n-sweep。
    gcache = fullfile(here, 'data', 'nmin_gain_maxwell10_single.mat');
    gB = arrayfun(@(R) F.gI_R(F.Rum == R), RSEL);           % 全取樣（既有）
    if exist(gcache,'file') && ~force
        C = load(gcache);
        if isequal(C.RSEL(:).', RSEL(:).') && isequal(C.Nuse(:).', Nuse(:).')
            gA = C.gA;   fprintf('loaded gain cache %s\n', gcache);   return;
        end
        fprintf('  [warn] gain 快取的 (R, N_min) 與本次不同 → 重算\n');
    end

    solver_path();
    cfg = model_config('long2016_hexapole_halfcut','tip40um');
    raw = extract_maxwell_data(cfg, 'all', cfg.default_variant);
    ad  = build_actuator_data(raw, cfg);
    ad  = interp_grid_sample(ad, cfg, max(RSEL)*1e-6, 10e-6);
    rr  = sqrt(ad.r2);
    Fm  = zeros(6, cfg.N_I);
    for j = 1:cfg.N_I, Fm(cfg.apdl_to_paper_idx(j), j) = 1; end

    MREP = size(S.ellall{1}, 1);
    gA = nan(1, numel(RSEL));
    fprintf('\n  [gain] 逐 R 重跑（%d 組取中位數）\n', MREP);
    for t = 1:numel(RSEL)
        sel = find(rr <= RSEL(t)*1e-6);   np = numel(sel);
        gd  = nan(1, MREP);
        for d = 1:MREP
            rng(d-1);   ord = sel(randperm(np));            % 與 l_hat 用同一組 seed/排列
            idx = ord(1:Nuse(t));
            P   = ad.Pa(idx,:);
            Bs  = zeros(3*Nuse(t), cfg.N_I);
            for j = 1:cfg.N_I
                Bs(:,j) = reshape(ad.Ba(idx,:,j).', [], 1);
            end
            [e, l_hat] = fitting(P, Bs, ad.Pc_base, 0.5e-3, false);
            [~, g]     = solve_current(l_hat, e, ad.Pc_base, P, Bs, Fm);
            gd(d) = g;
        end
        gA(t) = median(gd, 'omitnan');
        fprintf('   R=%3d  n=%5d  gI_med=%.4f  [%.4f ~ %.4f]\n', RSEL(t), Nuse(t), gA(t), min(gd), max(gd));
    end
    save(gcache, 'RSEL','Nuse','gA');
    fprintf('saved gain cache %s\n', gcache);
end

% ============================================================================
function solver_path()
    APDL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    CAL  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    warning('off','MATLAB:rmpath:DirNotFound');
    rmpath(fullfile(APDL,'function'));  rmpath(fullfile(APDL,'common_path'));
    warning('on','MATLAB:rmpath:DirNotFound');
    addpath(fullfile(CAL,'function'));  addpath(fullfile(CAL,'common_path'));
end

% ============================================================================
function adq = interp_grid_sample(ad, cfg, Rmax, h)
    g = -Rmax:h:Rmax;
    [X, Y, Z] = ndgrid(g, g, g);
    Pq = [X(:), Y(:), Z(:)];
    rq = vecnorm(Pq, 2, 2);
    in = rq <= Rmax;   Pq = Pq(in,:);   rq = rq(in);
    Pm  = (ad.R_act.' * Pq.').';
    air = filter_iron_nodes(Pm(:,1), Pm(:,2), Pm(:,3) + cfg.SPH_OFST, cfg);
    Pq  = Pq(air,:);   rq = rq(air);
    src = ad.r2 < (Rmax + 60e-6)^2;
    Xs  = ad.Pa(src,:);
    N_I = size(ad.Ba,3);
    Bq  = zeros(size(Pq,1), 3, N_I);
    Fi  = scatteredInterpolant(Xs(:,1), Xs(:,2), Xs(:,3), ad.Ba(src,1,1), 'linear', 'none');
    for j = 1:N_I
        for c = 1:3
            Fi.Values = ad.Ba(src,c,j);
            Bq(:,c,j) = Fi(Pq);
        end
    end
    ok = ~any(any(isnan(Bq),3),2);
    fprintf('  [interp] 格距 %.0f um：球內 %d 點 → 濾鐵後 %d → 有效 %d\n', ...
            h*1e6, nnz(in), size(Pq,1), nnz(ok));
    adq = ad;
    adq.Pa = Pq(ok,:);   adq.r2 = rq(ok).^2;   adq.Ba = Bq(ok,:,:);
end

% ============================================================================
function k = flat_one(nlist, ell, TOLPC, KWIN)
% 單一條 l_hat(n) 的平坦段起點：第一個 a，使其後 KWIN 點的相對變化率都 <= TOLPC%。
    dch = [NaN, abs(diff(ell))./ell(2:end)*100];
    k = NaN;
    for a = 1:numel(nlist)-KWIN
        if all(dch(a+1 : a+KWIN) <= TOLPC), k = nlist(a);  return; end
    end
end

% ============================================================================
function tk = inner_ticks(lo, hi, n)
    cand = [1 2 2.5 3 4 5 10];
    s = (hi-lo)/(n+1);   k = floor(log10(s));
    [~, i] = min(abs(cand*10^k - s));   s = cand(i)*10^k;
    ctr = round(((lo+hi)/2)/s)*s;
    tk  = ctr + (-(n-1)/2 : (n-1)/2)*s;
    tk  = tk(tk > lo & tk < hi);
end

% ============================================================================
function [lim, tk] = axlim_auto(lo, hi, nlist)
% 奇數個等距 tick，兩端留白 = tick 間距（figure-style）。
    cand = [1 2 2.5 3 4 5 10];
    mid  = (lo+hi)/2;   rng_ = max(hi-lo, realmin);
    best = {};   bestSpan = inf;
    for n = nlist
        k0 = floor(log10(rng_/(n+1)));
        for k = k0:(k0+4)
            hit = false;
            for c = cand
                s = c*10^k;
                t = round(mid/s)*s + (-(n-1)/2 : (n-1)/2)*s;
                L = [t(1)-s, t(end)+s];   clr = 0.15*s;
                if lo >= L(1)+clr && hi <= L(2)-clr
                    if (n+1)*s < bestSpan, bestSpan = (n+1)*s;  best = {L, t}; end
                    hit = true;  break;
                end
            end
            if hit, break; end
        end
    end
    if isempty(best)
        n = nlist(1);  s = rng_/(n+1);
        t = mid + (-(n-1)/2 : (n-1)/2)*s;   best = {[t(1)-s, t(end)+s], t};
    end
    lim = best{1};   tk = best{2};
end
