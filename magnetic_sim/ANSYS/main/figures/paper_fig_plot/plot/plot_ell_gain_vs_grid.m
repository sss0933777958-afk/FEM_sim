function plot_ell_gain_vs_grid(QTY, DESIGN, force)
% plot_ell_gain_vs_grid -- 校正參數隨「取樣點數」的變化（等測度網格族）
% =========================================================================
%   取樣半徑固定 R = 150 um，取樣點由 sphere_grid_sample 產生（格心取樣、三線性內插取場）。
%   四種掃描設計（DESIGN）：
%
%   'ratio'   三方向等比例放大，維持配比 N_r:N_phi:N_theta = 1:3:3pi
%             (1,3,9)…(5,15,47) → N = 27, 228, 756, 1824, 3525。水平軸 log。
%   'radial'  固定角向 (15,47)、只掃徑向 N_r=1..12 → N = 705..8460。水平軸 linear。
%             結論：徑向在 N_r>=2 之後完全飽和（中心空洞 119→52 um 對結果無影響）。
%   'angular' 固定徑向 N_r=2、只掃角向 (1,3)…(32,101) → N = 6..6464。水平軸 log。
%             結論：l_hat 與 g_I 都由**角向**決定；g_I 偏差 ~ N^-0.31 收斂極慢。
%   'inner'   **看小點數震盪**（使用者 2026-08-12 指定）：在 (5,15,47) 的最內殼
%             (r1 = 69.6 um) 以內自訂內層「環」，半徑 r = 5,10,…,50 um、環上點數
%             n = 5,10,…,50（一一對應），**由內往外累積** → N = 5,15,30,50,75,105,
%             140,180,225,275；最後一點再疊上 (5,15,47) 外層 → N = 3800。
%             另有 **N=0 錨點**：不做擬合，值 = 初值 l0（500 um），畫在 log 軸最左端、
%             標成 0（log 軸沒有 0，位置是視覺錨點、非座標值）。
%             ⚠ 內層環全在赤道面上 → 這族點集**共平面**，面外方向約束極弱，
%               小 N 的劇烈震盪正是由此而來（這就是本設計要呈現的東西）。
%
%   縱軸 QTY：'ell'（l_hat [um]，藍）| 'gain'（g_I_hat [mT/A]，紅）；
%   single(1-param) 實心圓、eighteen(18-param bias) 空心方。
%   基準（不畫在圖上、只印 console）= 球內全部原始 .fld 空氣格點（20um 格、1771 點）。
%
%   收斂判準（console 報告，不畫線）：**相對前一點**的變化率，連續後 5 點皆 < 0.5%。
%
%   風格①粗體框圖：FS 36 粗體、box on、grid off、刻度奇數等距、曲線首末點貼齊框邊、
%   起訖數字以 text 補、圖例照 figure-style 標準樣式。
%
%   輸出 → figures/paper_fig/Section2_E/<qty>_vs_<tag>_maxwell_R150.png
% =========================================================================
    clc;
    if nargin < 1 || isempty(QTY),    QTY    = 'ell';   end
    if nargin < 2 || isempty(DESIGN), DESIGN = 'ratio'; end
    if nargin < 3 || isempty(force),  force  = false;   end

    R     = 150e-6;                  % 取樣半徑 [m]
    l0    = 0.5e-3;                  % l_hat 初值 [m]（= N=0 錨點的值）
    TOL   = 0.005;                   % 收斂判準：相對前一點變化 < 0.5%
    KWIN  = 5;                       % 連續幾點都要滿足
    MODEL = 'long2016_hexapole_halfcut';   GEOM = 'tip40um';

    specs = {};
    switch lower(DESIGN)
        case 'ratio'
            for k = 1:5, specs{end+1} = struct('NRPT',[k, 3*k, round(3*pi*k)]); end %#ok<AGROW>
            tag = 'grid';    XSCALE = 'log';
        case 'radial'
            for k = 1:12, specs{end+1} = struct('NRPT',[k, 15, 47]); end %#ok<AGROW>
            tag = 'radial';  XSCALE = 'linear';
        case 'radial150'
            % [ADDED 2026-08-13] N_r 由 1 掃到 150（對數階梯 13 檔；1..150 全跑的話有 140 個
            %   設計彼此重複資訊）。角向在小 N_r 用「格子方正」比例 (3N_r, 3pi*N_r)，到
            %   N_r=5 剛好長成 **15x47** 後就固定 —— 故最後幾檔都是 15x47（使用者指定）。
            %   ⚠ 前 5 檔是「徑向+角向一起長」，N_r>=5 之後才是純徑向，圖說要分界。
            %   最內殼 r1 = 150*(0.5/N_r)^(1/3)：119.1 um (N_r=1) → 22.4 um (N_r=150)，
            %   會跨過單殼可用半徑門檻（l_hat 30 um ≈ N_r 62、gain-eighteen 65 um ≈ N_r 6）。
            % 低點數段：明確指定三元組（維持 N_phi : N_theta ~ 1 : pi 的球面方正），
            %   讓 N 落在 27 / 52 / 114 / 200 / 400 這種尺度 —— 純靠比例式 s 驅動會因
            %   N_r 取整而在 56 → 140 → 336 → 576 之間跳空，補不出這些點。
            % 總共 **150 個設計**，刻意把大部分配在低點數端。
            % 低段（84 檔）：N_r = 1..6 × N_phi = 2..15（N_theta = round(pi*N_phi)）
            %   → N = 12 … 4230，光是 N < 1000 就有 40 檔以上。
            %   ⚠ N_phi 從 2 起：N_phi=1 只有赤道一條帶 ⇒ 點全共平面、g_I 會翻負（實測過）。
            for k = 1:6
                for m = 2:15
                    specs{end+1} = struct('NRPT',[k, m, round(pi*m)]); %#ok<AGROW>
                end
            end
            % 高段（66 檔）：角向固定 15x47，N_r 由 7 取到 150 → N = 4935 … 105750
            for k = unique(round(linspace(7,150,66)))
                specs{end+1} = struct('NRPT',[k, 15, 47]); %#ok<AGROW>
            end
            tag = 'radial150'; XSCALE = 'log';
        case 'tri200'
            % [ADDED 2026-08-13] 三個參數同時遞增的單一序列（使用者指定，圖上 200 點）。
            %   規則 = 「按比例輪流 +1」（Bresenham 配額）：起點 [1 2 3]，每步比較三個
            %   相對進度 N_r/1、N_phi/3、N_theta/(3pi)，把最落後的那個 +1（平手取第一個）。
            %   ⇒ 三數全程只增不減、始終貼著 1 : 3 : 3pi 的方正配比 ⇒ 一個 N 一個設計、無鋸齒。
            %   ⚠ 起點不可用 [1 2 2]：N_theta=2 的方位角是 90/270 度，兩者都在 y-z 平面
            %     ⇒ 全部點共平面、面外無約束，g_I 會翻負（實測 -612 mT/A）。
            specs{end+1} = struct('anchor',l0);          % N=0 錨點（值 = 初值 500 um）
            w = [1, 3, 3*pi];   t = [1 2 3];
            for q = 1:200
                specs{end+1} = struct('NRPT',t); %#ok<AGROW>
                [~, j] = min(t ./ w);   t(j) = t(j) + 1;
            end
            tag = 'tri200';  XSCALE = 'log';
        case 'angular'
            specs = {struct('NRPT',[1 2 1]), struct('NRPT',[1 1 3]), struct('NRPT',[1 2 2])};
            for m = [1 2 3 4 5 6 7 8 9 10 12 14 16 18 21 24 28 32]
                specs{end+1} = struct('NRPT',[2, m, round(pi*m)]); %#ok<AGROW>
            end
            tag = 'angular'; XSCALE = 'log';
        case 'inner'
            % 混合版（使用者拍板 2026-08-13）：內層自訂 + 外層等體積，**兩段都逐層累積**。
            %   內層：r = 5,10,…,65 um（13 層，止於 (5,15,47) 的最內殼 69.6 um 之前），
            %         第 k 層 n = 5k 點，**Fibonacci 球面撒點**（不是赤道環——環會讓整組
            %         點共平面、面外無約束，g_I 會翻成負值且爆到 -600 mT/A，已實測）。
            %   外層：(5,15,47) 的 5 顆等體積殼（r = 69.6/100.4/119.1/133.2/144.8 um、
            %         每顆 705 點）**一顆一顆加**，不要一次全上（否則 455 → 3980 一步跳完）。
            %   → 曲線 19 個點：0, 5, 15, 30, 50, 75, 105, 140, 180, 225, 275, 330, 390,
            %      455, 1160, 1865, 2570, 3275, 3980。
            tag = 'inner';   XSCALE = 'log';
        otherwise, error('DESIGN 必為 ''ratio'' | ''radial'' | ''angular'' | ''inner''');
    end

    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    cachef = fullfile(here, 'data', sprintf('ell_gain_%s_maxwell_R150.mat', tag));

    %% ---- 計算（快取）------------------------------------------------------
    if exist(cachef,'file') && ~force
        S = load(cachef);   fprintf('由快取載入 %s\n', cachef);
    elseif strcmpi(DESIGN,'inner')
        S = sweep_hybrid(MODEL, GEOM, R, l0);
        save(cachef, '-struct', 'S');   fprintf('已存 %s\n', cachef);
    else
        S = sweep_designs(MODEL, GEOM, R, specs, l0);
        save(cachef, '-struct', 'S');   fprintf('已存 %s\n', cachef);
    end

    switch lower(QTY)
        case 'ell'
            L1 = S.ell_single;   L2 = S.ell_bias;
            b1 = S.ell_ref_single;   b2 = S.ell_ref_bias;
            cq = [0.05 0.10 0.95];   unit = 'um';   fmt = '%.1f';
            ylab = '$\mathbf{\hat{\ell}\;(micro\;meter)}$';
        case 'gain'
            L1 = S.gI_single;    L2 = S.gI_bias;
            b1 = S.gI_ref_single;    b2 = S.gI_ref_bias;
            cq = [0.85 0.10 0.10];   unit = 'mT/A';  fmt = '%.4f';
            ylab = '$\mathbf{{}^{B}\hat{g}_{I}\;(mT/A)}$';
        otherwise, error('QTY 必為 ''ell'' | ''gain''');
    end
    N = S.N;

    % gain：N<=4 的值無意義（N<=2 電荷是精確解、值隨初值變；N=4 的 eighteen 用 4 點擬
    %   18 參數給出 350 mT/A，一點就撐掉整張圖）→ 不畫，值仍留在 .mat。
    %   但 **N=0 的槽位保留**（值為 NaN、不畫 marker），讓水平軸一樣從 0 起
    %   —— gain 沒有「初值」這種東西（電荷是被解出來的，不是猜的），所以 0 那點不畫值。
    if strcmpi(QTY,'gain')
        m6 = (N >= 6) | (N == 0);   N = N(m6);   L1 = L1(m6);   L2 = L2(m6);
    end
    ok = isfinite(L1) | isfinite(L2) | (N == 0);
    N = N(ok);   L1 = L1(ok);   L2 = L2(ok);
    % [ADDED] 依 N 由小到大排序。設計清單的產生順序不保證 N 遞增（radial150 的低段是
    %   N_r 外迴圈、N_phi 內迴圈 ⇒ N 會跳回頭），不排序的話折線會在圖上來回穿梭，
    %   而且「相對前一點」的收斂判準會變成沿著錯誤的順序在算。
    [N, iso] = sort(N);   L1 = L1(iso);   L2 = L2(iso);

    fprintf('\n  N        : %s\n', num2str(N, '%8d'));
    fprintf(['  single   : ' repmat([fmt ' '],1,numel(N)) ' %s（基準 ' fmt '）\n'], L1, unit, b1);
    fprintf(['  eighteen : ' repmat([fmt ' '],1,numel(N)) ' %s（基準 ' fmt '）\n'], L2, unit, b2);
    fprintf('  vs 基準 [%%]: single %s\n',   num2str(100*(L1/b1-1), '%+8.2f'));
    fprintf('              eighteen %s\n',  num2str(100*(L2/b2-1), '%+8.2f'));
    % [MODIFIED 2026-08-13] 收斂點 = **單一定義**：l_hat 穩定（後 KWIN 步相對前一步 < TOL）
    %   **且** K_I 符合物理（對角全正、off-diag 全負、對角占優），該點與其後 KWIN 步都合格。
    %   ⇒ ell 圖與 gain 圖畫**同一組虛線**，不再各自用自己的曲線算收斂點。
    %   ⚠ 只看 l_hat 會選到 K_I 壞掉的設計：R=150 的 single 在 N=8 就通過 l_hat 判準，
    %     但 30 個 off-diag 有 13 個是正的、對角線有一顆掉到 0.547（要到 N=80 才歸位）。
    EL1 = S.ell_single(:).';   EL2 = S.ell_bias(:).';
    KI1 = logical(S.ki_single(:).');   KI2 = logical(S.ki_bias(:).');
    [~, iso0] = sort(S.N);                                   % 與繪圖同一套排序
    Nall = S.N(iso0);  EL1 = EL1(iso0);  EL2 = EL2(iso0);  KI1 = KI1(iso0);  KI2 = KI2(iso0);
    m0 = Nall > 0;                                           % 錨點 N=0 不參與判準
    Nc1 = report_stable_ki('single  ', Nall(m0), EL1(m0), KI1(m0), TOL, KWIN);
    Nc2 = report_stable_ki('eighteen', Nall(m0), EL2(m0), KI2(m0), TOL, KWIN);

    %% ---- 畫圖 -------------------------------------------------------------
    FS = 36;   LW = 3.0;   MS = 10;
    fig = figure('Color','w','Position',[100 100 1120 820]);
    ax  = axes(fig);   hold(ax,'on');

    % log 軸畫不了 N=0 → 把錨點擺在最左端（視覺錨點，端點數字另標為 0）
    Np = N;   zero_anchor = strcmpi(XSCALE,'log') && any(N == 0);
    if zero_anchor
        % 錨點往左擺**一個十年**（不是半格）：log 軸的刻度是整十次冪，軸範圍太窄時
        % 扣掉兩端淨空後會只剩一個刻度（踩過：錨點放 min/2 → 範圍 25~7775 只剩 10^2）。
        pos = N(N > 0);   Np(N == 0) = min(pos)/10;
    end

    % [MODIFIED 2026-08-13 使用者拍板] 視覺編碼改回**顏色 = 模型**：single 藍、eighteen 紅
    %   （原本是「顏色 = 物理量、marker = 模型」，兩條同色只靠圓/方分辨，太不好認）。
    %   收斂虛線也跟著各自的模型上色，見下方。
    c1 = [0.05 0.10 0.95];   c2 = [0.85 0.10 0.10];
    h1 = plot(ax, Np, L1, '-o', 'Color',c1, 'LineWidth',LW, 'MarkerSize',MS, ...
              'MarkerFaceColor',c1, 'Clipping','off');
    h2 = plot(ax, Np, L2, '-s', 'Color',c2, 'LineWidth',LW, 'MarkerSize',MS, ...
              'MarkerFaceColor','w', 'Clipping','off');

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.02 .02]);
    ax.Toolbar.Visible = 'off';

    XL = [min(Np) max(Np)];   xlim(ax, XL);   set(ax,'XScale',XSCALE);
    if strcmpi(XSCALE,'log'), XT = log_ticks(XL(1), XL(2));
    else,                     XT = inner_ticks(XL(1), XL(2), 3);   end
    set(ax,'XTick',XT);

    [YL, YT] = axlim_auto(min([L1 L2]), max([L1 L2]), [3 5]);
    ylim(ax, YL);   set(ax,'YTick',YT);

    % 收斂點（後 KWIN 點相對前一點變化皆 < TOL）以虛線標出，**顏色跟著各自的模型**；
    % 兩模型同值時只畫一條、用深灰表示「兩者共用」。
    if ~isnan(Nc1) && ~isnan(Nc2) && Nc1 == Nc2
        xline(ax, Nc1, '--', 'Color',[0.25 0.25 0.25], 'LineWidth',2.5, 'HandleVisibility','off');
    else
        if ~isnan(Nc1), xline(ax, Nc1, '--', 'Color',c1, 'LineWidth',2.5, 'HandleVisibility','off'); end
        if ~isnan(Nc2), xline(ax, Nc2, '--', 'Color',c2, 'LineWidth',2.5, 'HandleVisibility','off'); end
    end

    % 水平軸起訖：只標數字、不畫 tick mark（起點若是錨點則標 0）
    yoff = YL(1) - 0.022*diff(YL);
    lbl  = {sprintf('%g',XL(1)), sprintf('%g',XL(2))};
    if zero_anchor, lbl{1} = '0'; end
    for q = 1:2
        text(ax, XL(q), yoff, lbl{q}, 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end

    xlabel(ax, '$\mathbf{Number\;of\;points}$', 'Interpreter','latex', 'FontSize',FS);
    ylabel(ax, ylab, 'Interpreter','latex', 'FontSize',FS);

    lg = legend(ax, [h1 h2], {'Single parameter', 'Eighteen parameters'}, ...
                'Interpreter','tex', 'Location','northoutside', 'NumColumns',2);
    lg.FontSize = 24;  lg.FontWeight = 'bold';  lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;
    hold(ax,'off');

    drawnow;
    axp = get(ax,'Position');   lgp = get(lg,'Position');
    lgw = lgp(3);   lgh = lgp(4);
    GAPN = 0.022;   newTop = 1 - lgh - GAPN - 0.006;
    axp(4) = newTop - axp(2);   set(ax,'Position',axp);
    if lgw < 0.70*axp(3)
        set(lg, 'Position', [axp(1) + (axp(3)-lgw)/2, newTop + GAPN, lgw, lgh]);
    else
        set(lg, 'Position', [axp(1), newTop + GAPN, axp(3), lgh]);
    end

    out = fullfile(figdir, sprintf('%s_vs_%s_maxwell_R150.png', lower(QTY), tag));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function S = sweep_designs(MODEL, GEOM, R, specs, l0)
% 逐個設計取點 → 擬合；另算「全原始格點」基準。
%   spec 三種：.anchor（不擬合，值 = 初值）｜.NRPT｜.ring_r/.ring_n（自訂環，可與 NRPT 併用）
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    addpath(fullfile(CAL,'function'), fullfile(CAL,'common_path'), fullfile(CAL,'utils'));
    addpath(fullfile(CAL,'utils','long2016_hexapole_halfcut'));

    cfg = model_config(MODEL, GEOM);
    assert(all(cfg.s_source == 1), '本腳本假設 raw 已 all-source（Maxwell）');
    F = zeros(6, cfg.N_I);
    for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end

    raw = extract_maxwell_data(cfg, 'all', cfg.default_variant);
    ad  = build_actuator_data(raw, cfg);
    [P0, B0, n0] = cfg.select_ball(ad, R);
    fprintf('\n[基準] R=%.0f um 內原始格點 %d\n', R*1e6, n0);
    [S.ell_ref_single, S.gI_ref_single] = fit_one(P0, B0, cfg.Pc_base, l0, false, F);
    [S.ell_ref_bias,   S.gI_ref_bias  ] = fit_one(P0, B0, cfg.Pc_base, l0, true,  F);
    S.N_ref = n0;

    nD = numel(specs);
    [S.N, S.ell_single, S.ell_bias, S.gI_single, S.gI_bias] = deal(zeros(1,nD));
    [S.ki_single, S.ki_bias] = deal(false(1,nD));
    for d = 1:nD
        sp = specs{d};
        if isfield(sp,'anchor')                      % N=0：殘差無定義，值 = 初值
            S.N(d) = 0;
            S.ell_single(d) = sp.anchor*1e6;   S.ell_bias(d) = sp.anchor*1e6;
            S.gI_single(d)  = NaN;             S.gI_bias(d)  = NaN;
            S.ki_single(d)  = false;           S.ki_bias(d)  = false;
            fprintf('  N=     0 （初值錨點）：l_hat = %.1f um\n', sp.anchor*1e6);
            continue
        end
        P = [];   B = [];
        if isfield(sp,'ring_r')                      % 自訂環（measure frame，赤道面）
            Q = ring_query(sp.ring_r, sp.ring_n);
            [x,y,z,Bq] = sphere_grid_sample([], [], struct('frame','actuator','query',Q));
            P = [x y z];   B = Bq;
        end
        if isfield(sp,'NRPT')                        % 等測度網格（可與環併用）
            [x,y,z,Bq] = sphere_grid_sample(R, [], ...
                struct('frame','actuator','NRPT',sp.NRPT,'N_nodes',n0));
            P = [P; x y z];   B = cat(1, B, Bq);
        end
        Bs = zeros(3*size(P,1), size(B,3));
        for j = 1:size(B,3), Bs(:,j) = reshape(B(:,:,j).', [], 1); end
        S.N(d) = size(P,1);
        [S.ell_single(d), S.gI_single(d), S.ki_single(d)] = fit_one(P, Bs, cfg.Pc_base, l0, false, F);
        [S.ell_bias(d),   S.gI_bias(d),   S.ki_bias(d)  ] = fit_one(P, Bs, cfg.Pc_base, l0, true,  F);
        fprintf('  N=%6d：single l=%8.1f g=%9.4f｜eighteen l=%8.1f g=%9.4f\n', ...
                S.N(d), S.ell_single(d), S.gI_single(d), S.ell_bias(d), S.gI_bias(d));
    end
    S.R = R;   S.specs = specs;   S.model = MODEL;   S.geom = GEOM;
end

% ============================================================================
function S = sweep_hybrid(MODEL, GEOM, R, l0)
% 混合設計：內層自訂 Fibonacci 球殼（逐層累積）+ 外層 (5,15,47) 等體積殼（逐殼累積）。
%   只呼叫取樣器兩次（內層一次 query、外層一次 NRPT），之後靠半徑分層取子集，
%   所以 19 個設計不必重複讀場、也不必在此複製產點公式。
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    addpath(fullfile(CAL,'function'), fullfile(CAL,'common_path'), fullfile(CAL,'utils'));
    addpath(fullfile(CAL,'utils','long2016_hexapole_halfcut'));

    cfg = model_config(MODEL, GEOM);
    F = zeros(6, cfg.N_I);
    for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end

    % ---- 基準（全原始格點）----
    raw = extract_maxwell_data(cfg, 'all', cfg.default_variant);
    ad  = build_actuator_data(raw, cfg);
    [P0, B0, n0] = cfg.select_ball(ad, R);
    fprintf('\n[基準] R=%.0f um 內原始格點 %d\n', R*1e6, n0);
    [S.ell_ref_single, S.gI_ref_single] = fit_one(P0, B0, cfg.Pc_base, l0, false, F);
    [S.ell_ref_bias,   S.gI_ref_bias  ] = fit_one(P0, B0, cfg.Pc_base, l0, true,  F);
    S.N_ref = n0;

    % ---- 內層：10 顆 Fibonacci 球殼 ----
    %   [MODIFIED 2026-08-13] 半徑照原規格停在 50 um（外層最內殼 55.3 um，剛好錯開不重疊）；
    %   每層點數改從 50 起、步長 5（n = 50,55,…,95）。原本 n = 5k 時第一個設計只有 5 點、
    %   eighteen 的 l_hat 會噴到 1147.8 um（近簡併）；50 點壓到 701.6，尖刺消失。
    %   ⚠ r=5 um 的值本身仍是**偏的**（654/702 vs 基準 874）——那顆球整個在一個 20 um
    %     網格單元內，場是三線性內插的產物、近乎線性。加點只會更穩定地停在錯值
    %     （500 點 → 658.9）。對照：同樣 500 點放 r=50 um 給 874.0、r=65 um 給 877.1。
    rin = (5:5:50)*1e-6;   nin = 50:5:95;
    Q = zeros(sum(nin),3);   p = 0;
    for k = 1:numel(rin)
        Q(p+(1:nin(k)),:) = rin(k) * fib_sphere(nin(k));   p = p + nin(k);
    end
    [x,y,z,Bi] = sphere_grid_sample([], [], struct('frame','actuator','query',Q));
    Pi = [x y z];   ki = round(vecnorm(Pi,2,2) / 5e-6);        % 每點屬於第幾層（半徑判定）

    % ---- 外層：(5,15,47) 等體積 5 殼 ----
    [x,y,z,Bo] = sphere_grid_sample(R, [], ...
        struct('frame','actuator','NRPT',[10 15 47],'N_nodes',n0));   % 10 顆殼（角向不動）
    Po = [x y z];   ro = vecnorm(Po,2,2);
    rshell = unique(round(ro*1e9))/1e9;                        % 5 個殼半徑
    [~, ko] = min(abs(ro - rshell(:).'), [], 2);
    fprintf('[外層] 殼半徑 [um] = %s（每殼 %d 點）\n', ...
            num2str(rshell(:).'*1e6,'%.1f '), nnz(ko==1));

    % ---- 19 個設計：錨點 + 內層累積 13 + 外層累積 5 ----
    nD = 1 + numel(rin) + numel(rshell);
    [S.N, S.ell_single, S.ell_bias, S.gI_single, S.gI_bias] = deal(zeros(1,nD));
    [S.ki_single, S.ki_bias] = deal(false(1,nD));
    S.N(1) = 0;   S.ell_single(1) = l0*1e6;   S.ell_bias(1) = l0*1e6;
    S.gI_single(1) = NaN;   S.gI_bias(1) = NaN;   S.ki_single(1) = false;   S.ki_bias(1) = false;
    fprintf('  N=     0 （初值錨點）：l_hat = %.1f um\n', l0*1e6);

    for d = 2:nD
        j = d - 1;
        if j <= numel(rin)                                     % 內層累積
            m = ki <= j;   P = Pi(m,:);   B = Bi(m,:,:);
        else                                                   % 內層全上 + 外層累積
            g = j - numel(rin);
            P = [Pi; Po(ko<=g,:)];   B = cat(1, Bi, Bo(ko<=g,:,:));
        end
        Bs = zeros(3*size(P,1), size(B,3));
        for q = 1:size(B,3), Bs(:,q) = reshape(B(:,:,q).', [], 1); end
        S.N(d) = size(P,1);
        [S.ell_single(d), S.gI_single(d), S.ki_single(d)] = fit_one(P, Bs, cfg.Pc_base, l0, false, F);
        [S.ell_bias(d),   S.gI_bias(d),   S.ki_bias(d)  ] = fit_one(P, Bs, cfg.Pc_base, l0, true,  F);
        fprintf('  N=%6d：single l=%8.1f g=%9.4f｜eighteen l=%8.1f g=%9.4f\n', ...
                S.N(d), S.ell_single(d), S.gI_single(d), S.ell_bias(d), S.gI_bias(d));
    end
    S.R = R;   S.rin = rin;   S.nin = nin;   S.rshell = rshell;
    S.model = MODEL;   S.geom = GEOM;
end

% ============================================================================
function U = fib_sphere(n)
% n 個準均勻分布的**單位球面**方向（Fibonacci / 黃金角螺旋）。任意 n 皆可，無共平面問題。
    i  = (1:n).';
    z  = 1 - (2*i - 1)/n;                       % 等分 z（等面積帶的格心）
    th = 2*pi*i / ((1+sqrt(5))/2);              % 黃金角遞增方位
    s  = sqrt(max(1 - z.^2, 0));
    U  = [s.*cos(th), s.*sin(th), z];
end

% ============================================================================
function Q = ring_query(r, n)
% 赤道面上的同心環（measure frame、球心為原點）：第 k 環半徑 r(k)、n(k) 個點。
%   每環自己的 theta 偏移 (i-0.5)/n(k) 不同 → 各環不會徑向對齊成輻條。
    Q = zeros(sum(n), 3);   p = 0;
    for k = 1:numel(r)
        th = 2*pi*((1:n(k)) - 0.5)/n(k);
        Q(p+(1:n(k)), :) = [r(k)*cos(th(:)), r(k)*sin(th(:)), zeros(n(k),1)];
        p = p + n(k);
    end
end

% ============================================================================
function [ell_um, gI, kiok] = fit_one(P, Bstack, Pc_base, l0, USE_BIAS, F)
% 極小/退化設計可能讓 S 秩虧 → NaN。吞例外、回 NaN，掃描不中斷（圖上該點缺席）。
%   [ADDED 2026-08-13] 一併回傳 K_I 是否符合物理（對角全正、off-diag 全負、對角占優）。
    try
        [e, l_hat] = fitting(P, Bstack, Pc_base, l0, USE_BIAS);
        [K, gI]    = solve_current(l_hat, e, Pc_base, P, Bstack, F);
        ell_um     = l_hat * 1e6;   kiok = ki_ok(K);
    catch ME
        warning('fit_one:degenerate','N=%d 擬合失敗（%s）→ NaN', size(P,1), ME.identifier);
        ell_um = NaN;   gI = NaN;   kiok = false;
    end
end

% ============================================================================
function tf = ki_ok(K)
% K_I 符合物理：對角全正、off-diagonal 全負、對角占優（每列最大絕對值在對角）。
%   列和 ~ 0（電荷中性）是 emergent、不列硬條件（全格點版也只到 0.017~0.023）。
    if isempty(K) || any(~isfinite(K(:))), tf = false;  return; end
    M = ~eye(6);   [~, im] = max(abs(K), [], 2);
    tf = all(diag(K) > 0) && all(K(M) < 0) && isequal(im(:).', 1:6);
end

% ============================================================================
function Nc = report_stable_ki(name, N, v, ok, tol, K)
% 收斂點（雙條件）：第一個 i，使 ① 其後連續 K 步「相對前一點」變化率 < tol；
%                            ② 該點與其後 K 步的 K_I 都合格。未達判準回 NaN。
    rel = abs(diff(v) ./ v(1:end-1));
    Nc  = NaN;
    for i = 1:numel(rel)-K+1
        if all(rel(i:i+K-1) < tol) && all(ok(i:min(i+K, numel(ok))))
            Nc = N(i);
            fprintf('  收斂點（%s）：N = %d 起（l_hat 後 %d 步 < %.1f%% 且 K_I 合格）\n', ...
                    name, Nc, K, tol*100);
            return
        end
    end
    fprintf('  收斂點（%s）：未達雙條件判準\n', name);
end

% ============================================================================
function Nc = report_stable(name, N, v, tol, K)
% 收斂點：第一個 i，使其後連續 K 個「相對前一點」的變化率都 < tol。回傳該 N（未達判準回 NaN）。
    rel = abs(diff(v) ./ v(1:end-1));
    i0  = [];
    for i = 1:numel(rel)-K+1
        if all(rel(i:i+K-1) < tol), i0 = i;  break; end
    end
    if isempty(i0)
        Nc = NaN;
        fprintf('  收斂點（%s）：未達判準（後 %d 點相對變化皆 < %.1f%%）\n', name, K, tol*100);
    else
        Nc = N(i0);
        fprintf('  收斂點（%s）：N = %d 起，其後 %d 點相對變化皆 < %.1f%%\n', ...
                name, Nc, K, tol*100);
    end
end

% ============================================================================
function XT = log_ticks(lo, hi)
% log 軸內部刻度：等比、奇數個、與兩端留足淨空（端點另用 text 標數字，會互相疊到）。
%   [MODIFIED] 淨空改用「軸總跨距的 15%」而非固定一倍：固定倍數在長軸上不夠
%   （tri200 跨 5.2 個十年，20000 離右端 5.3 倍，標籤仍與端點數字 105248 疊在一起），
%   在短軸上又太嚴（inner 只跨 3.2 個十年，一刀切會只剩一個刻度）。
    gap  = 10^(0.15 * log10(hi/lo));                  % 兩端各留的淨空倍率
    best = [];
    for mult = [1 2 5]
        e = floor(log10(lo)) : ceil(log10(hi));
        t = mult * 10.^e;
        t = t(t >= gap*lo & t <= hi/gap);
        % [MODIFIED] log 軸**不套「刻度數量取奇數」**。那條規則是給 2D 線性軸求對稱用的
        %   （讓 0 落在正中），對數軸的刻度是整十次冪、位置由資料範圍決定，本來就不對稱，
        %   硬套只會平白砍掉一個刻度（tri200 的 10^4 離右端有 10.5 倍距離、根本不會撞標籤，
        %   卻因為湊出 4 個是偶數而被丟掉）。淨空檢查保留。
        if numel(t) > numel(best), best = t; end
    end
    XT = best;
end

% ============================================================================
function tk = inner_ticks(lo, hi, n)
% 線性軸的 n 個等距內部刻度（不含首末端點；n 取奇數）。
    cand = [1 2 2.5 3 4 5 10];
    s = (hi - lo)/(n+1);
    k = floor(log10(s));
    [~, i] = min(abs(cand*10^k - s));
    s = cand(i)*10^k;
    ctr = round(((lo+hi)/2)/s)*s;
    tk  = ctr + (-(n-1)/2 : (n-1)/2)*s;
    tk  = tk(tk > lo & tk < hi);
end

% ============================================================================
function [lim, tk] = axlim_auto(lo, hi, nlist)
% 奇數個等距 tick、兩端留白 = tick 間距（house style：各腳本自帶一份）。
    cand = [1 2 2.5 3 4 5 10];
    mid  = (lo+hi)/2;   rng_ = max(hi-lo, realmin);
    best = {};   bestSpan = inf;
    for n = nlist
        k0 = floor(log10(rng_/(n+1)));
        for k = k0:(k0+4)
            hit = false;
            for c = cand
                s   = c*10^k;
                ctr = round(mid/s)*s;
                T = { ctr + (-(n-1)/2 : (n-1)/2)*s };
                if lo >= 0 && lo < s, T = [{(0:n-1)*s}, T]; end %#ok<AGROW>
                for it = 1:numel(T)
                    t   = T{it};
                    L_  = [t(1)-s, t(end)+s];
                    clr = 0.15*s;
                    if lo >= L_(1)+clr && hi <= L_(2)-clr
                        if (n+1)*s < bestSpan, bestSpan = (n+1)*s;  best = {L_, t}; end
                        hit = true;  break;
                    end
                end
                if hit, break; end
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
