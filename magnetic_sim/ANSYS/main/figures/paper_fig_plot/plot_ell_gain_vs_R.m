function plot_ell_gain_vs_R(USE_BIAS, SRC, SAMPLING)
% plot_ell_gain_vs_R -- long2016 半切六極 paper 圖:ℓ̂ 與 ĝ_I 隨取樣半徑 R 的趨勢
% =========================================================================
%   graded FEM 資料 + current_base 校正(fix 單一 ℓ,USE_BIAS=false)。
%   掃描取樣半徑 R = 40~500 µm,每點:select_ball → fitting(fix) → solve_current。
%   上 panel = ℓ̂(R) [µm]、下 panel = ĝ_I(R) [mT/A](尺度差大 → 分兩 panel)。
%   前段(讀 .dat + build_actuator_data)只做一次;loop R 只重跑 select_ball/fitting/solve_current。
%   風格 = 選項①粗體框 @ paper scale(font 粗體、LineWidth3、box on、grid off、tick 減半、單位 ())。
%   輸出 → figures/paper_fig/Section2_E/ell_gain_vs_R.png(覆蓋迭代)。
%
%   ⚠ SRC='maxwell' 預設 SAMPLING='interp'：場先內插到 10µm 細格再掃描。
%     **此內插純為圖面可讀性**（消除硬球切過 20µm 晶格造成的鋸齒，二階差分 RMS 降約 3×），
%     **不改變任何結論**——ℓ̂ 與趨勢形狀不變、R=140 的 ĝ_I 僅差 0.7%(fix)/0.14%(bias)。
%     校正主結果（data/.mat + results/*.pdf）一律仍用**原始格點**，不吃內插偏差。
%     論文 caption 須揭露場經內插取樣（專案規則：內插必標示）。
% =========================================================================
    clc;
    if nargin < 1, USE_BIAS = false; end   % false = fix → ell_gain_vs_R.png；true = 18-param bias → ell_gain_vs_R_bias.png
    if nargin < 2, SRC = 'apdl';    end    % [ADDED] 求解器來源:'apdl'(ANSYS .dat) | 'maxwell'(.fld) → 檔名加 _maxwell
    % [ADDED] 取樣方式:'raw'   = 直接用原始節點/格點(Maxwell 原始格距 20µm)
    %                  'interp'= 先把場**內插到 10µm 細格**(見 local interp_grid_sample)再照舊掃描
    %                            → 球面切過晶格的顆粒感變細,ĝ_I(R) 抖動按格距等比下降。maxwell 預設。
    %   ⚠ 內插是「重新取樣」,不增加資訊;論文圖說須標明場經內插取樣。
    if nargin < 3 || isempty(SAMPLING)
        if strcmpi(SRC,'maxwell'), SAMPLING = 'interp'; else, SAMPLING = 'raw'; end
    end
    here   = fileparts(mfilename('fullpath'));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    CAL = solver_path(SRC);                                          % [ADDED] 依 SRC 切分支(內含 rmpath 另一分支防遮蔽)

    bstr = ''; if USE_BIAS, bstr = '_bias'; end
    sstr = ''; if strcmpi(SRC,'maxwell'), sstr = '_maxwell'; end     % [ADDED] 來源後綴
    % [ADDED] 掃描結果快取:R-sweep 要 24 次擬合(每次還得先讀 656k 節點 .dat)、數分鐘;
    %   純改樣式(軸標題/tick/顏色)不需重算 → 有快取就載入。要強制重算把此 .mat 刪掉。
    istr = ''; if strcmpi(SAMPLING,'interp'), istr = '_interp'; end    % [ADDED] 取樣方式後綴(只進快取名,圖檔名不變→覆蓋)
    cachef = fullfile(here, sprintf('ell_gain_sweep%s%s%s.mat', sstr, istr, bstr));

    if exist(cachef, 'file')
        S = load(cachef);   Rum = S.Rum;  ell_R = S.ell_R;  gI_R = S.gI_R;
        fprintf('loaded cache %s  (delete it to force recompute)\n', cachef);
    else
        % ---- 前段 pipeline(只做一次):graded / current / actuator frame ----
        cfg = model_config('long2016_hexapole_halfcut','tip40um');
        raw = load_raw(SRC, cfg);                           % [MODIFIED] 依 SRC:apdl='graded' .dat / maxwell=.fld
        ad  = build_actuator_data(raw, cfg);   Pc_base = ad.Pc_base;
        F   = zeros(6, cfg.N_I);
        for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end
        if strcmpi(SAMPLING,'interp')                       % [ADDED] 先把場內插到更細的均勻格(10µm),再照舊掃描
            ad = interp_grid_sample(ad, cfg, 500e-6, 10e-6);
        end

        % ---- 掃描 R(fix 模型;每 20µm 一點;從 ĝ_I 正值起 R=40µm;R≤20 病態非物理)----
        Rum = 40:20:500;                                    % 取樣半徑 [µm](固定 20µm 間距)
        nR  = numel(Rum);
        ell_R = nan(1,nR);  gI_R = nan(1,nR);  npts_R = zeros(1,nR);  J_R = nan(1,nR);
        fprintf('\n  R[um]   npts    ell_hat[um]   gI_hat[mT/A]      J[mT^2]\n');
        for i = 1:nR
            [P, Bstack, npts] = cfg.select_ball(ad, Rum(i)*1e-6);
            npts_R(i) = npts;
            if npts < 3, fprintf('  %5d  %6d   (skip: npts<3)\n', Rum(i), npts); continue; end
            [e, l_hat, J] = fitting(P, Bstack, Pc_base, 0.5e-3, USE_BIAS);   % fix / bias 由旗標切
            [~, gI_hat]   = solve_current(l_hat, e, Pc_base, P, Bstack, F);
            ell_R(i) = l_hat*1e6;  gI_R(i) = gI_hat;  J_R(i) = J;   % [ADDED] J = 擬合殘差平方和 [mT²]
            fprintf('  %5d  %6d   %10.2f   %11.4f   %12.4g\n', Rum(i), npts, ell_R(i), gI_R(i), J_R(i));
        end
        save(cachef, 'Rum', 'ell_R', 'gI_R', 'npts_R', 'J_R', 'USE_BIAS');   % [MODIFIED] 加存 J_R
        fprintf('saved cache %s\n', cachef);
    end

    % ---- 繪圖範圍(SRC 相依) --------------------------------------------------
    % [ADDED] maxwell 是 0.02mm **均勻格**:R=40/60µm 球內只有 28/106 個格點 → 擬合病態
    %   (ĝ_I 衝到 27.6/11.3,非物理)。故 maxwell 從 **R=80µm**(252 點)起畫;快取仍存全範圍。
    %   apdl 是 graded 網格,近場極密,40µm 起即穩定 → 維持 40。
    if strcmpi(SRC,'maxwell'), Rmin = 80; else, Rmin = 40; end
    m = Rum >= Rmin;   Rum = Rum(m);  ell_R = ell_R(m);  gI_R = gI_R(m);
    XT = unique([Rmin 100 200 300 400 500]);            % 起點 + 百位 + 終點(照 APDL 樣式)
    XT(XT > Rmin & XT < Rmin+40) = [];                  % [ADDED] 起點與鄰近百位太近會壓字(80 vs 100)→ 去掉

    % ---- 畫圖(2×1 panel,選項①粗體框)----
    cL = [0.05 0.10 0.95];   cR = [0.85 0.10 0.10];     % ℓ̂ 亮藍 / ĝ_I 紅
    FS = 30;                                            % 對齊 circuit 圖字體
    fig = figure('Color','w','Position',[100 80 980 940]);
    t = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

    % 上:ℓ̂(R)
    ax1 = nexttile(t);  hold(ax1,'on');
    plot(ax1, Rum, ell_R, '-o', 'Color',cL, 'LineWidth',3.0, 'MarkerSize',8, 'MarkerFaceColor',cL);
    style_panel(ax1, FS);
    xlim(ax1,[Rmin 500]);   set(ax1,'XTick',XT);   % 含起點 + 終點 500
    [yl,yt] = ylim_pick(SRC, USE_BIAS, 'ell', ell_R);   % [MODIFIED] apdl 沿用定案值;maxwell 自動等距 4 內縮 tick
    ylim(ax1,yl);  set(ax1,'YTick',yt);
    ax1.XTickLabel = {};                                % 上 panel x 數字隱藏(共用軸)
    ylabel(ax1, '$\mathbf{\hat{\ell}\;(\mu m)}$', 'Interpreter','latex', 'FontSize',36);   % [ADDED] 上 panel y 軸標題（標準數學字體 \mathbf CM，同下 panel）
    hold(ax1,'off');

    % 下:ĝ_I(R)
    ax2 = nexttile(t);  hold(ax2,'on');
    plot(ax2, Rum, gI_R, '-s', 'Color',cR, 'LineWidth',3.0, 'MarkerSize',8, 'MarkerFaceColor',cR);
    style_panel(ax2, FS);
    xlim(ax2,[Rmin 500]);   set(ax2,'XTick',XT);   % 含起點 + 終點 500
    [yl,yt] = ylim_pick(SRC, USE_BIAS, 'g', gI_R);      % [MODIFIED] 同上
    ylim(ax2,yl);  set(ax2,'YTick',yt);
    ylabel(ax2, '$\mathbf{{}^{B}\hat{g}_{I}\;(mT/A)}$', 'Interpreter','latex', 'FontSize',36);   % 下 panel y 軸標題（標準數學字體 \mathbf CM）
    xlabel(ax2, '$\mathbf{Sampling\;range\;(\mu m)}$', 'Interpreter','latex', 'FontSize',36);   % [ADDED] 共用 x 軸標題（掃描半徑 R）
    hold(ax2,'off');

    out = fullfile(figdir, sprintf('ell_gain_vs_R%s%s.png', bstr, sstr));   % bstr/sstr 已在頂部定義
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('\nwrote %s\n', out);
end

% ============================================================================
function adq = interp_grid_sample(ad, cfg, Rmax, h)
% [ADDED] 把場內插到「更細的均勻格」，取代原始格點取樣（其餘流程完全不變：一樣硬球選點）。
%   Maxwell 原始格距 20µm；內插到 h（預設 10µm）後，球面切過晶格的顆粒感變細
%   → ĝ_I(R) 抖動約按 h 等比下降。
%   ⚠ 內插不增加資訊，是重新取樣；論文圖說須標明場經內插取樣。
%
%   步驟：① 在 actuator frame 建 h 間距立方格、只留 |p| ≤ Rmax
%         ② 濾鐵（轉回 measure frame，套 filter_iron_nodes，與 build_actuator_data 同判準）
%         ③ 用 ad 的空氣點雲建 scatteredInterpolant（linear/none 不外插），
%            **三角化只建一次**、之後只換 .Values（6 激發 × 3 分量）

    % ① 細格 → 球內
    g  = -Rmax:h:Rmax;
    [X, Y, Z] = ndgrid(g, g, g);
    Pq = [X(:), Y(:), Z(:)];
    rq = vecnorm(Pq, 2, 2);
    in = rq <= Rmax;
    Pq = Pq(in,:);   rq = rq(in);

    % ② 濾鐵
    Pm  = (ad.R_act.' * Pq.').';                         % actuator → measure(WP frame)
    air = filter_iron_nodes(Pm(:,1), Pm(:,2), Pm(:,3) + cfg.SPH_OFST, cfg);
    Pq = Pq(air,:);   rq = rq(air);

    % ③ 內插（三角化只建一次）
    src = ad.r2 < (Rmax + 60e-6)^2;                      % 源點只取球外一點點（省 Delaunay）
    Xs  = ad.Pa(src,:);
    N_I = size(ad.Ba,3);
    Bq  = zeros(size(Pq,1), 3, N_I);
    Fi  = scatteredInterpolant(Xs(:,1), Xs(:,2), Xs(:,3), ad.Ba(src,1,1), 'linear', 'none');
    for j = 1:N_I
        for c = 1:3
            Fi.Values = ad.Ba(src,c,j);                  % 重用三角化，只換值
            Bq(:,c,j) = Fi(Pq);
        end
    end

    ok = ~any(any(isnan(Bq),3),2);                       % 落在源網格外者剔除
    fprintf('  [interp] 格距 %.0f um：球內 %d 點 → 濾鐵後 %d → 有效 %d（源點 %d）\n', ...
            h*1e6, nnz(in), size(Pq,1), nnz(ok), nnz(src));
    adq = ad;
    adq.Pa = Pq(ok,:);   adq.r2 = rq(ok).^2;   adq.Ba = Bq(ok,:,:);
end

% ============================================================================
function CAL = solver_path(SRC)
% [ADDED] 依 SRC 掛上對應分支、並「移除另一分支」——兩分支有同名函式
%   (model_config/fitting/solve_current/emit_*…)，只 addpath 不 rmpath 會被先進 path 的那份遮蔽。
    APDL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    MW   = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    switch lower(SRC)
        case 'apdl',    CAL = APDL;  OTHER = MW;
        case 'maxwell', CAL = MW;    OTHER = APDL;
        otherwise, error('SRC 必為 ''apdl'' | ''maxwell''');
    end
    warning('off','MATLAB:rmpath:DirNotFound');
    rmpath(fullfile(OTHER,'function'));  rmpath(fullfile(OTHER,'common_path'));
    warning('on','MATLAB:rmpath:DirNotFound');
    addpath(fullfile(CAL,'function'));   addpath(fullfile(CAL,'common_path'));
    fprintf('[SRC=%s] model_config -> %s\n', SRC, which('model_config'));   % 自證用哪一份 code
end

% ============================================================================
function raw = load_raw(SRC, cfg)
% [ADDED] 依 SRC 載入原始場：apdl = graded .dat / maxwell = .fld（cfg.default_variant='maxwell'）
    if strcmpi(SRC,'maxwell')
        raw = extract_maxwell_data(cfg, 'all', cfg.default_variant);
    else
        raw = extract_ansys_data(cfg, 'all', 'graded');
    end
end

% ============================================================================
function [yl, yt] = ylim_pick(SRC, USE_BIAS, which_q, v)
% [ADDED] y 軸範圍 / 刻度：apdl 沿用既定案值（逐字不變）；maxwell 依實際資料自動取
%   「等距 4 個內縮 tick」（照 figure-style：不含端點、等距、數量固定）。
    if strcmpi(SRC,'apdl')
        switch which_q
            case 'ell', yl = [750 870];  yt = [770 800 830 860];
            case 'g'
                if USE_BIAS, yl = [9.5 10.5];  yt = [9.7 9.9 10.1 10.3];
                else,        yl = [7.7 8.5];   yt = [7.8 8.0 8.2 8.4];   end
        end
        return;
    end
    lo = min(v(:));  hi = max(v(:));  sp = hi - lo;      % maxwell：自動
    yl = [lo - 0.15*sp, hi + 0.22*sp];                   % 上方留多一點:首點常是最高點,別貼到上框
    s  = nice_step(diff(yl)/5);                          % 等距步長（4 內縮 tick → 分 5 段）
    c  = round(mean(yl)/s)*s;
    yt = c + (-1.5:1:1.5)*s;                             % 4 個等距 tick，對稱於中心
    yt = yt(yt > yl(1) & yt < yl(2));                    % 保證不含端點
end

function s = nice_step(x)
% [ADDED] 取 nice 等距步長（1/2/2.5/3/4/5/10 × 10^k）
    k = floor(log10(x));   m = x/10^k;
    cand = [1 2 2.5 3 4 5 10];
    [~,i] = min(abs(cand - m));
    s = cand(i)*10^k;
end

% ============================================================================
function style_panel(ax, FS)
% 選項①粗體框 @ paper scale(tick 由呼叫端明設,不減半)。
    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.02 .02]);
    ax.Toolbar.Visible = 'off';
end
