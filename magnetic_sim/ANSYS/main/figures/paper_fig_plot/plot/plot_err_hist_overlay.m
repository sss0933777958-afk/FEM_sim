function plot_err_hist_overlay(SRC, Rsel_um, CONV)
% plot_err_hist_overlay -- long2016 半切六極 R=150µm：fix vs bias 絕對殘差「疊圖」
% =========================================================================
%   同一取樣球(R=150µm)、同一絕對殘差定義 |b_FEM − S·ĝ_I·K̄·I| (mT)：
%     藍 = 無 bias(fix 單一 ℓ)、紅 = 有 bias(18-param)。
%   共用 fix 的 x 刻度(xlim[0 1]、tick 0.2~0.8 + 起終點)；共用 edges(nb=180)。
%   percentage 縱軸、圖例標顏色對應、**不標平均線**、無軸標題。
%   輸出 → figures/paper_fig/Section2_E/err_hist_overlay.png(覆蓋迭代)。
% =========================================================================
    clc;
    if nargin < 1, SRC = 'apdl'; end       % [ADDED] 'apdl'(ANSYS .dat) | 'maxwell'(.fld) → 檔名加 _maxwell
    if nargin < 2 || isempty(Rsel_um), Rsel_um = 150; end   % [ADDED] 取樣球半徑 [µm]；≠150 時檔名加 _R<NNN>
    % [ADDED 2026-08-13] CONV=true：**改用該 R 的 N_c 減量設計校正**（判準＝l_hat 穩定
    %   ∧ g_I 穩定 ∧ K_I 符合物理，各持續 10 步），再拿該 R 內全部真實格點當評估集。
    %   兩個模型各自有自己的 N_c。檔名加 _conv。
    if nargin < 3 || isempty(CONV), CONV = false; end
    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    CAL = solver_path(SRC);                                          % [ADDED] 依 SRC 切分支(內含 rmpath 另一分支防遮蔽)
    sstr = ''; if strcmpi(SRC,'maxwell'), sstr = '_maxwell'; end     % [ADDED] 來源後綴

    % ---- 前段 pipeline(只做一次)----
    cfg = model_config('long2016_hexapole_halfcut','tip40um');
    raw = load_raw(SRC, cfg);                                        % [MODIFIED] apdl='graded' .dat / maxwell=.fld
    ad  = build_actuator_data(raw, cfg);   Pc_base = ad.Pc_base;
    F   = zeros(6, cfg.N_I);  for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end
    Rum = Rsel_um;                                          % [MODIFIED] 取樣半徑可調（原寫死 150）
    [P, Bstack, npts] = cfg.select_ball(ad, Rum*1e-6);
    fprintf('取樣範圍 R <= %d um：npts = %d\n', Rum, npts);

    % ---- 兩模型各擬合一次 → 逐節點×激發 絕對殘差 (mT) ----
    if CONV
        addpath(fullfile(CAL,'utils'), fullfile(CAL,'utils','long2016_hexapole_halfcut'));
        [Pa, Ba, na] = conv_set(Rum, cfg, F, false, here);
        [Pb, Bb, nb_] = conv_set(Rum, cfg, F, true,  here);
        fprintf('N_c：single %d 點｜eighteen %d 點（評估集 = %d 個真實格點）\n', na, nb_, npts);
        err0 = fit_abs_resid(P, Bstack, Pc_base, F, npts, false, Pa, Ba);
        err1 = fit_abs_resid(P, Bstack, Pc_base, F, npts, true,  Pb, Bb);
    else
        err0 = fit_abs_resid(P, Bstack, Pc_base, F, npts, false);   % single_param (無 bias)
        err1 = fit_abs_resid(P, Bstack, Pc_base, F, npts, true);    % eighteen_param (有 bias)
    end
    cv0 = std(err0)/mean(err0)*100;   cv1 = std(err1)/mean(err1)*100;   % CV = σ/μ [%]（不受 bin 寬影響）
    fprintf('single_param  : mean=%.4f mT  max=%.4f mT  CV=%.1f%%\n', mean(err0), max(err0), cv0);
    fprintf('eighteen_param: mean=%.4f mT  max=%.4f mT  CV=%.1f%%\n', mean(err1), max(err1), cv1);

    % ---- 各自用自己的 edges(不共用;各自然解析度、nb=180),只是畫在同一 x 軸 ----
    nb  = 180;
    edg0 = linspace(0, max(err0), nb+1);  ctr0 = (edg0(1:end-1)+edg0(2:end))/2;   % single 0~0.88
    edg1 = linspace(0, max(err1), nb+1);  ctr1 = (edg1(1:end-1)+edg1(2:end))/2;   % eighteen 0~0.21（自己的細 bin）
    pct0 = histcounts(err0, edg0) / numel(err0) * 100;
    pct1 = histcounts(err1, edg1) / numel(err1) * 100;

    [cB, cR] = pick_bar_colors(Rum);   % [MODIFIED] 長條配色依取樣半徑分組（同 plot_err_hist 的 pick_colors）
    FS = 28;
    fig = figure('Color','w','Position',[100 100 1100 830]);  ax = axes(fig);  hold(ax,'on');   % [MODIFIED] 加高補償外置兩列圖例
    h0 = bar(ax, ctr0, pct0, 1, 'FaceColor',cB, 'FaceAlpha',0.60, 'EdgeColor','k', 'LineWidth',0.3);   % 亮藍(single) + 黑細邊
    h1 = bar(ax, ctr1, pct1, 1, 'FaceColor',cR, 'FaceAlpha',0.60, 'EdgeColor','none');                 % 紅(eighteen) 無邊；帶透明看得到交錯
    % 紅色：沿輪廓畫連續黑線(silhouette)+ 挑幾根 bar 從上到下畫垂直黑邊
    stairs(ax, edg1, [pct1 pct1(end)], 'Color','k', 'LineWidth',0.5);   % 頂部輪廓(粗度同垂直邊 0.5)
    stepV = 4;
    for i = 1:stepV:numel(pct1)
        if pct1(i) > 0
            line(ax, [edg1(i) edg1(i)], [0 pct1(i)], 'Color','k', 'LineWidth',0.5);   % 全高垂直邊
        end
    end
    % 兩個分布的 mean 虛線(深藍=single、深紅=eighteen)
    mu0 = mean(err0);  mu1 = mean(err1);
    ml0 = xline(ax, mu0, '--', 'Color',[0.00 0.00 0.00], 'LineWidth',2.8);   % single mean 0.250（黑，高對比）
    ml1 = xline(ax, mu1, '--', 'Color',[0.00 0.60 0.00], 'LineWidth',2.8);   % eighteen mean 0.038（深綠，高對比）

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.015 .015],'TickDir','out');
    [xr_, xt_] = xlim_pick(SRC, max(err0));   % [MODIFIED] 用 single(fix) 的殘差定刻度;apdl 沿用定案值
    xlim(ax, xr_);  set(ax,'XTick', xt_);
    % [REVERTED 2026-08-03] 縱軸維持原作法（使用者拍板；同 plot_err_hist）：
    %   ylim 貼齊資料最大值(進位到 0.1)、tick 取 linspace 的 3 個內縮值。
    ymax = max([pct0 pct1]);  ytop = ceil(ymax*10)/10;  ylim(ax,[0 ytop]);
    yd = round(linspace(0, ytop, 5), 1);  set(ax,'YTick', yd(2:4));  % 3 內縮 tick(首末端點不標)
    % 起點 0 + 終點：只標數字、不畫 tick(左右角落補字)
    xr = xlim(ax);
    text(ax, xr(1), -0.022*ytop, sprintf('%g',xr(1)), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS,'FontWeight','bold','Clipping','off');
    text(ax, xr(2), -0.022*ytop, sprintf('%g',xr(2)), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS,'FontWeight','bold','Clipping','off');
    hLegR = bar(ax, NaN, NaN, 'FaceColor',cR, 'FaceAlpha',0.60, 'EdgeColor','k', 'LineWidth',1.2);   % 只給圖例：紅框加黑邊(實際紅 bar 仍無邊)
    % [MODIFIED] 圖例移到圖框外「上方」、**兩列**(NumColumns=2 為欄優先填充 → 條目 1,2 = 左欄兩列、
    %   3,4 = 右欄兩列，故 [h0 hLegR ml0 ml1] 恰為「第一列 Single、第二列 Eighteen」)。
    %   欄距不用白字墊 —— 下面會強制圖例寬度 = 座標軸寬度，兩欄自動平均散開。
    lg = legend([h0 hLegR ml0 ml1], {'Single parameter', 'Eighteen parameters', ...
                sprintf('Single mean = %.3f mT', mu0), sprintf('Eighteen mean = %.3f mT', mu1)}, ...
                'Interpreter','tex', 'Location','northoutside', 'NumColumns',2);
    lg.FontSize = 24;  lg.FontWeight = 'bold';  lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;
    xlabel(ax, '$\mathbf{Residual\;(mT)}$', ...
           'Interpreter','latex', 'FontSize',36);      % 絕對殘差軸標題(標準數學字體、36 粗)
    ylabel(ax, '$\mathbf{Percentage\;(\%)}$', ...
           'Interpreter','latex', 'FontSize',36);      % [ADDED] 百分比縱軸標題(同 err_hist 風格)
    ax.Toolbar.Visible = 'off';  hold(ax,'off');

    % [ADDED] 圖例框:①左右邊界對齊圖框(寬度=ax 寬、x 起點同 ax，兩欄自動平均散開)
    %                ②壓低 ax 上緣，讓圖例底線與圖框上緣「留明確間隙」、不相黏
    drawnow;
    axp = get(ax,'Position');  lgh = get(lg,'Position');  lgh = lgh(4);
    GAPN   = 0.022;                          % 圖例底線 ↔ 圖框上緣的間隙(normalized)
    newTop = 1 - lgh - GAPN - 0.006;         % ax 新上緣(頂端留 0.006 給圖例外框)
    axp(4) = newTop - axp(2);                % 只壓上緣、下緣不動
    set(ax, 'Position', axp);
    set(lg, 'Position', [axp(1), newTop + GAPN, axp(3), lgh]);

    % [MODIFIED 2026-08-03] 檔名一律明確標「求解器_半徑」（overlay 本身含 single+eighteen，不標模型）
    cstr = ''; if CONV, cstr = '_conv'; end
    out = fullfile(figdir, sprintf('err_hist_overlay%s_%s_R%d.png', cstr, lower(SRC), Rum));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function [c0, c1] = pick_bar_colors(Rum)
% [ADDED 2026-08-03] 長條配色每個取樣半徑一組（使用者拍板，與 plot_err_hist 的 pick_colors 同組）：
%   R150：藍(single) / 紅(eighteen)
%   其他(R300…)：紫 #7B52AB (single) / 橘 #E69F00 (eighteen)
%   ⚠ mean 虛線不隨之改——一律黑(single)/深綠(eighteen)，中性色在兩組長條上都看得清。
    if Rum == 150
        c0 = [0.10 0.35 1.00];      c1 = [0.85 0.10 0.10];
    else
        c0 = [0.482 0.322 0.671];   c1 = [0.902 0.624 0.000];
    end
end

% ============================================================================
function CAL = solver_path(SRC)
% [ADDED] 依 SRC 掛上對應分支、並「移除另一分支」——兩分支有同名函式，只 addpath 會被遮蔽。
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
    fprintf('[SRC=%s] model_config -> %s\n', SRC, which('model_config'));
end

% ============================================================================
function raw = load_raw(SRC, cfg)
% [ADDED] apdl = graded .dat / maxwell = .fld
    if strcmpi(SRC,'maxwell')
        raw = extract_maxwell_data(cfg, 'all', cfg.default_variant);
    else
        raw = extract_ansys_data(cfg, 'all', 'graded');
    end
end

% ============================================================================
function [xr, xt] = xlim_pick(SRC, maxE)
% [ADDED] 橫軸(殘差 mT)：apdl 沿用定案 fix 刻度；maxwell 依 single 殘差最大值自動等距 4 內縮 tick。
    if strcmpi(SRC,'apdl'), xr = [0 1.0];  xt = [0.2 0.4 0.6 0.8];  return; end
    % [FIXED 2026-08-03] 原本用 `while xr(2)<maxE, s=nice_step(s*1.3); end` 找上界 —— **無窮迴圈**：
    %   nice_step 會把 s*1.3 吸附回同一格（1→1.3→1、2.5→3.25→2.5、5→6.5→5），s 永遠不變。
    %   R=300 踩到此路徑 → MATLAB 空轉數小時、無輸出。改成在候選清單上直接往上找，保證終止。
    cand = [1 2 2.5 5 10];
    k0 = floor(log10(max(maxE, realmin)/5));   s = [];
    for k = k0:(k0+3)
        for c = cand
            if 5*c*10^k >= maxE, s = c*10^k; break; end
        end
        if ~isempty(s), break; end
    end
    if isempty(s), s = maxE/5; end             % 保險（理論上不會走到）
    xr = [0, 5*s];
    xt = (1:4)*s;
    % [ADDED 2026-08-03] 長尾資料時 5×步長 會讓右端留白過大（違反 figure-style「留白不可大於間距」）：
    %   上界 >1.25×maxE 才收緊 → 改「n 個等距 tick + 半格留白」。R150 等既有圖不受影響（比值僅 1.05）。
    if xr(2) > 1.25*maxE
        s2 = nice_step(maxE/3.5);
        n  = max(3, floor(maxE/s2));
        up = (n + 0.5)*s2;
        while up < maxE, n = n + 1;  up = (n + 0.5)*s2;  end   % n 遞增，必終止
        xr = [0, up];   xt = (1:n)*s2;
    end
end

function s = nice_step(x)
% [ADDED] nice 等距步長（1/2/2.5/5/10 × 10^k）
    k = floor(log10(x));   m = x/10^k;
    cand = [1 2 2.5 5 10];
    [~,i] = min(abs(cand - m));
    s = cand(i)*10^k;
end


% ---- 擬合 + 逐節點×激發 絕對殘差 |S·G − Bstack| (mT) ----
function [Pq, Bq, np] = conv_set(Rum, cfg, F, USE_BIAS, here)
% [ADDED 2026-08-13] 取該 R 的 N_c 收斂設計 → 產減量取樣點與場（校正用）。
    tri = conv_design(Rum*1e-6, cfg, F, USE_BIAS, here);
    [x,y,z,B] = sphere_grid_sample(Rum*1e-6, [], struct('frame','actuator','NRPT',tri));
    Pq = [x y z];   np = size(Pq,1);
    Bq = zeros(3*np, size(B,3));
    for j = 1:size(B,3), Bq(:,j) = reshape(B(:,:,j).', [], 1); end
end

% ============================================================================
function err = fit_abs_resid(P, Bstack, Pc_base, F, npts, USE_BIAS, Pc_cal, Bc_cal)
% [MODIFIED 2026-08-13] 可選第 7/8 引數：**用另一組點校正**（N_c 減量設計），
%   再拿 (P, Bstack)（該 R 內全部真實格點）當評估集。不給就是原行為（同一組點校正+評估）。
    if nargin < 7 || isempty(Pc_cal), Pc_cal = P;  Bc_cal = Bstack; end
    [e, l_hat, ~] = fitting(Pc_cal, Bc_cal, Pc_base, 0.5e-3, USE_BIAS);
    [~, ~, G]     = solve_current(l_hat, e, Pc_base, Pc_cal, Bc_cal, F);
    Pc = make_Pc(e, Pc_base);
    pbar = P / l_hat;  S = zeros(3*npts, 6);
    for k = 1:6, d = pbar - Pc(:,k).'; r3 = sum(d.^2,2).^1.5; S(:,k) = reshape((d./r3).', 3*npts, 1); end
    resid = S*G - Bstack;
    err = zeros(npts, 6);
    for j = 1:6, rj = reshape(resid(:,j), 3, npts); err(:,j) = sqrt(sum(rj.^2,1)).'; end
    err = err(:);
end

% ---- 電荷格 Pc = Pc_base + E(e)(與 solve_current 內一致；fix e=0 → Pc_base)----
function Pc = make_Pc(e17, Pc_base)
    if isempty(e17) || all(e17(:)==0), Pc = Pc_base; return; end
    E = zeros(3, 6);
    E(:,1) = e17(1:3);     E(:,2) = e17(4:6);
    E(:,3) = e17(7:9);     E(:,4) = e17(10:12);
    E(:,5) = e17(13:15);
    E(1,6) = e17(16);      E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);
    Pc = Pc_base + E;
end
