function plot_err_hist(USE_BIAS, SRC, Rsel_um)
% plot_err_hist -- long2016 半切六極「絕對殘差直方圖」(取樣半徑可調，預設 R=150µm)
% =========================================================================
%   graded / current、R=150µm 校正球。逐節點(每 node×每激發)絕對殘差 (mT):
%     |b_FEM − Sᵢ·ᴮĝ_I·K̄·I_j| = ||S·G − Bstack||_node   (S 用擬合後電荷位置 Pc)。
%   直方圖:亮藍填色 + 細黑邊、mean 虛線、圖例(sampling range + mean)、percentage 縱軸、無軸標題。
%   USE_BIAS=false → fix(單一 ℓ)→ err_hist.png；true → 18-param bias → err_hist_bias.png。
%   同時印 RMSPE = sqrt(Σε²/Σb²)·100(與 PDF 一致)。
% =========================================================================
    clc;
    if nargin < 1, USE_BIAS = false; end   % false = fix → err_hist.png；true = bias → err_hist_bias.png
    if nargin < 2, SRC = 'apdl';    end    % [ADDED] 'apdl'(ANSYS .dat) | 'maxwell'(.fld) → 檔名加 _maxwell
    if nargin < 3, Rsel_um = 150;   end    % [ADDED] 取樣球半徑 [µm]（原寫死 150；同 plot_err_hist_overlay）
    here   = fileparts(mfilename('fullpath'));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    CAL = solver_path(SRC);                                          % [ADDED] 依 SRC 切分支(內含 rmpath 另一分支防遮蔽)
    sstr = ''; if strcmpi(SRC,'maxwell'), sstr = '_maxwell'; end     % [ADDED] 來源後綴

    % ---- 前段 + R=150 fix 擬合 ----
    cfg = model_config('long2016_hexapole_halfcut','tip40um');
    raw = load_raw(SRC, cfg);                                        % [MODIFIED] apdl='graded' .dat / maxwell=.fld
    ad  = build_actuator_data(raw, cfg);   Pc_base = ad.Pc_base;
    F   = zeros(6, cfg.N_I);  for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end
    Rum = Rsel_um;                                                   % [MODIFIED] 取樣半徑可調（原寫死 150）
    [P, Bstack, npts] = cfg.select_ball(ad, Rum*1e-6);
    fprintf('取樣範圍 R <= %d um：npts = %d\n', Rum, npts);
    [e, l_hat, ~] = fitting(P, Bstack, Pc_base, 0.5e-3, USE_BIAS);
    [~, ~, G] = solve_current(l_hat, e, Pc_base, P, Bstack, F);

    % ---- 重建 S(用擬合後電荷位置 Pc = Pc_base + E(e)；fix 時 e=0 → Pc=Pc_base)、模型場、殘差 ----
    Pc = make_Pc(e, Pc_base);
    pbar = P / l_hat;  S = zeros(3*npts, 6);
    for k = 1:6, d = pbar - Pc(:,k).'; r3 = sum(d.^2,2).^1.5; S(:,k) = reshape((d./r3).', 3*npts, 1); end
    resid = S*G - Bstack;                         % 3npts×6

    % ---- 逐節點×激發 絕對殘差 |b_FEM − S·ĝ_I·K̄·I| (mT) ----
    %   resid = S*G − Bstack；S·G ≡ S·(ĝ_I·K̄_I·I)（G = ĝ_I·K̄_I·F 的重排），故此殘差即 |b_FEM − 模型|
    err = zeros(npts, 6);
    for j = 1:6
        rj = reshape(resid(:,j), 3, npts);   % 殘差向量 (mT)
        err(:,j) = sqrt(sum(rj.^2,1)).';     % 每節點×激發 殘差大小 (mT)
    end
    err = err(:);
    RMSPE = sqrt(sum(resid(:).^2) / sum(Bstack(:).^2)) * 100;
    mu = mean(err);
    fprintf('npts=%d  RMSPE=%.3f%%  mean|resid|=%.4f mT  max=%.4f mT\n', npts, RMSPE, mu, max(err));

    % ---- 直方圖(percentage 縱軸;亮藍 + 粗黑邊)----
    nb  = 180;  maxE = max(err);                       % 全範圍、密(依直方圖規則密度)
    edg = linspace(0, maxE, nb+1);
    [cnt, edg] = histcounts(err, edg);
    pct = cnt / numel(err) * 100;
    ctr = (edg(1:end-1) + edg(2:end)) / 2;

    FS = 28;
    fig = figure('Color','w','Position',[100 100 1100 790]);  ax = axes(fig);  hold(ax,'on');   % [MODIFIED] 加高補償外置圖例佔的縱向空間
    [cBar, cMean] = pick_colors(USE_BIAS, Rum);   % [MODIFIED] 配色依取樣半徑分組（見 pick_colors）
    hb = bar(ax, ctr, pct, 1, 'FaceColor',cBar, 'FaceAlpha',0.95, ...
             'EdgeColor','k', 'LineWidth',0.3);                       % + 細黑邊(密)
    ml = xline(ax, mu, '--', 'Color',cMean, 'LineWidth',3.0);         % mean 虛線

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.015 .015],'TickDir','out');
    [xr_, xt_] = xlim_pick(SRC, USE_BIAS, maxE);   % [MODIFIED] apdl 沿用定案值;maxwell 依實際殘差自動(等距 4 內縮 tick)
    xlim(ax, xr_);  set(ax,'XTick', xt_);
    % [REVERTED 2026-08-03] 縱軸維持原作法（使用者拍板：原本的刻度就可以）：
    %   ylim 貼齊資料最大值(進位到 0.1)、tick 取 linspace 的 3 個內縮值。
    %   ⚠ 曾改成「等距 nice 步長 + ylim 由 tick 反推」，雖然完全等距，卻讓長條縮小、上方空一大塊 → 已退回。
    ymax = max(pct);  ytop = ceil(ymax*10)/10;  ylim(ax,[0 ytop]);
    yd = round(linspace(0, ytop, 5), 1);  set(ax,'YTick', yd(2:4));  % 3 內縮 tick(首末端點不標)
    % 起始點 0 + 終點：只標數字、不畫 tick（左右角落補字，text）
    xr = xlim(ax);
    text(ax, xr(1), -0.022*ytop, sprintf('%g',xr(1)), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS,'FontWeight','bold','Clipping','off');   % 起點
    text(ax, xr(2), -0.022*ytop, sprintf('%g',xr(2)), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS,'FontWeight','bold','Clipping','off');   % 終點
    % [MODIFIED] 圖例移到圖框外「上方」、水平一列（比照 Bode_Surface 參考圖）
    %   兩則之間的間隔:用白色(隱形)字元墊寬第一則。MATLAB 會修掉尾端空白、
    %   dummy 條目又固定佔 ~105px 太寬;白字每個 'M' 約 +27px,可精準調(現用 MM ≈ +53px)。
    GAP = '{\color{white}MM}';
    lg = legend([hb ml], {[sprintf('Sampling range \\leq %d {\\mu}m', Rum) GAP], sprintf('Mean = %.3f mT', mu)}, ...
                'Interpreter','tex', 'Location','northoutside', 'Orientation','horizontal');
    lg.FontSize = 24;  lg.FontWeight = 'bold';  lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;
    xlabel(ax, '$\mathbf{Residual\;(mT)}$', ...
           'Interpreter','latex', 'FontSize',36);      % 絕對殘差軸標題(標準數學字體、36 粗)
    ylabel(ax, '$\mathbf{Percentage\;(\%)}$', ...
           'Interpreter','latex', 'FontSize',36);      % [ADDED] 百分比縱軸標題(同風格)
    ax.Toolbar.Visible = 'off';  hold(ax,'off');       % 刻度數字保留

    % [MODIFIED 2026-08-03] 檔名一律明確標「求解器_模型_半徑」，不留隱含預設（使用者拍板）
    out = fullfile(figdir, sprintf('err_hist_%s_%s_R%d.png', ...
                   lower(SRC), ternary_str(USE_BIAS,'eighteen','single'), Rum));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function s = ternary_str(cond, a, b)
    if cond, s = a; else, s = b; end
end

% ============================================================================
function [cBar, cMean] = pick_colors(USE_BIAS, Rum)
% [ADDED 2026-08-03] 每個取樣半徑一組配色（使用者拍板），組內 bar/mean 兩色互換：
%   R150：single = 亮藍 bar / 紅 mean； eighteen = 紅 bar / 深綠 mean（與 err_hist_overlay 對齊）
%   其他(R300…)：single = 紫 bar / 橘 mean； eighteen = 橘 bar / 紫 mean
%     紫 #7B52AB / 橘 #E69F00（Okabe-Ito 系）——色相與 R150 的藍/紅分得最開、色盲安全
    PUR = [0.482 0.322 0.671];      % #7B52AB
    ORG = [0.902 0.624 0.000];      % #E69F00
    if Rum == 150
        if USE_BIAS, cBar = [0.85 0.10 0.10];  cMean = [0.00 0.60 0.00];
        else,        cBar = [0.10 0.35 1.00];  cMean = [0.85 0.10 0.10];
        end
    else
        if USE_BIAS, cBar = ORG;  cMean = PUR;
        else,        cBar = PUR;  cMean = ORG;
        end
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
% [ADDED] apdl = graded .dat / maxwell = .fld（cfg.default_variant='maxwell'）
    if strcmpi(SRC,'maxwell')
        raw = extract_maxwell_data(cfg, 'all', cfg.default_variant);
    else
        raw = extract_ansys_data(cfg, 'all', 'graded');
    end
end

% ============================================================================
function [xr, xt] = xlim_pick(SRC, USE_BIAS, maxE)
% [ADDED] 橫軸(殘差 mT)範圍/刻度：apdl 沿用定案值；maxwell 依實際 max 殘差自動取等距 4 內縮 tick。
    if strcmpi(SRC,'apdl')
        if USE_BIAS, xr = [0 0.25];  xt = [0.05 0.10 0.15 0.20];
        else,        xr = [0 1.0];   xt = [0.2 0.4 0.6 0.8];   end
        return;
    end
    % [FIXED 2026-08-03] 同 plot_err_hist_overlay：原 `while xr(2)<maxE, s=nice_step(s*1.3); end`
    %   會無窮迴圈（nice_step 把 s*1.3 吸附回同一格：1→1.3→1、2.5→3.25→2.5、5→6.5→5）。
    %   改成在候選清單上直接往上找，保證終止。
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
    % [ADDED 2026-08-03] 同 overlay：上界 >1.25×maxE 才收緊，避免長尾資料右端留白過大。
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


% ---- 電荷格 Pc = Pc_base + E(e)（含 e6z 約束；與 solve_current 內 make_Pc 一致）----
function Pc = make_Pc(e17, Pc_base)
    if isempty(e17) || all(e17(:)==0)   % fix：無偏移
        Pc = Pc_base;  return;
    end
    E = zeros(3, 6);
    E(:,1) = e17(1:3);     E(:,2) = e17(4:6);
    E(:,3) = e17(7:9);     E(:,4) = e17(10:12);
    E(:,5) = e17(13:15);
    E(1,6) = e17(16);      E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);
    Pc = Pc_base + E;
end
