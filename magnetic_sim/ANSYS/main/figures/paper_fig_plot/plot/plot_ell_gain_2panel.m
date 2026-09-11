function plot_ell_gain_2panel(SUF, NOLINE, AX, PW)
% plot_ell_gain_2panel -- 降取樣（收斂點設計）結果的 l_hat 與 g_I_hat vs 取樣半徑 R
% =========================================================================
%   上下兩個子圖、共用橫軸（取樣半徑 R）：
%     上：l_hat   [um]     single（藍）與 eighteen（紅）兩條疊在一起
%     下：g_I_hat [mT/A]   同上
%   **只畫降取樣（conv）那一組**——即每個 R 各自用「收斂點設計」的少量取樣點校正的結果；
%   全格點那一組（ell_f*/gI_f*）不畫（那是 full_vs_conv_vs_R_*.png 的內容）。
%
%   資料直接讀 plot_full_vs_conv_vs_R.m 的快取 data/full_vs_conv_vs_R_maxwell.mat
%   （欄位 ell_c1/ell_c2、gI_c1/gI_c2、R_um），**不重跑任何擬合**。
%
%   **橫軸（取樣範圍 R）自 0 起**（使用者指定）：資料只有 R = 80~500，
%   故在 R=0 補一個**初始值 0** 的錨點（兩個子圖皆然），曲線自原點拉起。
%   縱軸連帶自 0 起，照 figure-style「自 0 起的軸上緣只留 8% 裕度」。
%
%   風格①粗體框圖：FS 36 粗體、box on、grid off、刻度奇數等距、曲線首末點貼齊左右框邊、
%   端點數字以 text 補、圖例照 figure-style 標準樣式（共用一個、置於最上方）。
%   配色：single 藍 / eighteen 紅（顏色 = 模型）。
%
%   輸出 → figures/paper_fig/Section2_E/ell_gain_2panel_maxwell.png（覆蓋迭代）
% =========================================================================
    clc;
    % [ADDED 2026-08-27] SUF：快取與輸出檔名的共同後綴，讓不同判準/模型的版本並存。
    %   ''            -> full_vs_conv_vs_R_maxwell.mat        -> ell_gain_2panel_maxwell.png
    %   '_fullref'    -> ..._fullref.mat                      -> ..._fullref.png
    %   '_fullref_zhi'-> ..._fullref_zhi.mat                  -> ..._fullref_zhi.png
    if nargin < 1 || isempty(SUF), SUF = ''; end
    % [ADDED 2026-08-28 使用者指定] NOLINE：true = **只畫點、不用線連**（marker only），
    %   輸出檔名加 '_pts' 後綴、不覛掛有線版。資料、軸範圍、刻度、圖例全部相同，
    %   只把 LineStyle 拿掉（含 R=0 那個錨點，不畫線時它會成為孤點）。
    if nargin < 2 || isempty(NOLINE), NOLINE = false; end
    % [ADDED 2026-08-31] AX: optional [ax_top ax_bottom] to draw into. Given, the
    %   panels are placed in the caller's figure and nothing is exported here --
    %   this is how plot_ell_gain_err_merged reuses them without duplicating code.
    if nargin < 3, AX = []; end
    if NOLINE, LS1 = 'o';  LS2 = 's';  else, LS1 = '-o'; LS2 = '-s'; end
    here   = fileparts(fileparts(mfilename('fullpath')));       % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    df = fullfile(here, 'data', sprintf('full_vs_conv_vs_R_maxwell%s.mat', SUF));
    assert(exist(df,'file')==2, '找不到 %s —— 先跑 plot_full_vs_conv_vs_R(true)', df);
    D = load(df);
    R = D.R_um(:).';
    % [ADDED 2026-08-27] 丟掉整列都是 NaN 的 R（球內點太少、掃描時直接跳過，例如 R=20）。
    ok = isfinite(D.ell_c1(:).') & isfinite(D.ell_c2(:).') & ...
         isfinite(D.gI_c1(:).')  & isfinite(D.gI_c2(:).');
    if ~all(ok)
        fprintf('  略過無資料的 R = %s um\n', num2str(R(~ok)));
        R = R(ok);
        for f = {'ell_c1','ell_c2','gI_c1','gI_c2','n_c1','n_c2','ell_f1','ell_f2','gI_f1','gI_f2'}
            if isfield(D,f{1}), v = D.(f{1})(:).';  D.(f{1}) = v(ok); end
        end
    end

    fprintf('R  = %s\n', num2str(R,'%6d'));
    fprintf('l_hat  single  : %.1f ~ %.1f um   (N_c %d ~ %d)\n', min(D.ell_c1), max(D.ell_c1), min(D.n_c1), max(D.n_c1));
    fprintf('l_hat  eighteen: %.1f ~ %.1f um   (N_c %d ~ %d)\n', min(D.ell_c2), max(D.ell_c2), min(D.n_c2), max(D.n_c2));
    fprintf('g_I    single  : %.3f ~ %.3f mT/A\n', min(D.gI_c1), max(D.gI_c1));
    fprintf('g_I    eighteen: %.3f ~ %.3f mT/A\n', min(D.gI_c2), max(D.gI_c2));

    % ---- 樣式 ----
    % [ADDED 2026-09-01 使用者要求] PW = 這張圖在論文裡的**最終擺放寬度 [in]**
    %   （[] = 維持原本 12 in 畫布，輸出逐位不變）。給了 PW 就以最終尺寸重畫：
    %   畫布 = 4*PW 見方、字級/線寬 = 目標點數 x4 -> 縮到 PW 時軸字 8 pt、圖例 7 pt、
    %   框線 1.0 pt、資料線 0.8 pt、marker 3 pt。與 plot_err_hist_shell 的 PW 同一把尺。
    if nargin < 4 || isempty(PW), PW = []; end
    PAPER = ~isempty(PW);
    % [MODIFIED 2026-09-01 套繪圖規則 1/2/3] 刻度 FS / 軸標題 FSLAB / 圖例 FSLEG 拆開；
    %   線寬與 marker 隨刻度字級 36->60 同量級放大，框線 2.5->5.0。
    FS = 60;  FSLAB = 36;  FSLEG = 45;  LW = 5.0;  MS = 12;  LWBOX = 5.0;
    % [MODIFIED 2026-09-01 使用者指定] 起始點由粉色改**綠色**。
    CINIT = [0.00 0.65 0.20];      % 起始點（x=0）專用色
    MSP   = 22;                    % [ADDED 2026-09-01 使用者反映太小] 粉點 marker 尺寸
    if PAPER, FS = 32;  LW = 3.2;  MS = 12;  LWBOX = 4.0; end
    c1 = [0.05 0.10 0.95];      % single   藍
    c2 = [0.85 0.10 0.10];      % eighteen 紅

    % 手動 axes 定位（**不用 tiledlayout**：它不允許事後設 Position，
    %   圖例就無法照 figure-style 切齊座標框）。
    % [MODIFIED 2026-08-31 使用者指定] 改回**上下排列**：上 = l_hat、下 = g_I，兩格等寬等高，
    %   各自有完整的橫軸刻度與軸標題（上格的 xlabel 佔用兩格之間的 GAPV）。
    if isempty(AX)
        % [MODIFIED 2026-08-31 使用者指定] 正方形畫布。這張與 err_hist_conv_* 會在論文裡
        %   左右並排，兩張長寬比相同才能等寬等高擺在一起、誰都不會被壓扁。
        %   用英吋而非像素：像素尺寸超過螢幕時 MATLAB 會靜默縮小視窗，匯出比例就跑掉。
        CANV = 16.0;   if PAPER, CANV = 4*PW; end
        fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
        % [MODIFIED 2026-09-01] 版面常數重算：刻度字級 36->60 之後，3 位數的 y 刻度
        %   數字約 1.3 in 寬、加上 ylabel 與間隙需要約 3 in 左緣（原本 1.62 in 會把
        %   ylabel 整個擠出畫布）；右緣也要留半個 '500' 的寬度免得端點數字被裁。
        % [MODIFIED 2026-09-01] 畫布 12 -> 14.5 in（與單面板 axsh_* 同尺寸，60pt 的
        %   刻度數字相對大小才一致）。左緣 = 實測 TightInset 2.36 in -> 0.17；
        %   右緣留半個 '500'。框寬 0.79x14.5 = 11.5 in **大於圖例自然寬 10.66 in**
        %   -> 圖例終於可以切齊座標框（12 in 畫布時框只有 8.4 in，切不齊只能置中）。
        % [MODIFIED 2026-09-01] 45 pt 之下 'Single parameter'(5.09 in) +
        %   'Eighteen parameters'(6.19 in) = 11.28 in，加標記與間距約需 12.5 in。
        %   左緣 2.36 in（實測 TightInset）+ 右緣 0.55 in（'500' 半寬）
        %   -> 畫布至少 15.4 in，取 16.0；框寬 0.815x16 = 13.04 in，留有餘裕。
        AXL  = 0.150;   AXW = 0.800;                   % 左緣 / 每格寬（單欄）
        % [ADDED 2026-09-01] PAPER 模式的版面常數。畫布由 12 in 縮到 4*PW=6.88 in
        %   但字級只從 36 降到 32 -> 字相對大 1.63 倍，原本 13.5% 的左緣塞不下
        %   「ℓ̂ (micro meter)」+「900」（實測整條 ylabel 被切出畫布）。重算：
        %     左緣 = ylabel 行高 0.71 in + 刻度字寬 0.79 in + 留白 -> 1.6 in = 23.5%
        %     右緣 = 末端刻度「500」半寬 0.40 in + 留白 -> 0.45 in = 6.5%
        %     GAPV = 上格 xlabel 0.7 in + 其刻度數字 0.55 in -> 1.27 in = 18.5%
        % [MODIFIED 2026-09-01] 拿掉軸標題後，左緣只需放刻度數字（'900' 約 0.73 in
        %   = 10.6% 圖寬），右緣只需末端刻度半寬 -> 繪圖框大幅放大。
        if PAPER, AXL = 0.145;  AXW = 0.800; end
        % [MODIFIED 2026-08-31 使用者指定] 繪圖區總高（上格上緣→下格下緣）= 圖幅的 65%，
        %   與 err_hist_conv_*_sq.png 的 1348/2073 = 65.0% 對齊 —— 兩張在論文裡等寬並排時
        %   座標框高度才會一樣。補白無法把直方圖的 65% 撐到原本的 79%（那是它框高/框寬的
        %   固有比值），只能反過來把這張縮到 65%。
        SPANV = 0.650;                                 % 繪圖區總高（正規化）
        GAPV = 0.105;                                  % 兩格垂直間距（放上格的 xlabel）
        % [MODIFIED 2026-09-01] PAPER 模式：**上格不重複畫橫軸**（兩格共用同一條
        %   取樣半徑軸，只有下格標刻度與軸標題）。理由是垂直預算 ——
        %   1.72 in 高扣掉圖例 2 列 0.32 + 兩格各自的 xlabel+刻度 0.60 + 留白 0.08
        %   只剩 0.72 in 給兩個繪圖框（每格 9 mm），字必然互撞、刻度會被自動轉 45 度。
        %   砍掉上格那份 0.30 in 之後每格回到 0.51 in。BOT 也不再上下置中
        %   （上方要留給圖例、下方只要留給 xlabel，兩者需求不同）。
        PAPNOX = false;
        if PAPER, SPANV = 0.720;  GAPV = 0.030; end
        % [ADDED 2026-09-01 規則 8] 上下堆疊圖：**上格水平軸不標軸標題、也不標刻度數字**
        %   （刻度線照畫）。兩格共用同一條 x 軸，下格標一次就夠 -> 所有模式都套。
        %   圖例由一列變兩列（多了 Initial value），實測高 0.120、寬 0.888（12 in 畫布）。
        PAPNOX = true;   GAPV = 0.045;   SPANV = 0.720;
        BOT  = 0.130;                                  % 下方留 x 刻度數字(60pt)+xlabel
        if PAPER, BOT = 0.100; end                     % 下方只需 x 刻度數字的高度
        TOP1 = BOT + SPANV;
        PH   = (SPANV - GAPV) / 2;                     % 每格高
        BOT2 = BOT + PH + GAPV;                        % 上格下緣
        ax1 = axes(fig, 'Position', [AXL, BOT2, AXW, PH]);   % 上：l_hat
        ax2 = axes(fig, 'Position', [AXL,  BOT, AXW, PH]);   % 下：g_I
    else
        ax1 = AX(1);   ax2 = AX(2);   fig = ancestor(ax1,'figure');   PAPNOX = true;   % [規則8]
    end

    % [MODIFIED 2026-08-14] 每個 R 都是**真的跑過一次校正**（快取 full_vs_conv_vs_R_maxwell.mat，
    %   由 plot_full_vs_conv_vs_R 逐 R 搜尋收斂設計後擬合）。不補任何假錨點。
    %   橫軸自 0 起（使用者指定）；R=0 球內 0 點無法擬合，故 0 ~ 第一個 R 之間沒有資料點。
    % [MODIFIED 2026-08-14] 曲線自 **R=0 的起始點**拉起（使用者指定）：
    %   l_hat 的起始點 = 擬合初值 **500 um**；g_I 的起始點 = **0**。
    %   R>0 的每一點都是真的跑過一次校正（R = 20:20:500）。
    Rp = [0,   R];
    e1 = [500, D.ell_c1];   e2 = [500, D.ell_c2];
    g1 = [0,   D.gI_c1];    g2 = [0,   D.gI_c2];
    % 水平軸：整數刻度（使用者指定）
    XL = [0 R(end)];
    % [MODIFIED 2026-08-27 使用者指定] **水平軸起點不標數字**（拿掉 0 這根刻度），
    %   縱軸起點才標 —— 原點角落只留一個數字（縱軸的），不再兩個擠在一起。
    % [MODIFIED 2026-09-01 規則 4/5] 內部刻度四根 100/200/300/400（等距 100，與兩端
    %   0 / 500 的間距同為 100）；**起點與終點都改用 text 標數字、不畫刻度線**。
    XT = 100:100:400;

    % ===================== 上：l_hat =====================
    hold(ax1,'on');
    % [MODIFIED 2026-09-01 使用者指定] 第一個點（x=0 的擬合初值）**不分 single / eighteen**，
    %   兩條系列的該點 marker 用 MarkerIndices 略過，改畫一顆粉色點代表 Initial value。
    %   線仍從 x=0 拉起（MarkerIndices 只影響 marker，不影響線）。
    MIDX = 2:numel(Rp);
    h1 = plot(ax1, Rp, e1, LS1, 'Color',c1, 'LineWidth',LW, 'MarkerSize',MS, ...
              'MarkerFaceColor',c1, 'MarkerIndices',MIDX, 'Clipping','off');
    h2 = plot(ax1, Rp, e2, LS2, 'Color',c2, 'LineWidth',LW, 'MarkerSize',MS, ...
              'MarkerFaceColor','w', 'MarkerIndices',MIDX, 'Clipping','off');
    hp = plot(ax1, Rp(1), e1(1), 'o', 'Color',CINIT, 'MarkerFaceColor',CINIT, ...
              'MarkerSize',MSP, 'LineStyle','none', 'Clipping','off');
    % [MODIFIED 2026-08-14] 使用者指定：l_hat 縱軸 **自 500 起**，且刻度要把曲線夾在中間
    %   （資料 797~891 → 上下各有 700 / 900 兩個刻度）
    % [ADDED 2026-08-27] 縱軸範圍**依模型分支**。原本兩軸都寫死長飛的值，志鵬的
    %   g_I 最大 133.2 mT/A（R=40 退回全格點的近簡併值）遠超 [0 30] -> 會畫到框外。
    %   長飛的值維持逐字不變（使用者今日已逐輪定案）。
    isZhi  = contains(SUF, 'zhi');
    % [ADDED 2026-08-27] SUF='' 是舊 step 判準的快取，含 R=20 的 g_I = 40.13，
    %   g_I 軸上緣必須 > 40（用 fullref 的 [0 30] 會畫到框外）。
    isStep = isempty(SUF);
    % [MODIFIED 2026-08-27 使用者指定] l_hat 軸**沿用長飛的刻度間距 100**（跨圖同尺度、
    %   可直接比大小），但志鵬**上緣壓到 800**（資料最大 751.4）：
    %     長飛 500:100:900、上緣 930   志鵬 500:100:700、上緣 800
    %   志鵬的上緣 800 = 末刻度 700 + 100 -> **上方留白正好等於一個刻度間距**。
    %   g_I 軸無法共用（志鵬最大 133.2，是長飛的 4.8 倍）。
    if isZhi
        % 上緣 800 + **5 根 tick**：500/560/620/680/740，步長 60，
        % 上方留白 800-740 = 60 = 一個刻度間距。資料最大 751.4，不被切。
        YL1 = [500 800];   YT1 = 560:60:740;      % 4 根：[500,800] 平分 5 段
    else
        % [MODIFIED 2026-09-01 規則 4/5] 500 不再是刻度；[500,1000] 平分 5 段
        %   -> 內部四根 600/700/800/900，與框底 500、框頂 1000 的間距同為 100。
        % [MODIFIED 2026-09-01 使用者指定四根] 畫布放大到 16 in 之後每格高 5.4 in，
        %   四根刻度間距 1.08 in > 數字高 0.83 in，不再疊字（12 in 時只有 0.76 in 會撞）。
        %   [500,1000] 平分 5 段 -> 600/700/800/900，框底 500 與框頂 1000 都不標。
        YL1 = [500 1000];  YT1 = 600:100:900;
    end
    style_panel(ax1, XL, XT, YL1, YT1, FS, LWBOX);
    % [MODIFIED 2026-09-01 使用者指定] PAPER(_col) 模式**不畫軸標題**，只留刻度數字
    %   與圖例。在 PW=1.72 in 下，'Sampling range (micro meter)' 這串本身就有 1.71 in
    %   長 —— 等於整張圖的寬度，硬放只會被裁掉兩端並把 x 刻度擠成 45 度；縱向的
    %   'l-hat (micro meter)' 0.98 in 也比每格框高 0.49 in 還長，兩格的 ylabel 會疊在一起。
    if ~PAPER
        ylabel(ax1, '$\mathbf{\hat{\ell}\;(micro\;meter)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    end
    if PAPNOX                                     % 上格與下格共用橫軸：不重複標
        set(ax1, 'XTickLabel', {});
    else
        if ~PAPER
            xlabel(ax1, '$\mathbf{Sampling\;range\;(micro\;meter)}$', 'Interpreter','latex', 'FontSize',FSLAB);
        end
    end
    hold(ax1,'off');

    % ===================== 下：g_I_hat =====================
    hold(ax2,'on');
    % [MODIFIED 2026-08-14] g_I 縱軸自 **0** 起（起始點也是 0），刻度用**整數**
    %   （使用者指定：不要小數點）。上緣見下方 2026-08-17 的修訂。
    % [MODIFIED 2026-08-14] **所有點都畫**（使用者指定）—— 含 single 在 R=20/40 的
    %   小球近簡併值（六顆電荷近簡併、最小平方解的 G 下極對角整層翻負）。
    % [MODIFIED 2026-08-17] solve_current 改成 ĝ_I=(6/5)|H_I(1,1)| → R=20/40 由
    %   −40.13 / −14.85 變成 **+40.13 / +14.85**，資料範圍改為 0 ~ 40.13。
    %   縱軸改自 0 起（figure-style「自 0 起的軸上緣只留 8% 裕度」）：
    %   top = 1.08×40.13 = 43.3；刻度取**整數、奇數個、等距**的 8:8:40。
    % [MODIFIED 2026-08-27 使用者指定] 上緣 43.4 -> **32**、**內部 4 根 tick、上緣不標數字**。
    %   43.4 是當年為了容納 R=20 的 40.13 而設；R=20 的全格點只有 2 點、現已整點略過，
    %   實際資料最大是 R=40 的 27.606（single 退回全格點的值）-> 上緣 32 已足夠。
    %   照 figure-style 縱軸規則：[0,T] 平分 N+1 段 -> T = 5s（4 根內部刻度）。
    %   要「整數刻度」就得讓 T 是 5 的整數倍：T=32 -> s=6.4（6.4/12.8/19.2/25.6，難讀）
    %   -> 使用者要求換掉小數，取最接近的乾淨解 **T = 30、s = 6 -> 6/12/18/24**。
    %   資料最大 27.606（R=40 single 退回全格點）< 30，不會被切；填充率 92%（優於 T=32 的 86%）。
    %   30 不進 YTick -> 上緣自然無刻度線、無數字。
    % [ADDED 2026-08-28 使用者指定] axes-shell 掃描（SUF 含 _axsh）的 g_I 最大只有
    %   10.83（長飛）—— 用 fullref 的 [0 30] 會讓上方空掉三分之二。
    % [MODIFIED 2026-08-28 使用者指定] **中間要四根 tick**：原本 [0 12] + (0:3)*3 的
    %   0/3/6/9 只有三根在中間。改成 [0,T] 平分 5 段（T = 5s），s 取**整數 3**
    %   -> 刻度 0/3/6/9/12、上緣 15 不標數字，中間恰好四根。
    %   （s=2.5 -> 2.5/5/7.5/10、上緣 12.5，填充率 87% 比較滿，但帶小數點，
    %     與本 panel 2026-08-14 拍板的「刻度用整數、不要小數點」衝突，故不採。）
    if contains(SUF, 'axsh') && ~isZhi
        YL2 = [0 15];   YT2 = (1:4)*3;    % 4 根：3/6/9/12，上緣 15 不標
    elseif contains(SUF, 'axsh')
        % 志鵬 axes-shell：g_I 最大 37.41（fullref 那版是 133.2，上緣 150 會空掉四分之三）。
        % 同樣 [0,T] 平分 5 段、中間四根整數刻度：s=9 -> 0/9/18/27/36、上緣 45 不標數字。
        YL2 = [0 50];   YT2 = (1:4)*10;   % 4 根：10/20/30/40，上緣 50 不標
    elseif isZhi
        YL2 = [0 150];  YT2 = (1:4)*30;   % [規則4] 0 不再是刻度            % 志鵬 g_I 範圍 0(錨點)~133.2
    elseif isStep
        YL2 = [0 45];   YT2 = (1:4)*9;    % [規則4] 0 不再是刻度             % 舊 step 快取：資料最大 40.13（R=20）
    else
        YL2 = [0 30];   YT2 = (1:4)*6;    % [規則4] 0 不再是刻度
    end
    plot(ax2, Rp, g1, LS1, 'Color',c1, 'LineWidth',LW, 'MarkerSize',MS, ...
         'MarkerFaceColor',c1, 'MarkerIndices',MIDX, 'Clipping','off');
    plot(ax2, Rp, g2, LS2, 'Color',c2, 'LineWidth',LW, 'MarkerSize',MS, ...
         'MarkerFaceColor','w', 'MarkerIndices',MIDX, 'Clipping','off');
    plot(ax2, Rp(1), g1(1), 'o', 'Color',CINIT, 'MarkerFaceColor',CINIT, ...
         'MarkerSize',MSP, 'LineStyle','none', 'Clipping','off');
    style_panel(ax2, XL, XT, YL2, YT2, FS, LWBOX);
    if ~PAPER
        ylabel(ax2, '$\mathbf{{}^{B}\hat{g}_{I}\;(mT/A)}$', 'Interpreter','latex', 'FontSize',FSLAB);
        xlabel(ax2, '$\mathbf{Sampling\;range\;(micro\;meter)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    end

    % [MODIFIED 2026-09-01 規則 4] 水平軸**起點 0 與終點 500 都標數字、都不畫刻度線**
    %   （用 text 補在框下方）；縱軸起點與終點都不標。只加在**下格** —— 上格依規則 8
    %   連刻度數字都不標，端點數字自然也不能留。
    yr2 = ylim(ax2);   yoff2 = yr2(1) - 0.022*diff(yr2);
    for xv = XL
        text(ax2, xv, yoff2, sprintf('%g',xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    hold(ax2,'off');

    % ===================== 共用圖例（置於最上方）=====================
    % [FIXED 2026-08-31] 圖例掛在**位於目標位置的隱形 axes** 上，不掛 ax1。
    %   掛 ax1 時 MATLAB 的版面管理會在每次 render 把圖例的 y 拉回 'northoutside' 的
    %   自動值（0.799），正好跨在上格上框線（0.825）上把它蓋掉 —— 設 Location='none'、
    %   反覆設 Position、匯出前再設一次都無效（x 會照設定走，y 一定被拉回）。單一 axes
    %   的最小重現不會發生，是雙 axes 版面才有的行為。
    %   掛在隱形 axes 上就沒有可拉回的錨點，位置由該 axes 的 Position 決定。
    % [ADDED 2026-09-01] PAPER 模式圖例改**一欄兩列**：兩則文字加標記一行需要約
    %   2.9 in，而最終面板只有 PW=1.72 in -> 一行絕對放不下（實測溢出畫布兩側）。
    %   保留原本用詞，只改排列。
    % [ADDED 2026-09-01 使用者指定] 圖例多一則說明粉色起始點；NumColumns=2 之下
    %   三則會排成「兩則一列 + 一則一列」，正好是使用者要的「多加一行」。
    %   ⚠ 代理 handle 用**面板一那顆真的粉點 hp**（在 hold on 期間畫的）——
    %     在此處另外 plot 一顆會因為 ax1 已 hold off 而**清空整個座標軸**（踩過）。
    % [MODIFIED 2026-09-01 使用者指定] Single 與 Eighteen 要在**同一列**。
    %   MATLAB 的 legend 是**欄優先**填格：NumColumns=2、三則 ->
    %   第1欄放前兩則、第2欄放第三則。故順序給 [Single, Initial, Eighteen]，
    %   顯示出來第一列就是 Single | Eighteen，第二列是 Initial value。
    % [REWRITTEN 2026-09-01 使用者指定「文字與說明擺到中間」]
    %   MATLAB 的 legend **無法控制每列的水平對齊** —— 第二列只有一則時一定靠左，
    %   右邊留一大塊空白。故改成**手繪圖例**：隱形座標軸 + box on 當外框，
    %   逐列量文字寬度後把該列整組置中。順帶擺脫 legend 的兩個老問題
    %   （northoutside 會壓縮 ax1、以及不肯把框縮到比內容窄）。
    a1p = get(ax1,'Position');
    LFT = a1p(1);   SPAN = a1p(3);   TOPREF = a1p(2) + a1p(4);
    GAPN = 0.022;                                 % 上格框線與圖例之間的淨空
    NROW = 2;                                     % 第1列兩則、第2列一則
    lgh  = (NROW*1.15*FSLEG/72 + 0.25) / CANV;    % 列高 x 列數 + 內距（英吋 -> 正規化）
    lgpos = [LFT, TOPREF + GAPN, SPAN, lgh];      % 左右緣切齊座標框

    LG(1) = struct('mk','o', 'col',c1,   'face',c1,   'msz',MS,  'line',~NOLINE, ...
                   'txt','Single parameter');
    LG(2) = struct('mk','s', 'col',c2,   'face','w',  'msz',MS,  'line',~NOLINE, ...
                   'txt','Eighteen parameters');
    LG(3) = struct('mk','o', 'col',CINIT,'face',CINIT,'msz',MSP, 'line',false, ...
                   'txt','Initial value');
    ROWOF = [1 1 2];
    axL = draw_legend_box(fig, lgpos, LG, ROWOF, FSLEG, LW, LWBOX);

    if ~isempty(AX), return; end                  % 畫進呼叫端的座標軸：不自己出檔


    psuf = '';  if NOLINE, psuf = '_pts'; end
    if PAPER, psuf = [psuf '_col']; end            % col = 以單欄最終尺寸重畫
    out = fullfile(figdir, sprintf('ell_gain_2panel_maxwell%s%s.png', SUF, psuf));
    % [MODIFIED 2026-08-31] print 而非 exportgraphics：後者會把畫布周圍的空白**裁掉**，
    %   所以 12x12 吋的正方形畫布匯出後並不是正方形（實測 0.96 / 1.11）。這兩張要在
    %   論文裡等寬並排，長寬比必須一致，故保留整張畫布。
    CANV = 16.0;   if PAPER, CANV = 4*PW; end
    RES  = 200;  if PAPER, RES  = 300; end
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
    % [ADDED 2026-09-01] 匯出前再鎖一次圖例位置：建立當下 MATLAB 會照設定擺（已實測），
    %   但後續 render 仍可能把它往下挪。這裡 drawnow 後重設，並印出實際值供核對。
    drawnow;   set(axL, 'Position', lgpos);   drawnow;
    fprintf(['  [legend] set y=%.4f  actual y=%.4f  panel top=%.4f' newline], ...
            lgpos(2), getfield(get(axL,'Position'),{2}), TOPREF);
    print(fig, out, '-dpng', sprintf('-r%d', RES));
    fprintf('wrote %s\n', out);
    close(fig);
end

% ============================================================================
function style_panel(ax, XL, XT, YL, YT, FS, LWBOX)
    box(ax,'on');  grid(ax,'off');
    set(ax, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, 'TickLength',[.02 .02]);
    xlim(ax, XL);   set(ax, 'XTick', XT);
    ylim(ax, YL);   set(ax, 'YTick', YT);
    ax.Toolbar.Visible = 'off';
end

% ============================================================================
function [lim, tk] = axlim_from_zero(maxv)
% 自 0 起的軸（figure-style）：下緣固定 0、上緣只留 8% 裕度；
%   刻度間距由**資料範圍**取 nice（不隨上緣壓縮而變小）；刻度數取不超出上緣的最大奇數。
    cand = [1 2 2.5 3 4 5 10];
    x = maxv/4;   k = floor(log10(x));
    s = cand(find(cand*10^k >= x, 1)) * 10^k;
    top = 1.08 * maxv;
    n = floor(top/s);   if mod(n,2)==0, n = n - 1; end
    lim = [0 top];   tk = (1:n)*s;
end

% ============================================================================
function [lim, tk] = axlim_span(lo, hi)
% 涵蓋 lo~hi 的軸：奇數個等距 nice 刻度、兩端留白 = 刻度間距（figure-style）。
    cand = [1 2 2.5 3 4 5 10];
    mid  = (lo+hi)/2;   rng_ = max(hi-lo, realmin);
    for n = [3 5]
        for k = floor(log10(rng_/(n+1))) + (0:4)
            for c = cand
                s   = c*10^k;
                ctr = round(mid/s)*s;
                t   = ctr + (-(n-1)/2 : (n-1)/2)*s;
                L   = [t(1)-s, t(end)+s];
                if lo >= L(1)+0.15*s && hi <= L(2)-0.15*s
                    lim = L;   tk = t;   return;
                end
            end
        end
    end
    s = rng_/4;   t = mid + (-1:1)*s;   lim = [t(1)-s, t(end)+s];   tk = t;
end

% ============================================================================
function tk = ticks_lin_inner(lo, hi, n)
% 線性軸：n 個（奇數）等距 nice 刻度，嚴格落在 (lo,hi) 內部、不含端點。
%   n 是目標個數；湊不出就依序退到 n+2 / n-2 / n+4（都維持奇數）。
    cand = [1 1.5 2 2.5 3 4 5 10];   % 含 1.5：0~500 才湊得出 3 個等距（150/300/450）
    for nn = [n, n+2, n-2, n+4]
        if nn < 1, continue; end
        k0 = floor(log10(max(hi-lo, realmin)/(nn+1)));
        for k = k0:(k0+3)
            for c = cand
                s = c*10^k;
                t = s*ceil((lo + 0.5*s)/s) : s : (hi - 0.5*s);
                if numel(t) == nn, tk = t;  return; end
            end
        end
    end
    tk = lo + (1:n)/(n+1)*(hi-lo);
end

% ============================================================================
function axL = draw_legend_box(fig, pos, LG, ROWOF, FS, LW, LWBOX)
% 手繪圖例：隱形座標軸 + box on 當外框；每一列的「標記 + 文字」整組**水平置中**。
%   LG    : struct 陣列，欄位 mk/col/face/msz/line/txt
%   ROWOF : 每則屬於第幾列（1 = 最上列）
    axL = axes(fig, 'Position', pos, 'XLim',[0 1], 'YLim',[0 1], ...
               'Box','on', 'XTick',[], 'YTick',[], ...
               'XColor','k', 'YColor','k', 'LineWidth',LWBOX);
    drawnow;                                        % 必要：未 render 前 text Extent 會偏小
    hold(axL,'on');
    NR   = max(ROWOF);
    axw_in = pos(3) * get(fig,'Position')*[0;0;1;0];   % 本圖例框的實際寬度 [in]
    % [MODIFIED 2026-09-01 使用者：線段再長一點] 標記樣本 0.35 -> 0.50 in；
    %   預算從 padding / 間距擠出來（0.12->0.10、0.30->0.25），框寬同步 0.800->0.808。
    MW0   = 0.50 / axw_in;   % 標記樣本（線段）寬度 0.50 in -> 資料單位
    PADT0 = 0.10 / axw_in;   % 標記與文字之間
    GAPI0 = 0.25 / axw_in;   % 同列兩則之間
    for r = 1:NR
        idx = find(ROWOF == r);
        w   = zeros(1, numel(idx));
        for k = 1:numel(idx)                        % 先量文字寬（畫在視野外再刪）
            t = text(axL, 0, -1, LG(idx(k)).txt, 'FontSize',FS, 'FontWeight','bold');
            e = get(t,'Extent');   w(k) = e(3);   delete(t);
        end
        % 文字寬度是量出來的、不能改；裝飾（標記樣本 + 間距）若讓整列超出框寬，
        % 就等比縮小裝飾直到塞得下。不縮的話 x 會變負 -> 兩端被裁。
        % ⚠ 但**標記樣本不能縮到比 marker 本身還短**，否則線段會被 marker 完全蓋住、
        %   圖例看起來只剩一顆點（2026-09-01 踩過：MW 被縮到 0.14 in、marker 0.17 in）。
        %   故 MW 設下限 = marker 直徑 + 0.18 in，縮排只從 PADT / GAPI 取。
        MW = MW0;   PADT = PADT0;   GAPI = GAPI0;
        dec = numel(idx)*(MW+PADT) + (numel(idx)-1)*GAPI;
        if sum(w) + dec > 0.98
            MWmin = 0;
            for k = 1:numel(idx)                        % 有線段者的標記樣本下限
                if LG(idx(k)).line
                    MWmin = max(MWmin, (LG(idx(k)).msz/72 + 0.18)/axw_in);
                end
            end
            fix_ = numel(idx)*max(MW, MWmin);            % 不可縮的部分
            var_ = numel(idx)*PADT + (numel(idx)-1)*GAPI;% 可縮的部分
            kf   = max(0.05, (0.98 - sum(w) - fix_) / max(var_, eps));
            MW   = max(MW, MWmin);   PADT = PADT*kf;   GAPI = GAPI*kf;
            dec  = numel(idx)*(MW+PADT) + (numel(idx)-1)*GAPI;
        end
        tot = sum(w) + dec;
        x   = max(0.005, (1 - tot)/2);              % <- 整列置中
        assert(tot <= 0.995, 'legend row %d 放不下（tot=%.3f）—— 需加寬框或縮字級', r, tot);
        y   = 1 - (r - 0.5)/NR;
        for k = 1:numel(idx)
            it = LG(idx(k));
            if it.line
                plot(axL, [x, x+MW], [y y], '-', 'Color',it.col, 'LineWidth',LW);
            end
            plot(axL, x+MW/2, y, it.mk, 'Color',it.col, 'MarkerFaceColor',it.face, ...
                 'MarkerSize',it.msz, 'LineWidth',LW, 'LineStyle','none');
            text(axL, x+MW+PADT, y, it.txt, 'FontSize',FS, 'FontWeight','bold', ...
                 'VerticalAlignment','middle');
            x = x + MW + PADT + w(k) + GAPI;
        end
    end
    hold(axL,'off');
end
