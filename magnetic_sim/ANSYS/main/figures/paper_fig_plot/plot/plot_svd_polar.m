function plot_svd_polar(USE_BIAS, R_um, SRC, MODEL, GEOM, VARIANT)
% plot_svd_polar -- 控制指標 C（gain）與 κ（iso）在三個 actuator 參考切面上的極座標熱圖
% =========================================================================
%   **Maxwell 版**（資料源 = matlab/Maxwell 的 long2016 校正結果）。
%   對照組是 APDL 版 plot/long2016_hexapole_halfcut/current/plot_svd_heatmaps_2d.m
%   （R=500 µm、輸出到 Calibration 的 figures/single_param）。本檔差異：
%     ① 資料改用 Maxwell 校正的 (l_hat, e, g_I, K_I_bar)
%     ② R = 150 µm（校正的定案取樣半徑）
%     ③ 輸出到 figures/paper_fig/Section4_C/，字級與刻度照 figure-style
%
%   物理：C、κ 是**位置相關的解析函數**（點電荷模型導出，非 raw FEM 場）：
%     S(p) = (p̄ − Pc) ./ |p̄ − Pc|³          (3×6，p̄ = p/l_hat 無因次)
%     sv   = svd( S(p) · ᴮĤ_I )              ᴮĤ_I = ĝ_I · K̄_I  [mT/A]
%     C    = ∏ sv  [(mT/A)³]                 κ = sv₃/sv₁  [無因次]
%   故可在極座標網格上逐點連續評估（不受 FEM 節點位置限制）。
%
%   座標：Maxwell 的校正**本來就在 actuator frame**（Pc_base 是整數 canonical
%   ±ex/±ey/±ez），所以三個參考切面就是該框的 xy / yz / xz 平面：
%     x_a-y_a：面內極 P1(0°) P3(90°) P2(180°) P4(270°)
%     x_a-z_a：面內極 P1(0°) P5(90°) P2(180°) P6(270°)
%     y_a-z_a：面內極 P3(0°) P5(90°) P4(180°) P6(270°)
%   （APDL 版因為在 measure frame 算，才需要 ea_x/ea_y/ea_z = dhat(:,[1 3 5]) 轉換。）
%
%   風格：極座標圖無傳統軸 → 不放軸標題；角度 / 半徑標籤與 colorbar 皆粗體大字，
%   半徑環取**整數** 50/100/150 µm（figure-style「刻度用整數」）。
%
%   用法：plot_svd_polar            → eighteen、R=150（預設）
%         plot_svd_polar(false)     → single
%   輸出 → figures/paper_fig/Section4_C/{gain,iso}_polar_{xaya,xaza,yaza}_maxwell.png
% =========================================================================
    clc;
    % [預設 single] 對應參考圖（APDL 版 plot_svd_heatmaps_2d 用 fit_fixl = single）。
    %   ⚠ eighteen 的下極電荷偏移大（P1 的 Pc = [1.018, -0.118, +0.118]，離面 0.114），
    %     會破壞切面的四重對稱、且面內極判定門檻要放寬才標得出來。
    if nargin < 1 || isempty(USE_BIAS), USE_BIAS = false; end
    if nargin < 2 || isempty(R_um),     R_um     = 500;   end   % 使用者拍板：畫到 500 µm
    if nargin < 3 || isempty(SRC),      SRC      = 'maxwell'; end   % 'maxwell' | 'apdl'（對照）
    % [ADDED 2026-08-15] 參數化 model：long2016 之外也能用（如 zhi_peng 平面六極）。
    %   非 long2016 時輸出檔名加 _<model>，不覆蓋既有 long2016 的圖。
    if nargin < 4 || isempty(MODEL),    MODEL = 'long2016_hexapole_halfcut'; end
    if nargin < 5,                      GEOM  = 'tip40um';               end
    if strcmp(MODEL,'long2016_hexapole_halfcut') && isempty(GEOM), GEOM = 'tip40um'; end
    assert(~(strcmpi(SRC,'apdl') && ~strcmp(MODEL,'long2016_hexapole_halfcut')), ...
           'SRC=''apdl'' 的對照組只有 long2016 有');
    msfx = ''; if ~strcmp(MODEL,'long2016_hexapole_halfcut'), msfx = ['_' MODEL]; end

    here   = fileparts(fileparts(mfilename('fullpath')));       % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section4_C');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    addpath(fullfile(CAL,'function'), fullfile(CAL,'utils'), fullfile(CAL,'common_path'));
    cfg = model_config(MODEL, GEOM);
    % [ADDED 2026-08-21] variant 參數化。'' → cfg.default_variant（既有呼叫端行為不變）。
    %   同一個 model 有多版場時（zhi_peng：maxwell / maxwell_split）**必須明給**，否則
    %   conv_design_ws 會退回 default_variant、**靜默**用到舊場。
    if nargin < 6 || isempty(VARIANT), VARIANT = cfg.default_variant; end
    % [ADDED 2026-08-24] 檔名補上 **variant** 後綴：同一個 model 有多版場（zhi_peng 的
    %   maxwell / maxwell_split / maxwell_gap）時，只帶 model 會**互相覆蓋**。
    %   非預設 variant 才加，且照 short-names 規則剝掉整棵樹都一樣的分支名 'maxwell'
    %   -> 'maxwell_gap' 變 '_gap'、'maxwell_split' 變 '_split'。
    if ~strcmpi(VARIANT, cfg.default_variant)
        msfx = [msfx '_' regexprep(VARIANT, '^maxwell_?', '')];
    end

    tag = 'single';  if USE_BIAS, tag = 'eighteen'; end
    if strcmpi(SRC, 'apdl')
        % [ADDED 2026-08-15] 對照組：改餵 **APDL** 的校正結果（就是參考圖 gain_polar_*.png
        %   用的那顆 fit_fixl_R150um_gap_200um.mat）。用來證明「同一支腳本、同一個 R，
        %   換資料源會長成什麼樣」——把「畫錯」與「資料/半徑差異」分離。
        CALA = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
        A = load(fullfile(CALA,'data','long2016_hexapole_halfcut','.mat','fit_fixl_R150um_gap_200um.mat'));
        ell_m = A.ell*1e-6;   Hhat = A.gB * A.Khat;   Pc = cfg.Pc_base;   % fit_fixl = single、無偏移
        fprintf('資料 APDL fit_fixl (參考圖用的同一顆)\n  l_hat=%.1f um  gB=%.4f mT/A\n', ell_m*1e6, A.gB);
    else
        % ⚠ **校正半徑與繪圖半徑是兩回事**：校正用定案的 R=150 µm，圖畫到 R_um（500）。
        %   參考圖也是這樣（fit_fixl_R150um 的結果畫到 500）。
        % [MODIFIED 2026-08-15] 校正改用 **N_c 降取樣**（使用者拍板）：不再讀 main.m 產的
        %   全格點 calib .mat（那是 R<=150 內全部 1771 個 .fld 格點），改成就地找該 R 的
        %   收斂點設計，用那 N_c 個等測度網格點擬合。
        R_FIT = 150;
        % [MODIFIED 2026-08-23] 取樣+取場已併入 conv_design_ws（上方已 addpath）
        F = zeros(6, cfg.N_I);
        for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end
        %   ki_gate：K̄_I「非對角全負」這道閘只對 long2016 成立；六極不等強的設計
        %   （zhi_peng：±x 兩根的場是另四根的 2.3 倍）該條恆不成立 → 關掉。
        %   非 long2016：K̄_I 的物理結構條件不適用（六極不等強）→ **完全不納入判準**，
        %   收斂點只由 l_hat 與 g_I 決定（使用者拍板 2026-08-15）。
        is_l2016 = strcmp(MODEL,'long2016_hexapole_halfcut');
        sg = struct('model',MODEL, 'geom',GEOM, 'variant',VARIANT, ...
                    'ki_gate', is_l2016, 'ki_req', is_l2016);
        % ---- 讀 main.m 產的收斂設計校正結果（不再自己重跑階梯 + 校正）--------
        %   [MODIFIED 2026-08-23 使用者拍板] 校正與收斂判準已搬回 main.m
        %   （conv_design_ws / conv_design_sensor 只負責決定內插點位置與取場），
        %   繪圖端改成**接收 main 產完的結果** -> 圖與結果 PDF 保證出自同一次校正。
        %   ⚠ 該組合必須先跑過 main.m（GRID_NRPT='auto'）；找不到就報錯，不猜。
        %   ⚠ 同一組合可能有多顆 convN 檔（舊實驗留下的，例如 long2016 R150
        %     eighteen 就有 convN6/convN80/convN88）-> 只認 conv_auto==true 那顆。
        md_ = fullfile(CAL, 'data', MODEL, '.mat');
        tg_ = 'single';   if USE_BIAS, tg_ = 'eighteen'; end
        dd_ = dir(fullfile(md_, sprintf('calib_current_%s_convN*_R%03d_%s.mat', ...
                                        VARIANT, round(R_FIT), tg_)));
        if numel(dd_) > 1
            ok_ = false(1, numel(dd_));
            for k_ = 1:numel(dd_)
                f_ = fullfile(md_, dd_(k_).name);   w_ = whos('-file', f_);
                if ismember('conv_auto', {w_.name})
                    r_ = load(f_, 'conv_auto');   ok_(k_) = logical(r_.conv_auto);
                end
            end
            dd_ = dd_(ok_);
        end
        assert(numel(dd_) == 1, ['找到 %d 顆收斂校正 .mat（需恰好 1 顆）。請先跑 ' ...
               'main.m：MODEL=''%s''、R_select=%ge-6、USE_BIAS=%d、GRID_NRPT=''auto''。'], ...
               numel(dd_), MODEL, R_FIT, USE_BIAS);
        cal = load(fullfile(md_, dd_(1).name));
        S = struct('l_hat',cal.l_hat, 'e',cal.e, 'gI_hat',cal.gI_hat, ...
                   'KI_bar',cal.KI_bar, 'npts',cal.npts, 'NMAE',cal.NMAE);
        tri = cal.GRID_NRPT;   Nc = cal.npts;
        fprintf('  N_c = %d 點，設計 (%d,%d,%d)\n', Nc, tri);
        ell_m = S.l_hat;                       % [m]
        Hhat  = S.gI_hat * S.KI_bar;           % ᴮĤ_I [mT/A]
        Pc    = make_Pc(S.e, cfg.Pc_base);     % 3×6 電荷格（actuator frame、無因次）
        fprintf('資料 Maxwell %s / %s / %s（N_c 降取樣校正 @R=%d um）\n  l_hat=%.1f um  g_I=%.4f mT/A  N_c=%d  NMAE=%.2f%%\n', ...
                MODEL, VARIANT, tag, R_FIT, ell_m*1e6, S.gI_hat, S.npts, S.NMAE);
    end

    % ---- 極座標網格（連續評估，非 FEM 節點）----
    r_um = linspace(0, R_um, 161);         % 徑向
    th   = (0:2:360)*pi/180;               % 方位
    [RR, TH] = meshgrid(r_um, th);

    ex = [1;0;0];   ey = [0;1;0];   ez = [0;0;1];    % actuator 三軸
    % [MODIFIED 2026-08-18 使用者拍板] 標籤拿掉 ^{1/3}（數值仍是三個奇異值的幾何平均、單位 mT/A）
    Clab = '$\mathbf{\mathcal{C}\;[mT/A]}$';
    Klab = '$\mathbf{\kappa}$';

    % 三個切面：(u, v, w) = 面內兩基底 + 面法向（判哪些極落在面內）
    %   （2026-08-18：曾加「面內兩軸名稱」兩欄，使用者拍板只留軸線、不標名 → 已收回）
    FACE = { 'xaya', ex, ey, ez, '$\mathbf{x_{a}\!-\!y_{a}\;plane}$';
             'xaza', ex, ez, ey, '$\mathbf{x_{a}\!-\!z_{a}\;plane}$';
             'yaza', ey, ez, ex, '$\mathbf{y_{a}\!-\!z_{a}\;plane}$' };
    % [MODIFIED 2026-08-18 使用者拍板] 版面重排：
    %   ① x_a-y_a 與 x_a-z_a **左右並排成一張圖**、共用一支 colorbar（放最右邊）
    %   ② y_a-z_a **獨立一張**
    %   ③ 但**三個切面共用同一個色階**（clim 取三面全域 min/max）→ 三張圖可直接互比
    %   故必須先把三個切面都算完、定出全域 clim，才能開始畫。
    NF = size(FACE,1);
    Cv = cell(NF,1);   Kv = cell(NF,1);
    for a = 1:NF
        [fname, u, v, ~, ~] = FACE{a,:};
        [Cv{a}, Kv{a}] = ck_grid(RR, TH, u, v, ell_m, Hhat, Pc);
        fprintf('  %s：C^(1/3) %.3f ~ %.3f mT/A   kappa %.3f ~ %.3f\n', ...
                fname, min(Cv{a}(:)), max(Cv{a}(:)), min(Kv{a}(:)), max(Kv{a}(:)));
    end
    clC = [min(cellfun(@(x)min(x(:)),Cv)), max(cellfun(@(x)max(x(:)),Cv))];
    clK = [min(cellfun(@(x)min(x(:)),Kv)), max(cellfun(@(x)max(x(:)),Kv))];
    fprintf('  共用色階：C^(1/3) %.3f ~ %.3f mT/A｜kappa %.3f ~ %.3f\n', clC, clK);


    % [FIXED 2026-08-15] 檔名必須帶 msfx，否則不同 model 會互相覆蓋。
    P12 = [1 2];   P3 = 3;                    % 1=xaya 2=xaza（並排）／3=yaza（獨立）
    render_polar(RR, TH, Cv, Clab, R_um, Pc, FACE, P12, clC, @jet, ...
        fullfile(figdir, sprintf('gain_polar_xaya_xaza_%s%s.png', lower(SRC), msfx)));
    render_polar(RR, TH, Cv, Clab, R_um, Pc, FACE, P3,  clC, @jet, ...
        fullfile(figdir, sprintf('gain_polar_yaza_%s%s.png', lower(SRC), msfx)));
    render_polar(RR, TH, Kv, Klab, R_um, Pc, FACE, P12, clK, @jet, ...
        fullfile(figdir, sprintf('iso_polar_xaya_xaza_%s%s.png', lower(SRC), msfx)));
    render_polar(RR, TH, Kv, Klab, R_um, Pc, FACE, P3,  clK, @jet, ...
        fullfile(figdir, sprintf('iso_polar_yaza_%s%s.png', lower(SRC), msfx)));
end

% ============================================================================
function [C, K] = ck_grid(RR, TH, u, v, ell_m, Hhat, Pc)
% 逐格解析評估：p = r(u cosθ + v sinθ) → S(p)·Ĥ 的奇異值 → C=∏σ、κ=σ₃/σ₁
    C = zeros(size(RR));   K = zeros(size(RR));
    for a = 1:numel(RR)
        p  = (RR(a)*1e-6) * (u*cos(TH(a)) + v*sin(TH(a)));   % 3×1 [m]
        Dk = p/ell_m - Pc;                                    % 3×6（無因次）
        sv = svd( (Dk ./ (vecnorm(Dk).^3)) * Hhat );
        % [MODIFIED 2026-08-18 使用者拍板] C 改回報**三個奇異值相乘後開三次方根**
        %   （幾何平均）→ 單位由 (mT/A)^3 變成 mT/A，量級與 g_I 同級、好解讀。
        C(a) = prod(sv)^(1/3);   K(a) = sv(3)/sv(1);
    end
end

% ============================================================================
function render_polar(RR, TH, vals, clab, R, Pc, FACE, idx, cl, cmapf, out)
% 多面板極座標熱圖。idx = 要畫的切面索引（[1 2] → 左右並排；純量 → 單張）。
%   cl = **共用色階** [lo hi]（由呼叫端跨三個切面算出），所有面板與 colorbar 都用它。
%   cmapf = colormap 函式 handle。gain / iso 目前**都用 @jet**
%     （2026-08-18 曾把 iso 改 @turbo 增強低值區對比，使用者拍板改回原色）。
%   colorbar 只畫一支、放**最右邊**。
% [MODIFIED 2026-08-18 使用者拍板]
%   ③ 半徑環標籤**不加白底**、黑色粗體
%   ④ P1..P6 標籤：一度移除後**又加回**（使用者要求）
    % [MODIFIED 2026-08-18] 字級全部統一 **36**，照 figure-style「字體大小統一 = 36
    %   （所有 paper 圖通用）：刻度數字、軸標題/colorbar 標題同 36」。
    %   極座標圖的「刻度數字」= 角度標籤 + 半徑環標籤，故一併 36。
    FSA = 36;   % 角度標籤（= 方位刻度數字）
    % [MODIFIED 2026-08-18 使用者拍板] 半徑環標籤 36 -> **30**：36 級時四個標籤沿
    %   120° 方位排下來會撐出圓盤、並與角度標籤打架（其餘字級維持規則的 36）。
    FSR = 30;   % 半徑環標籤（= 徑向刻度數字）
    FSP = 36;   % 磁極標籤
    FSC = 36;   % colorbar 刻度數字與標題
    % [ADDED 2026-08-20 使用者拍板] 顏色軸：**起點與終點都要有數值**、刻度等距、
    %   數字為整數或一位小數 0.N。做法：把 clim **往外擴到乾淨邊界**再均分。
    [cl2, ctk] = cbar_ticks(cl);
    np  = numel(idx);
    X = RR.*cos(TH);   Y = RR.*sin(TH);

    % [MODIFIED 2026-08-18 使用者拍板] ① 拿掉上方的切面標題方框 → 座標軸往上撐滿；
    %   ② colorbar 依 figure-style「Colorbar 寬度（paper 合併圖）」：寬度 = 佔圖寬的
    %      **固定比例 CBW_RATIO = 0.009**（細長條）。normalized 單位本身就是「佔圖寬比例」，
    %      故直接把寬度設成 CBW_RATIO —— 兩張圖（寬 2025 / 1180 px）縮到同寬時視覺同粗。
    %   ③ 並往內收緊：colorbar 左緣只離右側面板 0.010，不再留大片空白。
    CBW_RATIO = 0.009;
    if np == 1
        fig = figure('Color','w','Position',[80 80 1180 980]);
        % [MODIFIED 2026-08-18 使用者拍板] 單張與並排**圓盤大小/圖高要一致**：plot box 已鎖正方形，
        %   邊長 = min(軸寬, 軸高)。兩個版面的軸高都是 0.940*980 = 921 px，只要把軸寬都設成
        %   **大於 921 px**，邊長就同由高度決定 → 圓一樣大、匯出圖高一樣。
        %   （並排 0.440*2060 = 906 px、單張 0.780*1180 = 920 px，兩個都偏小且不相等。）
        AXP = {[0.0150 0.030 0.790 0.940]};
        CBP = [0.815 0.105 CBW_RATIO 0.790];
    else
        % [MODIFIED 2026-08-18 使用者拍板] 中間收緊：兩面板中心距 912 -> 817 px（各往中間 ~48 px），
        %   圖寬同步 2120 -> 2025、colorbar 左移，右側留白與單張版一致。
        %   收緊上限 = 左panel 的 0 度標籤與右panel 的 180 度標籤不可相碰；
        %   實測空白由 249 px 縮到 ~60 px（PNG 寬 4015 的 1.5%），仍未接觸。
        fig = figure('Color','w','Position',[80 80 2060 980]);
        AXP = {[0.0033 0.030 0.450 0.940], [0.4305 0.030 0.450 0.940]};
        CBP = [0.8772 0.105 CBW_RATIO 0.790];
    end

    for q = 1:np
        [~, u, v, w, ~] = FACE{idx(q),:};
        val = vals{idx(q)};
        ax  = axes(fig, 'Position', AXP{q});   hold(ax,'on');
        surf(ax, X, Y, zeros(size(X)), val, 'EdgeColor','none');
        view(ax,2);   shading(ax,'interp');
        colormap(ax, cmapf(256));   clim(ax, cl2);  % 三個切面共用色階（已擴到乾淨邊界）
        axis(ax,'off');

        % ---- 極座標網格 overlay（畫在高 z、view(2) 下蓋住填色）----
        zt  = cl(2) + 1;
        thg = linspace(0, 2*pi, 360);
        % [MODIFIED 2026-08-18 使用者拍板] 半徑標籤方位：107° → **120°**。
        %   ⚠ 120° 正好是 spoke（每 30° 一根）的方位，標籤會壓在那條線上。
        %   （曾短暫改成依底色亮度自動選角，已作廢。）
        rla = 120*pi/180;
        if     mod(R,4) == 0, sr = R/4;
        elseif mod(R,3) == 0, sr = R/3;
        else,                 sr = round(R/4);
        end
        for rr = sr:sr:R
            plot3(ax, rr*cos(thg), rr*sin(thg), zt*ones(size(thg)), '-', ...
                  'Color',[.35 .35 .35], 'LineWidth',1.2);
            % [MODIFIED 2026-08-18 使用者拍板] **半徑數字拿掉**（只留環線）——所有極座標圖通用。
        end
        for aa = 0:30:330
            ar = aa*pi/180;
            plot3(ax, [0 R*cos(ar)], [0 R*sin(ar)], [zt zt], '-', ...
                  'Color',[.35 .35 .35], 'LineWidth',1.0);
            % [MODIFIED 2026-08-18 使用者拍板] **角度數字拿掉**（只留 spoke 線）——所有極座標圖通用。
            %   （同日稍早的「只標 0/90/180/270」已被此次取代。）
        end

        % ---- 面內座標軸：穿過圓心的兩根軸（使用者拍板 2026-08-18）----
        %   由圓心指向 0 度（橫軸）與 90 度（縱軸）的黑色箭頭。**只有軸線、不標軸名**
        %   （軸名 x_a/y_a/z_a 一度加過，使用者拍板拿掉）；長度 0.90R -> **0.45R**。
        LAX = 0.45*R;
        draw_axis(ax, [1 0], LAX, zt+2);        % 橫軸 -> 0 度
        draw_axis(ax, [0 1], LAX, zt+2);        % 縱軸 -> 90 度

        %   [MODIFIED 2026-08-18 使用者拍板] 標籤整體外移：角度 1.28R→1.40R、
        %   P 標籤 1.13R→1.20R（兩者間距由 0.15R 拉開到 0.20R，下方 P4/P6 與 270° 不再貼合）。
        for k = 1:size(Pc,2)
            dk = Pc(:,k) / norm(Pc(:,k));
            if abs(dk.'*w) < 0.15
                tk = atan2(dk.'*v, dk.'*u);
                plot3(ax, R*cos(tk), R*sin(tk), zt+1, 'o', 'MarkerSize',13, ...
                      'MarkerFaceColor',[0.90 0.90 0.90], ...
                      'MarkerEdgeColor',[0.5 0 0.12], 'LineWidth',2.0);
                % [MODIFIED 2026-08-18] P 標籤加回（使用者要求）。放**圓外** 1.13R，
                %   並沿切線逆時針錯開 10° —— 六個極正好落在 0/90/180/270°，與角度
                %   標籤同方向，不錯開的話兩個文字框會沿同一半徑排隊互相擠壓。
                % [MODIFIED 2026-08-18] 角度數字沒了，10 度錯開（原為避開同半徑的角度標籤）也不需要
                %   → 錯開歸零、P 標籤正對磁極方位；半徑 1.20R -> 1.15R。
                tlab = tk;
                text(ax, 1.15*R*cos(tlab), 1.15*R*sin(tlab), zt+1, sprintf('P%d',k), ...
                     'FontSize',FSP, 'FontWeight','bold', 'Color',[0.5 0 0.12], ...
                     'HorizontalAlignment','center', 'BackgroundColor','w', 'Margin',1, ...
                     'EdgeColor',[0.5 0 0.12]);
            end
        end
        % [MODIFIED 2026-08-18 使用者拍板] 角度 / 半徑數字拿掉後，圓外只剩 P 標籤（1.15R）
        %   → 外緣 1.75R 收到 **1.35R**：圓盤放大、版面緊湊。
        xlim(ax, [-1.35*R 1.35*R]);   ylim(ax, [-1.35*R 1.35*R]);
        % [FIXED 2026-08-18] 圓被壓成橢圓的根因：`axis equal` 是在 surf 之後、格線 plot3
        %   （帶 z = zt）之前呼叫的，之後 ZLim 被撐開 → 當時定下的 plot box 比例失效，
        %   x/y 就各自被拉去填滿 Position（906 x 921 px → 橫向壓 1.6%）。
        %   正解：x/y 範圍本來就相等 → 直接把 plot box 鎖成正方形（z 是虛擬軸，view(2) 不管它）。
        set(ax, 'DataAspectRatioMode','auto', 'PlotBoxAspectRatio',[1 1 1]);
        ax.Toolbar.Visible = 'off';
        hold(ax,'off');

        % [MODIFIED 2026-08-18] 切面標題方框已移除（使用者拍板）。
        if q == np                                  % colorbar 只掛一支、放最右
            cb = colorbar(ax);
            cb.FontSize = FSC;   cb.FontWeight = 'bold';
            cb.Label.Interpreter = 'latex';          % 先設 Interpreter 再設 String
            cb.Label.String = clab;   cb.Label.FontSize = FSC;
            % [MODIFIED 2026-08-20] 不再「減半」（那會把上緣丟掉）：直接給均分刻度，
            %   首尾就是 clim 兩端 -> 起點與終點都有數字。
            cb.Limits = cl2;   cb.Ticks = ctk;
            set(ax, 'Position', AXP{q});             % colorbar 會擠壓 axes → 還原
            cb.Position = CBP;
        end
    end

    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
    close(fig);
end

% ============================================================================
% ============================================================================
function [cl2, tk] = cbar_ticks(cl)
% [ADDED 2026-08-20 使用者拍板] colorbar 刻度：
%   • **起點與終點都要有數值**（首尾刻度就坐在 clim 兩端）
%   • 刻度**等距**
%   • 每個數字是**整數或一位小數 0.N**
% 做法：選一個步長 s（一位小數以內的 nice 值），把 clim **往外**擴到 s 的倍數，
%   區間數 n 限制在 3~6。往外擴對顏色映射的影響是均勻拉伸，
%   三個切面共用同一個 cl2，互比性不受影響。
%
% [MODIFIED 2026-08-24] 挑選判準由「n 最接近 4」改成「**跨度最小**（= 資料覆蓋最大）」，
%   同分再取整數刻度。起因：舊判準只管刻度好不好看、不管浪費多少色階 ——
%   zhi_peng gap 的資料 19.90~30.83 被撐成 clim [15,35]，**只用到色條的 54.6%**，
%   而同一支腳本畫的無氣隙版（31.86~51.10 -> [30,55]）用了 77.0%。兩張圖的色階
%   利用率差 22 個百分點，看起來就像「分布不一樣」（實測空間形狀只差 3%）。
%   新判準給 gap [18,33]（覆蓋 72.8%、中心落在 34.2% vs 無氣隙的 32.2%）。
%   ⚠ 已逐一驗算：**既有的其他圖 clim 完全不變** —— no gap C [30,55]、
%     long2016 C [13,18]、long2016 kappa [0,1]、zhi_peng kappa [0,0.8] 都是原值。
%   ⚠ cand 補上 3（figure-style 的 nice 清單本來就有），[18,33] 才生得出來。
    cand = [0.1 0.2 0.3 0.5 1 2 3 5 10 20 30 50 100 200 500];
    bs = [];   bsp = inf;   bint = -1;
    for s = cand
        lo = floor(cl(1)/s + 1e-9)*s;
        hi = ceil( cl(2)/s - 1e-9)*s;
        n  = round((hi-lo)/s);
        if n < 3 || n > 6, continue; end
        sp   = hi - lo;                                   % 跨度：越小 = 資料覆蓋越大
        tkq  = lo:s:hi;
        isI  = double(all(abs(tkq - round(tkq)) < 1e-9)); % 1 = 刻度全整數（同分時優先）
        if sp < bsp - 1e-9 || (abs(sp - bsp) < 1e-9 && isI > bint)
            bsp = sp;   bint = isI;   bs = [lo hi s];
        end
    end
    if isempty(bs)                                        % 保險：直接四等分
        bs = [cl(1) cl(2) (cl(2)-cl(1))/4];
    end
    lo = bs(1);   hi = bs(2);   s = bs(3);
    tk  = round((lo:s:hi)/0.1)*0.1;                       % 抹掉浮點尾巴
    cl2 = [tk(1) tk(end)];                                % ← 用 tk 定 clim，首尾才不會差一個 ulp 被丟掉
end

% ============================================================================
function draw_axis(ax, d, Lx, z)
% 由圓心指向 d（面內單位向量）的黑色箭頭軸：桿 + 兩翼箭頭
%   （自己畫、不用 quiver：quiver 的頭大小綁資料尺度，換 R 就變形）
    d = d(:).'/norm(d);   n = [-d(2) d(1)];        % 面內法向（用來畫兩翼）
    tip = Lx*d;
    plot3(ax, [0 tip(1)], [0 tip(2)], [z z], '-', 'Color','k', 'LineWidth',3.0);
    hl = 0.15*Lx;   hw = 0.06*Lx;                  % 箭頭長 / 半寬（綁軸長，縮軸時等比縮）
    b  = tip - hl*d;
    plot3(ax, [tip(1) b(1)+hw*n(1)], [tip(2) b(2)+hw*n(2)], [z z], '-', 'Color','k', 'LineWidth',3.0);
    plot3(ax, [tip(1) b(1)-hw*n(1)], [tip(2) b(2)-hw*n(2)], [z z], '-', 'Color','k', 'LineWidth',3.0);
end
function Pc = make_Pc(e17, Pc_base)
% 電荷格 Pc = Pc_base + E(e)（含 e6z 約束；與 solve_current 內的版本一致）
    if isempty(e17) || all(e17(:) == 0), Pc = Pc_base;  return; end
    E = zeros(3, 6);
    E(:,1) = e17(1:3);     E(:,2) = e17(4:6);
    E(:,3) = e17(7:9);     E(:,4) = e17(10:12);
    E(:,5) = e17(13:15);
    E(1,6) = e17(16);      E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);
    Pc = Pc_base + E;
end
