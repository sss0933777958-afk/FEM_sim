function plot_svd_polar(USE_BIAS, R_um, SRC, MODEL, GEOM)
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
        addpath(fullfile(CAL,'utils','long2016_hexapole_halfcut'));   % sphere_grid_sample
        F = zeros(6, cfg.N_I);
        for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end
        %   ki_gate：K̄_I「非對角全負」這道閘只對 long2016 成立；六極不等強的設計
        %   （zhi_peng：±x 兩根的場是另四根的 2.3 倍）該條恆不成立 → 關掉。
        %   非 long2016：K̄_I 的物理結構條件不適用（六極不等強）→ **完全不納入判準**，
        %   收斂點只由 l_hat 與 g_I 決定（使用者拍板 2026-08-15）。
        is_l2016 = strcmp(MODEL,'long2016_hexapole_halfcut');
        sg = struct('model',MODEL, 'geom',GEOM, ...
                    'ki_gate', is_l2016, 'ki_req', is_l2016);
        [tri, Nc] = conv_design(R_FIT*1e-6, cfg, F, USE_BIAS, here, cfg.R_norm, sg);
        assert(all(tri > 0), '%s：R=%d um 找不到收斂點設計', MODEL, R_FIT);
        og = struct('frame','actuator', 'NRPT',tri, 'model',MODEL, 'geom',GEOM);
        [gx, gy, gz, gB] = sphere_grid_sample(R_FIT*1e-6, [], og);
        Pq = [gx gy gz];
        Bq = zeros(3*size(Pq,1), size(gB,3));
        for j = 1:size(gB,3), Bq(:,j) = reshape(gB(:,:,j).', [], 1); end
        [e_fit, l_fit] = fitting(Pq, Bq, cfg.Pc_base, cfg.R_norm, USE_BIAS);
        [KI, gI, ~, rm] = solve_current(l_fit, e_fit, cfg.Pc_base, Pq, Bq, F);
        S = struct('l_hat',l_fit, 'e',e_fit, 'gI_hat',gI, 'KI_bar',KI, ...
                   'npts',Nc, 'NMAE',rm.NMAE);
        fprintf('  N_c = %d 點，設計 (%d,%d,%d)\n', Nc, tri);
        ell_m = S.l_hat;                       % [m]
        Hhat  = S.gI_hat * S.KI_bar;           % ᴮĤ_I [mT/A]
        Pc    = make_Pc(S.e, cfg.Pc_base);     % 3×6 電荷格（actuator frame、無因次）
        fprintf('資料 Maxwell %s / %s（N_c 降取樣校正 @R=%d um）\n  l_hat=%.1f um  g_I=%.4f mT/A  N_c=%d  NMAE=%.2f%%\n', ...
                MODEL, tag, R_FIT, ell_m*1e6, S.gI_hat, S.npts, S.NMAE);
    end

    % ---- 極座標網格（連續評估，非 FEM 節點）----
    r_um = linspace(0, R_um, 161);         % 徑向
    th   = (0:2:360)*pi/180;               % 方位
    [RR, TH] = meshgrid(r_um, th);

    ex = [1;0;0];   ey = [0;1;0];   ez = [0;0;1];    % actuator 三軸
    Clab = '$\mathbf{\mathcal{C}\;[(mT/A)^{3}]}$';
    Klab = '$\mathbf{\kappa}$';

    % 三個切面：(u, v, w) = 面內兩基底 + 面法向（判哪些極落在面內）
    FACE = { 'xaya', ex, ey, ez, '$\mathbf{x_{a}\!-\!y_{a}\;plane}$';
             'xaza', ex, ez, ey, '$\mathbf{x_{a}\!-\!z_{a}\;plane}$';
             'yaza', ey, ez, ex, '$\mathbf{y_{a}\!-\!z_{a}\;plane}$' };
    for a = 1:size(FACE,1)
        [fname, u, v, w, ttl] = FACE{a,:};
        [C, K] = ck_grid(RR, TH, u, v, ell_m, Hhat, Pc);
        fprintf('  %s：C %.0f ~ %.0f (mT/A)^3   kappa %.3f ~ %.3f\n', ...
                fname, min(C(:)), max(C(:)), min(K(:)), max(K(:)));
        % [FIXED 2026-08-15] 檔名必須帶 msfx，否則不同 model 會互相覆蓋（已踩過：
        %   zhi_peng 的圖蓋掉 long2016 的六張）。long2016 的 msfx='' → 檔名不變。
        render_polar(RR, TH, C, Clab, R_um, Pc, u, v, w, ttl, ...
                     fullfile(figdir, sprintf('gain_polar_%s_%s%s.png', fname, lower(SRC), msfx)));
        render_polar(RR, TH, K, Klab, R_um, Pc, u, v, w, ttl, ...
                     fullfile(figdir, sprintf('iso_polar_%s_%s%s.png', fname, lower(SRC), msfx)));
    end
end

% ============================================================================
function [C, K] = ck_grid(RR, TH, u, v, ell_m, Hhat, Pc)
% 逐格解析評估：p = r(u cosθ + v sinθ) → S(p)·Ĥ 的奇異值 → C=∏σ、κ=σ₃/σ₁
    C = zeros(size(RR));   K = zeros(size(RR));
    for a = 1:numel(RR)
        p  = (RR(a)*1e-6) * (u*cos(TH(a)) + v*sin(TH(a)));   % 3×1 [m]
        Dk = p/ell_m - Pc;                                    % 3×6（無因次）
        sv = svd( (Dk ./ (vecnorm(Dk).^3)) * Hhat );
        C(a) = prod(sv);   K(a) = sv(3)/sv(1);
    end
end

% ============================================================================
function render_polar(RR, TH, val, clab, R, Pc, u, v, w, ttl, out)
    FSA = 26;   % 角度標籤
    FSR = 20;   % 半徑環標籤
    FSP = 26;   % 磁極標籤
    FSC = 32;   % colorbar
    X = RR.*cos(TH);   Y = RR.*sin(TH);

    fig = figure('Color','w','Position',[80 80 1180 980]);
    ax  = axes(fig);   hold(ax,'on');
    surf(ax, X, Y, zeros(size(X)), val, 'EdgeColor','none');
    view(ax,2);   shading(ax,'interp');
    colormap(ax, jet);   clim(ax, [min(val(:)) max(val(:))]);
    axis(ax,'equal');    axis(ax,'off');

    % ---- 極座標網格 overlay（畫在高 z、view(2) 下蓋住填色）----
    zt  = max(val(:)) + 1;
    thg = linspace(0, 2*pi, 360);
    rla = 107*pi/180;                                  % 半徑標籤角度（避開 spoke，且離 90° 標籤遠一點）
    % 半徑環：取整數刻度（figure-style）。R=500 → 125/250/375/500；R=150 → 50/100/150
    if     mod(R,4) == 0, s = R/4;
    elseif mod(R,3) == 0, s = R/3;
    else,                 s = round(R/4);
    end
    rings = s:s:R;
    for rr = rings
        plot3(ax, rr*cos(thg), rr*sin(thg), zt*ones(size(thg)), '-', ...
              'Color',[.35 .35 .35], 'LineWidth',1.2);
        text(ax, rr*cos(rla), rr*sin(rla), zt, sprintf('%d\\mum', rr), ...
             'FontSize',FSR, 'FontWeight','bold', 'Color','k', ...
             'HorizontalAlignment','center', 'BackgroundColor','w', 'Margin',1);
    end
    for a = 0:30:330
        ar = a*pi/180;
        plot3(ax, [0 R*cos(ar)], [0 R*sin(ar)], [zt zt], '-', ...
              'Color',[.35 .35 .35], 'LineWidth',1.0);
        text(ax, 1.28*R*cos(ar), 1.28*R*sin(ar), sprintf('%d\\circ', a), ...
             'HorizontalAlignment','center', 'FontSize',FSA, 'FontWeight','bold');
    end

    % ---- 面內磁極標記（|Pc_k · w| ≈ 0 → 該極落在此切面上）----
    for k = 1:size(Pc,2)
        dk = Pc(:,k) / norm(Pc(:,k));
        if abs(dk.'*w) < 0.15          % 門檻放寬到 0.15：eighteen 的下極離面 0.114（面外極是 ~1.0，仍分得開）
            tk = atan2(dk.'*v, dk.'*u);
            plot3(ax, R*cos(tk), R*sin(tk), zt+1, 'o', 'MarkerSize',13, ...
                  'MarkerFaceColor',[0.90 0.90 0.90], 'MarkerEdgeColor',[0.5 0 0.12], 'LineWidth',2.0);
            % [MODIFIED 2026-08-15] P 標籤移到**圓外**，且沿切線**逆時針錯開 10°**：
            %   六個極正好落在 0/90/180/270°，與角度標籤同方向 —— 不錯開的話兩個文字框
            %   必然沿同一半徑排隊、互相擠壓（放多遠都一樣）。錯開後離最近的角度標籤 20°。
            tlab = tk + 10*pi/180;
            text(ax, 1.13*R*cos(tlab), 1.13*R*sin(tlab), zt+1, sprintf('P%d',k), ...
                 'FontSize',FSP, 'FontWeight','bold', 'Color',[0.5 0 0.12], ...
                 'HorizontalAlignment','center', 'BackgroundColor','w', 'Margin',1, ...
                 'EdgeColor',[0.5 0 0.12]);
        end
    end
    xlim(ax, [-1.42*R 1.42*R]);   ylim(ax, [-1.42*R 1.42*R]);

    cb = colorbar(ax);
    cb.FontSize = FSC;   cb.FontWeight = 'bold';
    cb.Label.Interpreter = 'latex';                    % 先設 Interpreter 再設 String
    cb.Label.String = clab;   cb.Label.FontSize = FSC;
    if numel(cb.Ticks) >= 5, cb.Ticks = cb.Ticks(1:2:end); end
    ax.Toolbar.Visible = 'off';
    hold(ax,'off');

    % [MODIFIED 2026-08-15] 切面標題改成**方框樣式**（使用者指定）：黑粗框 + 白底，
    %   與圖例的外框一致；用 annotation 放在 axes 上方（不佔 axes 座標範圍）。
    drawnow;
    axp = get(ax,'Position');
    annotation(fig, 'textbox', [axp(1)+axp(3)/2-0.16, 0.925, 0.32, 0.062], ...
        'String', ttl, 'Interpreter','latex', 'FontSize',FSC, ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
        'EdgeColor','k', 'LineWidth',2.5, 'BackgroundColor','w', 'FitBoxToText','off');

    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
    close(fig);
end

% ============================================================================
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
