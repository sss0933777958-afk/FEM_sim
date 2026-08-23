function plot_sphere_grid_3d(R_um, SPEC, MONO, force)
% plot_sphere_grid_3d -- 球內「等測度網格」取樣點的 3D 分布（conv_design_ws 的配比）
% =========================================================================
%   直接呼叫 matlab/Maxwell/function/conv_design_ws.m
%   產點（不另抄一份配比公式），把該球內的查詢點畫成 3D 散點。
%
%   用法
%     plot_sphere_grid_3d(150, 2)              過取樣倍率 c = 2 → 內建配比 (5,15,47)
%     plot_sphere_grid_3d(150, [1 2 3], true)  **直接指定三元組** + 單色（全藍）
%
%   輸入
%     R_um  取樣球半徑 [um]（預設 150）
%     SPEC  純量 → 過取樣倍率 c；1x3 向量 → 直接給 [N_r N_phi N_theta]
%     MONO  true = 所有點同一個藍色；false（預設）= 依所屬殼上色（turbo）
%     force true = 重新產點（否則讀快取）
%
%   配比（見 conv_design_ws 檔頭）：N_r : N_phi : N_theta = 1 : 3 : 3pi
%     r     = R * u^(1/3)        u_k = (k-0.5)/N_r        等分 r^3   （等體積殼）
%     phi   = acos(1 - 2*w)      w_j = (j-0.5)/N_phi      等分 cos phi（等面積帶）
%     theta = 2*pi*v             v_i = (i-0.5)/N_theta    等分 theta
%
%   ⚠ 座標框 = measure frame（球心為原點、極軸 = 框的 z）。網格是**繞全域 z 建**的，
%     轉進 actuator frame 會讓極軸相對於框傾斜、配比結構看不出來 —— 此圖為幾何示意，
%     故不轉。要 actuator 版把 FRAME 改成 'actuator' 即可。
%
%   風格①粗體框圖 + 3D 變體 A（三軸同尺度）：daspect([1 1 1])、box off + 手動框邊
%   （省最近角 3 邊）、框線 4.0（ruler 同值）、grid off、字級 36 粗體、刻度數字轉正、
%   三軸同刻度且奇數個不含端點。純幾何示意圖 -> 依規則不放軸標題（刻度數字照留）。
%
%   輸出 → figures/paper_fig/Section2_E/sphere_grid_3d_R<R>_<tag>.png（覆蓋迭代）
% =========================================================================
    clc;
    if nargin < 1 || isempty(R_um), R_um  = 150;   end
    if nargin < 2 || isempty(SPEC), SPEC  = 2;     end
    if nargin < 3 || isempty(MONO), MONO  = false; end
    if nargin < 4 || isempty(force),force = false; end

    FRAME = 'measure';                                     % 'measure'（預設）| 'actuator'
    FS    = 36;                                            % 字級（paper 圖統一）
    LWBOX = 4.0;                                           % 框線寬（ruler 與手動框邊同值）
    AZ    = -30;   EL = 20;                                % 視角
    CBLUE = [0.05 0.10 0.95];   CRED = [0.85 0.10 0.10];   % 與收斂曲線同一組色（single 藍 / eighteen 紅）

    % MONO 解讀：false/[] = 依殼上色；true = 藍；'r'/'red' = 紅；1x3 = 自訂 RGB
    if ischar(MONO) || isstring(MONO)
        switch lower(char(MONO))
            case {'r','red'},  CMONO = CRED;
            case {'b','blue'}, CMONO = CBLUE;
            otherwise, error('MONO 顏色只支援 ''r''/''red''、''b''/''blue'' 或 1x3 RGB');
        end
        MONO = true;
    elseif isnumeric(MONO) && numel(MONO) == 3
        CMONO = MONO(:).';   MONO = true;
    else
        CMONO = CBLUE;       MONO = logical(MONO);
    end

    if isscalar(SPEC), tag = sprintf('c%g', SPEC);
    else,              tag = sprintf('%dx%dx%d', SPEC(1), SPEC(2), SPEC(3));   end

    here   = fileparts(fileparts(mfilename('fullpath')));           % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    % [MODIFIED 2026-08-21] 舊路徑 utils/long2016_hexapole_halfcut 已不存在；
    %   conv_design_ws 在 Maxwell/function/。
    CALMW = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    addpath(fullfile(CALMW,'function'), fullfile(CALMW,'utils'), fullfile(CALMW,'common_path'));

    %% ---- 產點（呼叫真正的取樣器，確保圖與程式一致）------------------------
    cachef = fullfile(here, 'data', sprintf('sphere_grid_R%d_%s.mat', R_um, tag));
    if exist(cachef,'file') && ~force
        S = load(cachef);   P = S.P;   info = S.info;
        fprintf('由快取載入 %s（%d 點）\n', cachef, size(P,1));
    else
        % [MODIFIED 2026-08-23] conv_design_ws 已併入 conv_design_ws：
        %   純量 SPEC = 過取樣倍率 c（由配比反解三元組）；1x3 = 直接給三元組。
        if isscalar(SPEC)
            [Pm, ~, ~, info] = conv_design_ws([], [], [], R_um*1e-6, ...
                                   struct('frame',FRAME, 'c',SPEC, 'quiet',false));
        else
            [Pm, ~, ~, info] = conv_design_ws(SPEC(1), SPEC(2), SPEC(3), R_um*1e-6, ...
                                   struct('frame',FRAME, 'quiet',false));
        end
        P = Pm * 1e6;                                      % m → um
        save(cachef, 'P', 'info');
        fprintf('已存 %s（%d 點）\n', cachef, size(P,1));
    end

    r_k = info.r_k * 1e6;                                  % 殼半徑 [um]
    rr  = vecnorm(P, 2, 2);
    [~, shell] = min(abs(rr - r_k(:).'), [], 2);           % 每點屬於哪一殼（由半徑判，不靠索引序）
    ns = accumarray(shell, 1, [numel(r_k) 1]);
    fprintf('(N_r,N_phi,N_theta) = (%d,%d,%d)｜%d 點｜h = %.2f um\n', ...
            info.N_r, info.N_phi, info.N_theta, size(P,1), info.h*1e6);
    fprintf('r_k [um] = %s｜每殼點數 = %s\n', num2str(r_k,'%.1f '), num2str(ns.','%d '));

    %% ---- figure ------------------------------------------------------------
    fig = figure('Color','w','Position',[80 80 1040 980]);
    ax  = axes(fig);   hold(ax,'on');

    % 取樣球邊界（淡線框，給尺度參考）
    [sx, sy, sz] = sphere(30);
    surf(ax, R_um*sx, R_um*sy, R_um*sz, 'FaceColor','none', ...
         'EdgeColor',[0.55 0.60 0.68], 'EdgeAlpha',0.20, 'LineWidth',0.4);

    % 點大小隨總點數縮放：3525 點時 26，點少時放大到看得見（上限 220）
    MS = min(220, 26 * (3525/max(size(P,1),1))^(1/3));

    % [ADDED] 畫「原點→取樣點」的輻條：3D 散點在透視投影下很難判讀深度，
    %   連到原點的線段一畫出來，方位與仰角立刻讀得出來。
    % [MODIFIED 2026-08-13] 原本只在 N<=100 畫（使用者反映 528 點那張少了連線）→ **一律畫**；
    %   改由點數決定線寬與透明度，點多時自動變細變淡，才不會糊成一團。
    np_ = size(P,1);
    LWs = max(0.4, 2.0 * (100/max(np_,1))^(1/3));      % 100 點→2.0、528 點→1.14、3525→0.63
    ALs = max(0.12, 0.45 * (100/max(np_,1))^(1/3));
    for q = 1:np_
        plot3(ax, [0 P(q,1)], [0 P(q,2)], [0 P(q,3)], '-', ...
              'Color',[CMONO ALs], 'LineWidth',LWs);
    end

    if MONO                                                % 單色（點與輻條同色）
        scatter3(ax, P(:,1), P(:,2), P(:,3), MS, 'filled', ...
                 'MarkerFaceColor', CMONO, 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.90);
    else                                                   % 依殼上色
        cmap = turbo(numel(r_k));
        for k = 1:numel(r_k)
            m = shell == k;
            scatter3(ax, P(m,1), P(m,2), P(m,3), MS, 'filled', ...
                     'MarkerFaceColor', cmap(k,:), 'MarkerEdgeColor','none', ...
                     'MarkerFaceAlpha', 0.90);
        end
    end

    % 原點（規則：黑色標記、不加文字）
    %   用**十字**而非實心圓點：圓點在點雲圖裡會被誤讀成一顆取樣點，
    %   但等測度網格刻意不取 r=0（最內殼在 R*0.1^(1/3)），中心必為空。
    plot3(ax, 0, 0, 0, 'k+', 'LineWidth',3.5, 'MarkerSize',18);

    %% ---- 框 / 刻度（變體 A：同尺度立方）------------------------------------
    bh = 1.2 * R_um;                                       % 框半寬（留一格裕度）
    XL = [-bh bh];   YL = XL;   ZL = XL;
    xlim(ax,XL);  ylim(ax,YL);  zlim(ax,ZL);
    daspect(ax,[1 1 1]);                                   % 不用 axis equal（會撐開 limits）
    grid(ax,'off');  box(ax,'off');
    view(ax, AZ, EL);                                      % 先設 view（draw_box_edges3 用 campos）

    tk = [-100 0 100];                                     % 三軸同刻度：奇數個、等距、不含端點
    set(ax, 'XTick',tk, 'YTick',tk, 'ZTick',tk, ...
            'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, 'TickLength',[.018 .018]);
    ax.XAxis.TickLabelRotation = 0;                        % 3D 預設會沿軸旋轉標籤 → 強制轉正
    ax.YAxis.TickLabelRotation = 0;
    ax.ZAxis.TickLabelRotation = 0;
    % 純幾何示意圖 → 不放軸標題（規則 figure-style「純幾何示意圖：不放軸標題」）

    draw_box_edges3(ax, XL, YL, ZL, LWBOX);                % 手動框邊（必須在 hold off 之前）
    hold(ax,'off');

    %% ---- 輸出 --------------------------------------------------------------
    outf = fullfile(figdir, sprintf('sphere_grid_3d_R%d_%s.png', R_um, tag));
    exportgraphics(fig, outf, 'Resolution', 150);
    fprintf('saved %s\n', outf);
end

% ============================================================================
function draw_box_edges3(ax, XL, YL, ZL, lw)
% 立方 12 邊黑框，省掉「離相機最近角」相連的 3 邊（含會橫穿內部那條）→ 剩 9 邊。
%   （el >= 0 的標準開口箱；本圖 el = 20 > 0。）
    C = [XL([1 1 1 1 2 2 2 2]).', YL([1 1 2 2 1 1 2 2]).', ZL([1 2 1 2 1 2 1 2]).'];
    E = [1 2; 1 3; 1 5; 2 4; 2 6; 3 4; 3 7; 4 8; 5 6; 5 7; 6 8; 7 8];
    cp = campos(ax);
    [~, near] = min(sum((C - cp).^2, 2));
    for m = 1:size(E,1)
        if any(E(m,:) == near), continue; end              % 省最近角的 3 邊
        plot3(ax, C(E(m,:),1), C(E(m,:),2), C(E(m,:),3), 'k-', 'LineWidth', lw);
    end
end
