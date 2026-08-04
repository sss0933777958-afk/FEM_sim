function plot_sensor_cylinder_B_3d(pidx, FACE, NPT_IN)
% plot_sensor_cylinder_B_3d -- SOFF=4.572mm 的 Hall sensor 圓柱 + 柱內 B 向量（3D）
% =========================================================================
%   圓柱**依自身軸向 n+ 傾斜繪製**（不是軸對齊）：底面圓心 = sensor_pos、沿 +n̂ 高 H。
%   幾何與 build_V_matrix 的 sensor_geometry 完全相同（ell*e1 + soff*e2 + AIR*n̂）。
%
%   ⚠ **本圖為內插**（使用者 2026-08-04 指定）：Maxwell 的 sensor 區匯出格距 0.1mm，
%     而 sensor 圓柱 R=0.15mm / H=0.10mm —— 柱內只有 ~7-14 個真格點，畫不出向量場。
%     故以 scatteredInterpolant（取源鄰域 1.5mm，與 build_V_matrix 的 'scattered' 同法）
%     內插到柱內 **100 個撒點**。撒點分布與 build_V_matrix 一致（a~U(0,H)、r=R√U、θ~U(0,2π)）。
%
%   場來源：Maxwell `B_voltage_p<k>.fld`（sensor 區粗格）。號誌套 cfg.s_source（全 source）。
%   風格：figure-style 選項①粗體框 + **3D 變體 A**（box off + 手動框邊省最近角 3 邊、
%         daspect([1 1 1])、三軸同刻度 5 個等距、兩端留白 = 間距、字級 36）。
%   輸出 → figures/paper_fig/Section3_A/sensor_cylinder_B_3d_P<k>.png（覆蓋迭代）。
% =========================================================================
    if nargin < 1
        plot_sensor_cylinder_B_3d(1);  plot_sensor_cylinder_B_3d(2);  return;
    end
    % [ADDED 2026-08-04] FACE：sensor 貼哪個面
    %   'cone'（預設）= 既有定案位置（下極貼半切件的底錐面、上極貼完整錐面）
    %   'flat'        = **下極的平切上表面**（法線 +z、沿水平極軸走 SOFF、離面 AIR）；僅下極有此面
    if nargin < 2 || isempty(FACE), FACE = 'cone'; end
    clc;
    here   = fileparts(mfilename('fullpath'));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    MW = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    addpath(fullfile(MW,'function'));  addpath(fullfile(MW,'common_path'));

    cfg = model_config('long2016_hexapole_halfcut','tip40um');
    SOFF = 4.572e-3;  RS = 0.15e-3;  HS = 0.10e-3;   % sensor 圓柱（與 main.m 預設一致）
    NPT  = 100;                                       % 柱內內插點數（預設；可由 NPT_IN 覆寫）
    if nargin >= 3 && ~isempty(NPT_IN), NPT = NPT_IN; end   % [ADDED 2026-08-04] 點數可調
    %   取樣誤差 = (sigma/mu)/sqrt(NPT)；本例 sigma/mu ≈ 2.3% → 100 點 0.22%、500 點 0.10%
    RLOC = 1.5e-3;                                    % 內插取源鄰域半徑（同 cfg.sensor_r_loc）
    FS   = 36;

    % ---- sensor 幾何（WP 框）：與 build_V_matrix/sensor_geometry 同式 ----
    [sc, nh] = sensor_geom(cfg, pidx, SOFF, FACE);
    ci = sc + [0;0;cfg.SPH_OFST];                     % WP 框 → Maxwell/ANSYS 框

    % ---- 柱內撒點（與 build_V_matrix 同分布）----
    rng(0);
    t1 = [-nh(2); nh(1); 0];  if norm(t1) < 1e-9, t1 = [1;0;0]; end
    t1 = t1/norm(t1);  t2 = cross(nh, t1);
    a  = HS * rand(NPT,1);
    rr = RS * sqrt(rand(NPT,1));
    th = 2*pi * rand(NPT,1);
    Q  = ci.' + a.*nh.' + (rr.*cos(th)).*t1.' + (rr.*sin(th)).*t2.';   % NPT×3（Maxwell 框）

    % ---- 載入 Maxwell sensor 區場 → 內插到撒點（有快取就用；.fld 一顆 ~2GB 讀很久）----
    cachef = fullfile(here, sprintf('sensor_cyl_B_P%d_%s_n%d.mat', pidx, lower(FACE), NPT));
    if exist(cachef,'file')
        L = load(cachef,'Q','B');  Q = L.Q;  B = L.B;
        fprintf('P%d loaded cache %s\n', pidx, cachef);
    else
        fdir = cfg.fld_dir;                                  % Maxwell 匯出根（同 extract_maxwell_data）
        if isfield(cfg,'fld_variant_subdir') && cfg.fld_variant_subdir
            fdir = fullfile(fdir, cfg.default_variant);
        end
        d = import_maxwell_fld(fullfile(fdir, cfg.fld_files_voltage{pidx}));
        m = abs(d.x-ci(1))<=RLOC & abs(d.y-ci(2))<=RLOC & abs(d.z-ci(3))<=RLOC;
        P = [d.x(m), d.y(m), d.z(m)];
        fprintf('鄰域(%.1fmm 立方)內格點: %d\n', RLOC*1e3, size(P,1));
        Fx = scatteredInterpolant(P, d.bx(m), 'linear', 'none');
        Fy = scatteredInterpolant(P, d.by(m), 'linear', 'none');
        Fz = scatteredInterpolant(P, d.bz(m), 'linear', 'none');
        clear d;
        s = cfg.s_source(pidx);                              % 全 source 號誌
        B = s * [Fx(Q), Fy(Q), Fz(Q)] * 1e3;                 % T → mT
        ok = all(isfinite(B),2);  Q = Q(ok,:);  B = B(ok,:);
        save(cachef, 'Q', 'B');   fprintf('P%d saved cache %s\n', pidx, cachef);
    end
    Bm = vecnorm(B,2,2);
    Bn = B*nh;                                               % 法向分量（sensor 讀的就是它）
    fprintf('P%d [%s]: n=%d 點  |B| %.4f~%.4f mT (mean %.4f)  <B.n>=%.4f mT   對準 %.2f deg  -> V=%.1f mV\n', ...
            pidx, FACE, size(Q,1), min(Bm), max(Bm), mean(Bm), mean(Bn), ...
            acosd(mean(Bn)/mean(Bm)), cfg.S_hall*mean(Bn));

    % ---- 座標改回 WP 框、換 mm ----
    Qw = (Q - [0 0 cfg.SPH_OFST]) * 1e3;
    cw = sc*1e3;   gc = cw + nh*(HS/2)*1e3;                  % 圓柱幾何中心（底面 + H/2）

    % ---- 軸範圍/刻度：**3 個**等距 tick、步長 0.1mm、兩端留白 = 間距（三軸同刻度）----
    %   tick 中心優先取到 0.1（標籤短、不擠）；若圓柱因此被框裁到才退回 0.01。
    %   圓柱沿 u 軸的半幅 = R·sqrt(1-(u·n)^2) + (H/2)·|u·n|
    st = 0.1;   LIM = zeros(3,2);  TK = zeros(3,3);
    for k = 1:3
        un  = nh(k);                                          % u·n（u = e_k）
        ext = (RS*sqrt(max(0,1-un^2)) + (HS/2)*abs(un))*1e3;  % 該軸半幅 [mm]
        c0  = round(gc(k), 1);
        if abs(gc(k)-c0) + ext > 2*st - 1e-9, c0 = round(gc(k), 2); end   % 裝不下才用細的
        TK(k,:)  = c0 + (-1:1)*st;
        LIM(k,:) = [TK(k,1)-st, TK(k,end)+st];
    end

    % ---- 畫（manual pixel 佈局，確保 colorbar 標題不被裁）----
    fig = figure('Color','w','Units','pixels','Position',[60 40 1320 980]);
    ax  = axes(fig,'Units','pixels','Position',[190 150 830 790]);  hold(ax,'on');

    draw_cylinder(ax, cw, nh, RS*1e3, HS*1e3);               % 依 n+ 傾斜的圓柱

    % [MODIFIED 2026-08-04] 色階**貼齊資料範圍**（原本向外進位到 0.2mT，兩端色階用不到）
    CLO  = min(Bm);   CLIM = max(Bm);
    nb = 28;  edg = linspace(CLO, CLIM, nb+1);  cmap = turbo(nb);
    lmin = 0.015;  lmax = 0.045;                             % 箭頭長度 [mm]（壓縮映射；不超過柱半徑 1/3）
    len  = lmin + (lmax-lmin)*(Bm/max(Bm)).^0.35;
    U    = B./Bm .* len;
    for k = 1:nb
        if k < nb, sel = Bm>=edg(k) & Bm<edg(k+1); else, sel = Bm>=edg(k); end
        if any(sel)
            quiver3(ax, Qw(sel,1),Qw(sel,2),Qw(sel,3), U(sel,1),U(sel,2),U(sel,3), ...
                    0, 'Color',cmap(k,:), 'LineWidth',2.0, 'MaxHeadSize',0.5);
        end
    end

    grid(ax,'off');  box(ax,'off');  daspect(ax,[1 1 1]);
    xlim(ax,LIM(1,:)); ylim(ax,LIM(2,:)); zlim(ax,LIM(3,:));
    view(ax, -37.5, 22);
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',3.0,'TickLength',[.018 .018]);
    set(ax,'XTick',TK(1,:),'YTick',TK(2,:),'ZTick',TK(3,:));
    % [MODIFIED 2026-08-04] 使用者要求拿掉「軸標題」（x/y/z (mm)）；**刻度數字保留**。
    colormap(ax,turbo);  clim(ax,[CLO CLIM]);
    cb = colorbar(ax,'Units','pixels');
    cb.Label.Interpreter = 'latex';                          % ⚠ Interpreter 要在 String 之前設
    cb.Label.String = '$\mathbf{\|b\|\;(mT)}$';  cb.Label.FontSize = FS;   % [MODIFIED] norm 雙豎線 ‖b‖
    cb.FontSize = FS;  cb.FontWeight = 'bold';
    cb.Position = [1055 150 26 790];                         % 明確定位，標題留在圖內
    draw_box_edges(ax, LIM, 3.0);                            % 變體 A：省最近角 3 邊（須在 hold off 前）
    ax.Toolbar.Visible = 'off';  hold(ax,'off');

    fsfx = ''; if strcmpi(FACE,'flat'), fsfx = '_flat'; end     % 'cone' 沿用原檔名
    out = fullfile(figdir, sprintf('sensor_cylinder_B_3d_P%d%s.png', pidx, fsfx));
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function [pos, n] = sensor_geom(cfg, i, SOFF, FACE)
% 'cone'：與 build_V_matrix 的 local sensor_geometry 同式（下極貼半切件底錐面、上極貼完整錐面）。
% 'flat'：[ADDED] 下極的**平切上表面** —— 半切面是通過極軸的水平面，故沿水平極軸走 SOFF、
%         法線 +z（由鋼往上出）。上極無此面，傳 'flat' 會報錯。
    if nargin < 4 || isempty(FACE), FACE = 'cone'; end
    be = atan2(cfg.POLE_R, cfg.POLE_CONE_LEN);      % 半錐角
    ps = atan2(cfg.R_norm_z, cfg.R_norm_xy);        % magic angle
    iu = cfg.upper_incline;                          % 上極錐軸傾角
    AIR = 0.41e-3;
    dr = @(el,az) [cos(el)*cos(az); cos(el)*sin(az); sin(el)];
    th = cfg.pole_angles(i)*pi/180;
    if strcmpi(FACE,'flat')
        assert(cfg.pole_is_lower(i), 'FACE=''flat'' 只適用下極（上極是完整錐、無平切面）');
        e1 = dr(-ps, th);                            % → 極尖
        e2 = dr(0, th);                              % 平切面內方向 = 水平極軸
        n  = [0; 0; 1];                              % 平切面外法線（鋼在下方）
    elseif cfg.pole_is_lower(i)
        e1 = dr(-ps, th);  e2 = dr(-be, th);        n = dr(-be-pi/2, th);
    else
        e1 = dr(+ps, th);  e2 = dr(iu+be, th);      n = dr(iu+be+pi/2, th);
    end
    pos = cfg.R_norm*e1 + SOFF*e2 + AIR*n;
end

% ============================================================================
function draw_cylinder(ax, c, n, R, H)
% 依軸向 n 傾斜的圓柱：側面 + 上下底圓（參考 NTU sensor_area_B_3d 樣式）。
    t1 = [-n(2); n(1); 0];  if norm(t1) < 1e-9, t1 = [1;0;0]; end
    t1 = t1/norm(t1);  t2 = cross(n, t1);
    ph = linspace(0, 2*pi, 80);
    ring = @(h) (c + h*n).' + R*(cos(ph).'*t1.' + sin(ph).'*t2.');
    r0 = ring(0);  r1 = ring(H);
    surf(ax, [r0(:,1) r1(:,1)], [r0(:,2) r1(:,2)], [r0(:,3) r1(:,3)], ...
         'FaceColor',[1 0.82 0.82], 'FaceAlpha',0.28, 'EdgeColor','none');
    for r = {r0, r1}
        plot3(ax, r{1}(:,1), r{1}(:,2), r{1}(:,3), '-', 'Color',[0.55 0 0], 'LineWidth',2.2);
    end
end

% ============================================================================
function draw_box_edges(ax, LIM, lw)
% 3D 變體 A：手動畫 12 邊、**省略離相機最近角**相連的 3 邊（留 9 邊全粗）。
    [X,Y,Z] = ndgrid(LIM(1,:), LIM(2,:), LIM(3,:));
    C = [X(:) Y(:) Z(:)];
    cp = campos(ax);
    [~, near] = min(sum((C - cp).^2, 2));
    for i = 1:8
        for j = i+1:8
            if nnz(abs(C(i,:)-C(j,:)) > 1e-12) ~= 1, continue; end   % 只連相鄰角
            if i == near || j == near, continue; end                 % 省最近角 3 邊
            plot3(ax, C([i j],1), C([i j],2), C([i j],3), 'k-', 'LineWidth', lw);
        end
    end
end
