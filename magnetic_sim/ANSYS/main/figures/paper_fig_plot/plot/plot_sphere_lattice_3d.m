function plot_sphere_lattice_3d()
% plot_sphere_lattice_3d -- long2016 半切六極 paper 圖:ℓ̂ 球殼 + 6 極尖端(actuator frame)
% =========================================================================
%   **全場景轉進 actuator frame**(R_act 旋轉)→ box 三軸 = 致動軸 x_a/y_a/z_a、
%   6 根極落在 ±x_a/±y_a/±z_a 上(P1=+x P2=-x P3=+y P4=-y P5=+z P6=-z)。view(65,25)。
%   內容:ℓ̂ 線框球殼 + 6 顆**小磁極尖端**(~40µm fillet 圓角+短錐,draw_tip,取代原紅點磁荷)
%          + WP→尖端 虛線 + 藍色 x_a/y_a/z_a 箭頭 + P1..P6 + WP 十字。
%   資料:calibration data/long2016.../fit_fixl_R150um_gap_200um.mat(ell)。
%   單位 mm;框 = box on + daspect;tick -0.5:0.5:0.5(奇數含原點、不擠角落);camlight 打尖端。
%   輸出 → figures/paper_fig/sphere_lattice_3d.png。
% =========================================================================
    clc;
    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';

    % ---- 載 fix(無 bias)fit → ell;轉 actuator frame;單位 mm ----
    S = load(fullfile(CAL,'data','long2016_hexapole_halfcut','.mat','fit_fixl_R150um_gap_200um.mat'),'ell');
    ellm = S.ell / 1000;                                 % mm(≈0.868)
    c   = mt_constants();
    tip  = [c.pole_tip_x; c.pole_tip_y; c.pole_tip_z_wp];
    dhat = tip ./ vecnorm(tip);                          % 3×6 measure 單位方向
    pc   = ellm * dhat;                                  % 6 尖端落在 ℓ̂ 殼上 [mm]（measure frame, 磁角位置）
    fprintf('ell=%.3f mm (measure frame, tips at magic angle)\n', ellm);

    % ---- figure ----
    fig = figure('Color','w','Position',[80 80 1000 940]);  ax=axes(fig);  hold on;

    % 固定 ℓ̂ 半透明球殼
    [sx,sy,sz] = sphere(40);
    surf(ellm*sx, ellm*sy, ellm*sz, 'FaceColor',[0.80 0.86 0.93], 'FaceAlpha',0.06, ...   % #63 淡藍透明殼
         'EdgeColor',[0.62 0.68 0.80], 'EdgeAlpha',0.45, 'LineWidth',0.3, 'FaceLighting','none');

    % 中央綠色網格小球(R150 取樣/WP 球;#107 風格)
    Rg = 0.15;  [gx,gy,gz] = sphere(24);
    surf(Rg*gx, Rg*gy, Rg*gz, 'FaceColor',[0.55 0.80 0.55], 'FaceAlpha',0.18, ...
         'EdgeColor',[0.20 0.60 0.32], 'EdgeAlpha',0.75, 'LineWidth',0.4, 'FaceLighting','none');

    % ---- 各極軸 + 尖端 apex(下極水平半切、上極傾斜全錐;龍飛模型)----
    tipcol = [0.50 0.53 0.60];
    beta = atan2(c.POLE_R, c.POLE_CONE_LEN);   % 半錐角 ≈ 11.3°
    inc  = c.upper_incline;                    % 上極錐軸傾角 ≈ 36.59°
    Lc   = 0.12;  retr = 0.05;                 % 尖端錐長 / 往外退量 [mm]
    islow = logical(c.pole_is_lower);
    A6=zeros(3,6); U6=A6; V6=A6; apex=A6;
    for k = 1:6
        th = c.pole_angles(k)*pi/180;
        if islow(k)
            a = [cos(th); sin(th); 0];                          % 下極:水平徑向軸
            v = [0;0;-1];  u = cross(a,v);  u = u/norm(u);       % 半切平頂在 z=tip_z
        else
            a = [cos(inc)*cos(th); cos(inc)*sin(th); sin(inc)];  % 上極:傾斜向上軸
            u = cross(a,[0;0;1]);  u = u/norm(u);  v = cross(a,u);
        end
        A6(:,k)=a; U6(:,k)=u; V6(:,k)=v;  apex(:,k) = pc(:,k)+retr*a;   % 尖端 apex
    end

    % WP → 磁極尖端 虛線(6 條,連到 apex)
    for k = 1:6
        plot3([0 apex(1,k)],[0 apex(2,k)],[0 apex(3,k)], '--','Color',[0.45 0.45 0.45],'LineWidth',1.2);
    end

    % actuator 三軸(x_a/y_a/z_a;藍;沿虛線方向指向 P1/P3/P5 尖端 apex)——中間座標系
    ca=[0.10 0.35 0.85];  apole=[1 3 5];
    for k = 1:3
        di = apex(:,apole(k))/norm(apex(:,apole(k)));   % 對齊虛線(指向該極尖端 apex)
        v = 0.50*ellm*di;
        quiver3(0,0,0, v(1),v(2),v(3), 0, 'Color',ca,'LineWidth',3.0,'MaxHeadSize',0.9);
    end

    % 6 顆磁極尖端(龍飛模型:下極 P1/P3/P6 水平半切、上極 P2/P4/P5 傾斜全錐)
    for k = 1:6
        draw_pole(apex(:,k), A6(:,k), U6(:,k), V6(:,k), 0.04, beta, Lc, islow(k), tipcol);
    end

    % WP 十字
    plot3(0,0,0,'k+','MarkerSize',13,'LineWidth',2.2);

    % ---- 框 / 視角(規則變體 A:同尺度立方 → box off + 手動 draw_box_edges)----
    bh = 0.9;                                            % mm（再壓緊、球殼幾乎貼框）
    grid off; box off; daspect([1 1 1]);
    xlim([-bh bh]); ylim([-bh bh]); zlim([-bh bh]);
    view(65,25);  camlight('headlight'); lighting gouraud; material dull;   % 先設 view(draw_box_edges 用 campos)
    set(ax,'FontSize',36,'FontWeight','bold','LineWidth',3.0,'TickLength',[0.025 0.025]);
    % tick 取 xlim 內部、奇數個(含原點),不含 ±極值端點 → box 角落不擠 tick/數字
    set(ax,'XTick',-0.5:0.5:0.5,'YTick',-0.5:0.5:0.5,'ZTick',-0.5:0.5:0.5);
    draw_box_edges(bh, 3.0);                             % 手動框邊:9 邊(省最遠角 3 邊)全黑等粗
    % 拿掉 box x/y/z 軸字母標籤(只保留刻度數字)
    ax.Clipping='off';  ax.Toolbar.Visible='off';  hold off;   % 貼邊標籤不被 box 切

    out = fullfile(figdir,'sphere_lattice_3d.png');
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function draw_box_edges(bh, lw)
% box off + 手動畫立方邊全黑等粗;省略「離相機最近角」相連 3 邊(前面 3 條→開口)→ 9 邊。
    s = [-bh bh];  [Xc,Yc,Zc] = ndgrid(s,s,s);  C = [Xc(:) Yc(:) Zc(:)];   % 8 角
    E = [];                                                                 % 12 邊(差一個座標)
    for i = 1:8, for j = i+1:8
        if nnz(abs(C(i,:)-C(j,:)) > 0) == 1, E = [E; i j]; end
    end, end
    cp = campos;  d = sum((C - cp).^2, 2);  [~,near] = min(d);              % 最近角(前面)
    E = E(~(E(:,1)==near | E(:,2)==near), :);                               % 省最近角 3 邊
    for k = 1:size(E,1)
        p1 = C(E(k,1),:);  p2 = C(E(k,2),:);
        plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)], 'k-','LineWidth',lw);
    end
end

% ============================================================================
function draw_pole(tip, a, u, v, rf, beta, L, half, col)
% 磁極尖端:rf 圓角(最前點=tip)+ 切線直錐 → L,合併單一曲面。half=true → 下極半切平頂。
%   tip=尖端(朝 WP), a=外向極軸(單位), u,v=⊥ 軸基底。單位 mm。(龍飛模型幾何)
    ts = rf*(1 + sin(beta));  rt = rf*cos(beta);          % 切點軸距 / 切圓半徑
    psi = linspace(0, pi/2 + beta, 16);                   % 圓角段
    nt  = 20;  t = linspace(ts, L, nt);                   % 錐段
    axoff = [rf*(1 - cos(psi)),  t(2:end)];
    radm  = [rf*sin(psi),        rt + (t(2:end)-ts)*tan(beta)];
    if half, phi = linspace(0, pi, 44); else, phi = linspace(0, 2*pi, 64); end
    N = numel(axoff);  X = zeros(N,numel(phi)); Y = X; Z = X;
    for j = 1:numel(phi)
        rad = u*cos(phi(j)) + v*sin(phi(j));
        P = tip + a*axoff + rad*radm;
        X(:,j) = P(1,:).';  Y(:,j) = P(2,:).';  Z(:,j) = P(3,:).';
    end
    surf(X,Y,Z, 'FaceColor',col, 'FaceAlpha',1.0, 'EdgeColor','none');   % 實心錐(不透明)
    rim = tip + a*axoff(end) + (u*cos(phi) + v*sin(phi))*radm(end);       % 外端 cap
    cen = tip + a*axoff(end);
    patch([cen(1) rim(1,:)],[cen(2) rim(2,:)],[cen(3) rim(3,:)], col, 'EdgeColor','none','FaceAlpha',1.0);
    if half                                                              % 半切平頂 cap
        e0 = tip + a*axoff + u*radm;  ep = tip + a*axoff - u*radm;  poly = [e0, fliplr(ep)];
        patch(poly(1,:),poly(2,:),poly(3,:), col, 'EdgeColor','none','FaceAlpha',1.0);
    end
end
