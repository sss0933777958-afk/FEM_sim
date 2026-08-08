function plot_ref_planes_3d()
% plot_ref_planes_3d -- paper 圖(Section4_C):R=150µm 取樣球殼 + 三 actuator 參考切面 + 6 極方向
% =========================================================================
%   全場景在 R=0.15mm(150µm 取樣殼)上:綠色線框球殼 + 三過球心大圓盤(actuator 座標系
%   x_a=P1, y_a=P3, z_a=P5,magic-angle 兩兩正交):ref_{x_ay_a}(藍)/ref_{y_az_a}(橘)/ref_{x_az_a}(紫)。
%   6 磁極方向落在 150µm 殼面(紅點 + P1..P6)。單位 mm。
%   規則:粗體框圖變體 A(同尺度立方 → box off + daspect + 手動 draw_box_edges 省最近角3邊、字體36)、
%         圖例移上方水平(含綠殼說明)、每軸 3 個內縮 tick。
%   輸出 → figures/paper_fig/Section4_C/ref_planes_3d.png。
% =========================================================================
    clc;
    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section4_C');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    % [MODIFIED 2026-08-08] 脫離 backup（規則 no-backup-data）→ live 樹。
    CALROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CALROOT,'function'), fullfile(CALROOT,'utils'), fullfile(CALROOT,'common_path'));
    c = model_config('long2016_hexapole_halfcut','tip40um');

    R = 0.15;                                                     % mm（= 150µm 取樣殼半徑）

    % actuator 座標軸(measure frame):x_a=P1, y_a=P3, z_a=P5 尖端方向(magic-angle 正交)
    tip  = [c.pole_tip_x; c.pole_tip_y; c.pole_tip_z_wp];
    dhat = tip ./ vecnorm(tip);                                  % 3×6 單位方向
    ea_x = dhat(:,1);  ea_y = dhat(:,3);  ea_z = dhat(:,5);

    fig = figure('Color','w','Position',[80 60 1040 1020]);  ax = axes(fig);  hold on;

    % 綠色 150µm 線框球殼(取樣區)
    [sx,sy,sz] = sphere(40);
    surf(R*sx, R*sy, R*sz, 'FaceColor',[0.45 0.78 0.55], 'FaceAlpha',0.10, ...
         'EdgeColor',[0.35 0.68 0.45], 'EdgeAlpha',0.30, 'LineWidth',0.3, 'FaceLighting','none');

    % 三 actuator 參考切面(過球心大圓盤,沿 actuator 軸)
    th = linspace(0, 2*pi, 90);  cth = cos(th);  sth = sin(th);
    c_xy=[0.30 0.50 0.85]; e_xy=[0.10 0.22 0.55];      % x_a-y_a 藍
    c_yz=[0.92 0.60 0.24]; e_yz=[0.78 0.42 0.10];      % y_a-z_a 橘
    c_xz=[0.58 0.40 0.75]; e_xz=[0.42 0.20 0.62];      % x_a-z_a 紫
    draw_ref_plane(R*(ea_x*cth + ea_y*sth), c_xy, e_xy);
    draw_ref_plane(R*(ea_y*cth + ea_z*sth), c_yz, e_yz);
    draw_ref_plane(R*(ea_x*cth + ea_z*sth), c_xz, e_xz);

    % 圖例 proxy(三面 rim 色粗線 + 綠殼)
    hL(1) = plot3(nan,nan,nan,'-','Color',e_xy,'LineWidth',6);
    hL(2) = plot3(nan,nan,nan,'-','Color',e_yz,'LineWidth',6);
    hL(3) = plot3(nan,nan,nan,'-','Color',e_xz,'LineWidth',6);

    % （6 磁極方向紅點 + P1–P6 編號皆已移除）

    % WP 中心十字
    plot3(0,0,0,'k+','MarkerSize',15,'LineWidth',2.6);

    % ---- 框/視角(變體 A:同尺度立方 → box off + daspect + 手動框邊)----
    bh = 1.2*R;                                                 % mm（收緊、球殼填滿框、減少上方空白）
    grid off; box off; daspect([1 1 1]);
    xlim([-bh bh]); ylim([-bh bh]); zlim([-bh bh]);
    view(120, 25);                                              % 先設 view(draw_box_edges 用 campos)
    set(ax,'Position',[0.11 0.07 0.79 0.79]);                 % 軸留邊(刻度數字不裁)、框往上撐縮短與圖例距離
    set(ax,'FontSize',36,'FontWeight','bold','LineWidth',3.0,'TickLength',[0.022 0.022]);   % 刻度小線標保留
    set(ax,'XTick',[-0.15 0 0.15],'YTick',[-0.15 0 0.15],'ZTick',[-0.15 0 0.15]);   % 每軸 3 個 tick：-0.15/0/0.15
    % （軸標題 x_m/y_m/z_m 皆移除；只留刻度數字）
    draw_box_edges(bh, 3.0);                                     % 手動框邊:省最近角 3 邊
    ax.Clipping='off';  ax.Toolbar.Visible='off';

    % 圖例:上方水平、手動貼近圖框上緣(不壓到框)
    lg = legend(hL(1:3), {'ref_{x_ay_a}','ref_{y_az_a}','ref_{x_az_a}'}, ...
                'Interpreter','tex', 'Orientation','horizontal', ...
                'FontSize',22, 'Box','on', 'EdgeColor',[0.30 0.30 0.30]);
    lg.FontWeight = 'bold';  lg.Units = 'normalized';  drawnow;
    lg.Position(1) = 0.5 - lg.Position(3)/2;                    % 水平置中
    lg.Position(2) = 0.905;                                     % 貼近圖框上緣、留小間距(不壓到框)
    hold off;

    out = fullfile(figdir,'ref_planes_3d.png');
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s (R=%.2f mm shell + 3 actuator ref planes + 6 pole dirs)\n', out, R);
    close(fig);
end

% ============================================================================
function draw_ref_plane(P, cface, cedge)
    fill3(P(1,:), P(2,:), P(3,:), cface, 'FaceAlpha',0.30, 'EdgeColor',cedge, 'LineWidth',2.6);
end

% ============================================================================
function draw_box_edges(bh, lw)
% box off + 手動立方邊全黑等粗;省「離相機最近角」相連 3 邊(前面開口)→ 9 邊。
    s = [-bh bh];  [Xc,Yc,Zc] = ndgrid(s,s,s);  C = [Xc(:) Yc(:) Zc(:)];
    E = [];
    for i = 1:8, for j = i+1:8
        if nnz(abs(C(i,:)-C(j,:)) > 0) == 1, E = [E; i j]; end
    end, end
    cp = campos;  d = sum((C - cp).^2, 2);  [~,near] = min(d);
    E = E(~(E(:,1)==near | E(:,2)==near), :);
    for k = 1:size(E,1)
        p1 = C(E(k,1),:);  p2 = C(E(k,2),:);
        plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)], 'k-','LineWidth',lw);
    end
end
