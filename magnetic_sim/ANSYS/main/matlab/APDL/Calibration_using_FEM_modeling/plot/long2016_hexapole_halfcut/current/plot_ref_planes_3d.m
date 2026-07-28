function plot_ref_planes_3d()
% PLOT_REF_PLANES_3D  R=500µm 取樣球 + 三個 actuator 座標參考切面 + 6 磁極尖端點標記。
%   綠色 R=500µm 線框球（gain/iso 取樣區）內含三過球心大圓盤，沿 actuator 座標系
%   （x_a=P1, y_a=P3, z_a=P5 尖端軸，magic-angle 兩兩正交）：
%   ref_{x_ay_a}(藍)、ref_{y_az_a}(橘)、ref_{x_az_a}(紫)。6 磁極尖端標紅點 + P1–P6 標籤。
%   框 = box on + daspect（cube 幾何、view-robust）。輸出 current_base/figures/ref_planes_3d.png。

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');   % mt_constants
    cnst = mt_constants();
    here   = fileparts(mfilename('fullpath'));
    fixdir = fileparts(fileparts(here));
    figdir = fullfile(fixdir, 'figures','shared');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    R  = 500;                                                      % µm（取樣球半徑 = R_norm）

    %% actuator 座標系軸（measure 座標）：x_a=P1, y_a=P3, z_a=P5 尖端方向（magic-angle 正交）
    tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp]*1e6;   % 3x6 µm（|tip|=500µm）
    dhat = tip ./ vecnorm(tip);
    ea_x = dhat(:,1);  ea_y = dhat(:,3);  ea_z = dhat(:,5);             % actuator 三軸

    fig = figure('Color','w','Position',[80 80 900 860]); hold on;

    %% 綠色 R=500µm 線框球（取樣區）
    [sx,sy,sz] = sphere(40);
    surf(R*sx, R*sy, R*sz, 'FaceColor',[0.45 0.78 0.55], 'FaceAlpha',0.10, ...
         'EdgeColor',[0.35 0.68 0.45], 'EdgeAlpha',0.30, 'LineWidth',0.3);

    th = linspace(0, 2*pi, 90);  cth = cos(th);  sth = sin(th);
    %% 三個 actuator 座標參考切面（過球心大圓盤，沿 actuator 軸對；標籤改右上角圖例）
    c_xy = [0.30 0.50 0.85];  e_xy = [0.10 0.22 0.55];    % x_a-y_a 藍
    c_yz = [0.92 0.60 0.24];  e_yz = [0.78 0.42 0.10];    % y_a-z_a 橘
    c_xz = [0.58 0.40 0.75];  e_xz = [0.42 0.20 0.62];    % x_a-z_a 紫
    draw_ref_plane(R*(ea_x*cth + ea_y*sth), c_xy, e_xy);
    draw_ref_plane(R*(ea_y*cth + ea_z*sth), c_yz, e_yz);
    draw_ref_plane(R*(ea_x*cth + ea_z*sth), c_xz, e_xz);
    % 圖例 proxy（rim 色粗線代表三面；真圖例在右上角）
    hL(1) = plot3(nan,nan,nan,'-','Color',e_xy,'LineWidth',5);
    hL(2) = plot3(nan,nan,nan,'-','Color',e_yz,'LineWidth',5);
    hL(3) = plot3(nan,nan,nan,'-','Color',e_xz,'LineWidth',5);

    %% 6 磁極尖端點（在 R_norm=500µm 球面）+ P1–P6 標籤（tip/dhat 上方已算；兩三角面已移除）
    plot3(tip(1,:), tip(2,:), tip(3,:), 'o', 'MarkerSize',9, ...
          'MarkerFaceColor',[0.85 0.15 0.35], 'MarkerEdgeColor','k', 'LineWidth',1.0);
    for k = 1:6
        t = tip(:,k) + 90*dhat(:,k);
        text(t(1),t(2),t(3), sprintf('P%d',k), 'FontSize',14,'FontWeight','bold', ...
             'Color',[0.6 0 0.2], 'HorizontalAlignment','center');
    end

    %% 中心點標記（WP）
    plot3(0,0,0,'k+','MarkerSize',13,'LineWidth',2.2);

    %% 框（box on + daspect：cube 幾何、view-robust）
    bh = 1.32*R;
    grid off; box on; daspect([1 1 1]);
    xlim([-bh bh]); ylim([-bh bh]); zlim([-bh bh]);
    view(120, 25);
    set(gca,'FontSize',15,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    set(gca,'XTick',-500:250:500,'YTick',-500:250:500,'ZTick',-500:250:500);
    xlabel('x_m (\mum)','FontWeight','bold');
    ylabel('y_m (\mum)','FontWeight','bold');
    zlabel('z_m (\mum)','FontWeight','bold');
    ax = gca; ax.Toolbar.Visible = 'off';
    lg = legend(hL, {'ref_{x_ay_a}','ref_{y_az_a}','ref_{x_az_a}'}, ...
                'Location','northeast', 'FontSize',15, 'Box','on', 'EdgeColor',[0.30 0.30 0.30]);
    lg.FontWeight = 'bold';
    hold off;

    out = fullfile(figdir,'ref_planes_3d.png');
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('3 actuator ref planes (x_ay_a / y_az_a / x_az_a) + 6 pole points（R=500µm）；saved %s\n', out);
    close(fig);
end

function draw_ref_plane(P, cface, cedge)
% 畫一個過球心大圓盤（3xN 點 P，actuator 面）。標籤改用右上角圖例。
    fill3(P(1,:), P(2,:), P(3,:), cface, 'FaceAlpha',0.30, 'EdgeColor',cedge, 'LineWidth',2.4);
end
