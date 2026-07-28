function plot_ref_planes_3d()
% PLOT_REF_PLANES_3D  R=500µm 取樣球 + 兩正交參考切面 ref_xy/ref_xz + 6 磁極尖端點標記。
%   綠色 R=500µm 線框球（gain/iso 取樣區）內含兩過球心大圓盤：ref_xy(z=0,藍)、ref_xz(y=0,橘)。
%   6 磁極尖端（在 R_norm=500µm 球面）標紅點 + P1–P6 標籤（frames_lattice 風格）。
%   框 = box on + daspect（cube 幾何、view-robust）。輸出 current_base/figures/ref_planes_3d.png。

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');   % hung mt_constants
    cnst = mt_constants();
    here   = fileparts(mfilename('fullpath'));
    fixdir = fileparts(fileparts(here));
    figdir = fullfile(fixdir, 'figures','shared');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    R  = 500;                                                      % µm（取樣球半徑 = R_norm）
    fig = figure('Color','w','Position',[80 80 900 860]); hold on;

    %% 綠色 R=500µm 線框球（取樣區）
    [sx,sy,sz] = sphere(40);
    surf(R*sx, R*sy, R*sz, 'FaceColor',[0.45 0.78 0.55], 'FaceAlpha',0.10, ...
         'EdgeColor',[0.35 0.68 0.45], 'EdgeAlpha',0.30, 'LineWidth',0.3);

    th = linspace(0, 2*pi, 90);
    %% ref_xy：z=0 平面大圓盤（藍）+ 引線 + 標籤
    fill3(R*cos(th), R*sin(th), zeros(size(th)), [0.30 0.50 0.85], ...
          'FaceAlpha',0.32, 'EdgeColor',[0.10 0.22 0.55], 'LineWidth',2.4);
    plot3([R 1.19*R],[0 0],[0 0], '-', 'Color',[0.10 0.22 0.55], 'LineWidth',1.6);
    text(1.23*R, 0, 0, 'ref_{xy}', 'FontSize',17,'FontWeight','bold', ...
         'Color',[0.10 0.22 0.55], 'HorizontalAlignment','left', ...
         'BackgroundColor','w','EdgeColor',[0.10 0.22 0.55],'Margin',3);

    %% ref_xz：y=0 平面大圓盤（橘）+ 引線 + 標籤
    fill3(R*cos(th), zeros(size(th)), R*sin(th), [0.92 0.60 0.24], ...
          'FaceAlpha',0.32, 'EdgeColor',[0.78 0.42 0.10], 'LineWidth',2.4);
    plot3([0 0],[0 0],[R 1.19*R], '-', 'Color',[0.78 0.42 0.10], 'LineWidth',1.6);
    text(0, 0, 1.24*R, 'ref_{xz}', 'FontSize',17,'FontWeight','bold', ...
         'Color',[0.85 0.45 0.10], 'HorizontalAlignment','center', ...
         'BackgroundColor','w','EdgeColor',[0.78 0.42 0.10],'Margin',3);

    %% 6 磁極尖端點（在 R_norm=500µm 球面）+ P1–P6 標籤（frames_lattice 風格）
    tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp]*1e6;   % 3x6 µm（|tip|=500µm）
    dhat = tip ./ vecnorm(tip);
    lo = [1 3 6];  up = [2 4 5];                                          % 上下極三角平面（frames_lattice 風格）
    fill3(tip(1,lo),tip(2,lo),tip(3,lo), [0.90 0.28 0.28], 'FaceAlpha',0.22, 'EdgeColor',[0.75 0.10 0.10], 'LineWidth',1.3);   % 下極 P1/P3/P6（紅）
    fill3(tip(1,up),tip(2,up),tip(3,up), [0.45 0.65 0.95], 'FaceAlpha',0.20, 'EdgeColor',[0.25 0.45 0.85], 'LineWidth',1.3);   % 上極 P2/P4/P5
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
    hold off;

    out = fullfile(figdir,'ref_planes_3d.png');
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('ref_xy + ref_xz + 6 pole points（R=500µm）；saved %s\n', out);
    close(fig);
end
