function plot_ref_planes_3d()
% PLOT_REF_PLANES_3D  獨立的 R=150µm 校正球 + 由下而上水平參考切面 ref_1..ref_N。
%   綠色 R=150µm 線框球（校正取樣區），內部切一疊水平圓盤（球橫切面）：
%     z_i = -150:30:150（N=11），半徑 ρ(z)=√(150²−z²)（赤道 150→兩極 0）。
%     ref_1 = 最底（z=−150）… ref_N = 最頂（z=+150），由下而上編號。
%   每個 ref_i 面 = 之後在該 z 高度做一張 gain/iso 山丘的取樣平面（svd_gain/iso 的 z=0 推廣）。
%   框 = box on + daspect（cube 幾何、view-robust）。輸出 fix_dir/figures/ref_planes_3d.png。

    here   = fileparts(mfilename('fullpath'));
    fixdir = fileparts(fileparts(here));
    figdir = fullfile(fixdir, 'figures');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    R  = 150;                                                      % µm

    fig = figure('Color','w','Position',[80 80 900 860]); hold on;

    %% 綠色 R=150µm 線框球（同 frames 綠球樣式）
    [sx,sy,sz] = sphere(40);
    surf(R*sx, R*sy, R*sz, 'FaceColor',[0.45 0.78 0.55], 'FaceAlpha',0.10, ...
         'EdgeColor',[0.35 0.68 0.45], 'EdgeAlpha',0.30, 'LineWidth',0.3);

    %% 只留中心（赤道 z=0）那一片參考切面（填色+粗邊+引線）+ 紅色 ref 標籤
    th = linspace(0, 2*pi, 80);
    zk = 0;  rk = sqrt(R^2 - zk^2);
    fill3(rk*cos(th), rk*sin(th), zk*ones(size(th)), [0.30 0.50 0.85], ...
          'FaceAlpha',0.45, 'EdgeColor',[0.10 0.22 0.55], 'LineWidth',2.4);
    plot3([rk 178],[0 0],[zk zk], '-', 'Color',[0.10 0.22 0.55], 'LineWidth',1.4);   % 引線
    text(181, 0, zk, 'ref', 'FontSize',17,'FontWeight','bold', ...
         'Color',[0.90 0.10 0.10], 'HorizontalAlignment','left');   % 紅色 ref

    %% 中心點標記（同 frames 的 WP + 標記，不標字）
    plot3(0,0,0,'k+','MarkerSize',13,'LineWidth',2.2);

    %% 框（box on + daspect：cube 幾何、view-robust；視角對齊 frames_lattice_R150）
    bh = 185;
    grid off; box on; daspect([1 1 1]);
    xlim([-bh bh]); ylim([-bh bh]); zlim([-bh bh]);
    view(120, 25);
    set(gca,'FontSize',15,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    set(gca,'XTick',-150:75:150,'YTick',-150:75:150,'ZTick',-150:75:150);
    xlabel('x_m (\mum)','FontWeight','bold');
    ylabel('y_m (\mum)','FontWeight','bold');
    zlabel('z_m (\mum)','FontWeight','bold');
    ax = gca; ax.Toolbar.Visible = 'off';
    hold off;

    out = fullfile(figdir,'ref_planes_3d.png');
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('只留中心切面 + ref；saved %s\n', out);
    close(fig);
end
