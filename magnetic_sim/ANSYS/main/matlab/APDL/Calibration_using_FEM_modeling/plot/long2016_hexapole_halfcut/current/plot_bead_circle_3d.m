function plot_bead_circle_3d()
% plot_bead_circle_3d — 確認幾何：磁珠 R=50µm 圓軌跡 @ measure/WP frame xm-ym 平面(z=0)
%   + R≤150µm 球殼（磁珠活動範圍）。純幾何、不需 FEM 場。measure 座標、單位 µm。
%   風格：粗體框 3D 立方（daspect[1 1 1] + box on + 粗體 + 單位 (µm)），view 似附圖。
%   輸出：current_base/figures/shared/bead_circle_confirm.png（figure-output：覆蓋迭代）。

    here   = fileparts(mfilename('fullpath'));            % .../current_base/code/plot
    cbase  = fileparts(fileparts(here));                 % .../current_base
    figdir = fullfile(cbase,'figures','shared');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    outpng = fullfile(figdir,'bead_circle_confirm.png');

    Rc  = 50;    % 圓軌跡半徑 µm
    Rsh = 150;   % 球殼半徑 µm（磁珠活動範圍 R<=150）
    th  = linspace(0,2*pi,300);
    cx  = Rc*cos(th); cy = Rc*sin(th); cz = zeros(size(th));

    f = figure('Color','w','Position',[100 100 1050 950]); hold on;

    % R<=150µm 球殼（半透明綠、磁珠活動範圍）
    [sx,sy,sz] = sphere(46);
    surf(Rsh*sx,Rsh*sy,Rsh*sz,'FaceColor',[0.45 0.78 0.55],'FaceAlpha',0.10, ...
        'EdgeColor',[0.55 0.72 0.58],'EdgeAlpha',0.15,'LineWidth',0.3);

    % R=50µm 圓軌跡（xm-ym 平面 z=0，藍粗線）
    plot3(cx,cy,cz,'-','Color',[0 0 1],'LineWidth',2.6);
    % z=0 平面淡圓盤示意（可選）
    patch(cx,cy,cz,[0 0 1],'FaceAlpha',0.05,'EdgeColor','none');

    % 原點標記（去掉 WP 字眼，只留點）
    plot3(0,0,0,'k.','MarkerSize',20);
    text(Rc*cos(pi/4)+6, Rc*sin(pi/4)+6, 12,'R=50\mum','Color','b','FontSize',13,'FontWeight','bold');

    % ---- 風格：粗體框 3D 立方 ----
    ax=gca; bh=175;
    xlim([-bh bh]); ylim([-bh bh]); zlim([-bh bh]);
    daspect([1 1 1]); grid off; box on;
    view(-37,15);
    set(ax,'FontSize',15,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    xlabel('x_m (\mum)','FontWeight','bold'); ylabel('y_m (\mum)','FontWeight','bold'); zlabel('z_m (\mum)','FontWeight','bold');
    set(ax,'XTick',-150:75:150,'YTick',-150:75:150,'ZTick',-150:75:150);
    camlight('headlight'); lighting gouraud; material dull;
    ax.Toolbar.Visible='off';

    exportgraphics(f,outpng,'Resolution',150); close(f);
    fprintf('wrote %s\n',outpng);
end
