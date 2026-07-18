function plot_poles_sensors_3d()
% plot_poles_sensors_3d — NTU_hexapole 六極各出 2 張 3D 圖：全極輪廓+sensing area 圓柱 / zoom。
%   純幾何（自 pole_params.md 9 段 profile + sensor 慣例建），世界座標 measure/world，單位 mm。
%   6 極（paper 名 / 方位 / 層）：P1 0°下、P3 120°下、P6 240°下、P5 60°上、P2 180°上、P4 300°上。
%   每極 sensing area 圓柱 R0.15×H0.10 @ pole-local(13.15,-0.99) 旋轉；z=板底+0.66（板頂上 0.41mm，比照 mt_constants）。
%   每極 2 張：<P>_full（全極輪廓+圓柱）、<P>_zoom（放大到圓柱）。共 12 張。
%   風格（粗體框規則）：扁平幾何(z<<xy)依規則用 **box on 標準 3D 框** + pbaspect 讀性 +
%     粗體字/LineWidth2 + grid off + tick 減半 + 單位 ()。輸出實檔覆蓋迭代。

    here   = fileparts(mfilename('fullpath'));
    vbase  = fileparts(fileparts(here));
    figdir = fullfile(vbase, 'figures', 'shared');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    addpath(fullfile(vbase,'code','function'));      % [ADDED] 共用 helper：ntu_pole_profile/ntu_draw_plate

    B2 = ntu_pole_profile();                         % Nx2 磁極邊界（pole-local；共用 helper）

    % 6 極：paper 名 / 方位角 / 板底 z
    name = {'P1','P3','P6','P5','P2','P4'};
    ang  = [0 120 240 60 180 300];
    zb   = [8.80 8.80 8.80 9.55 9.55 9.55];
    TH   = 0.25;
    SX=13.15; SY=-0.99; SR=0.15; SH=0.10;            % sensor pole-local + 尺寸

    colLo=[0.30 0.55 0.95]; colHi=[0.95 0.60 0.25]; colS=[0.85 0.10 0.10];

    for p = 1:6
        th = ang(p)*pi/180;  R=[cos(th) -sin(th); sin(th) cos(th)];
        Pw = (R*B2.').';                             % 世界座標輪廓
        z0 = zb(p); z1 = z0+TH;
        sx = SX*cos(th)-SY*sin(th);  sy = SX*sin(th)+SY*cos(th);  sz = z0+0.66;   % 板頂上 0.41mm（= plate_zb+0.66，比照 mt_constants）
        if z0<9.0; colP=colLo; else; colP=colHi; end

        % ---- 圖 1：全極 + sensing area 圓柱 ----
        f=figure('Color','w','Position',[100 100 1180 940]); hold on;
        ntu_draw_plate(Pw(:,1),Pw(:,2),z0,z1,colP,0.30);
        draw_cylinder(sx,sy,sz,SR,SH,colS);
        style_box(gca, ...
            [min(Pw(:,1))-1.5 max(Pw(:,1))+1.5], [min(Pw(:,2))-1.5 max(Pw(:,2))+1.5], [z0-0.3 sz+0.4], ...
            -35, 28, 0.5);                           % full 仰角再拉低（偏側視）；z 軸短→0.5 間距 tick
        text(sx,sy,sz+0.12,'sensing area','FontSize',15,'FontWeight','bold','Color',colS, ...
            'HorizontalAlignment','center','VerticalAlignment','bottom');
        export_fig(f, fullfile(figdir,[name{p} '_full.png']));

        % ---- 圖 2：zoom 到 sensing area 圓柱 ----
        f=figure('Color','w','Position',[100 100 1180 940]); hold on;
        ntu_draw_plate(Pw(:,1),Pw(:,2),z0,z1,colP,0.30);
        draw_cylinder(sx,sy,sz,SR,SH,colS);
        style_box(gca, [sx-1.3 sx+1.3], [sy-1.4 sy+1.4], [z1-0.15 sz+0.22], -42, 16, 0.5);
        text(sx+0.30,sy,sz,' sensing area','FontSize',15,'FontWeight','bold','Color',colS);
        export_fig(f, fullfile(figdir,[name{p} '_zoom.png']));
    end
    fprintf('wrote 12 figs (P1..P6 x full/zoom) to %s\n', figdir);
end

% ============================================================================
function style_box(ax, xl, yl, zl, az, el, zstep)
% 粗體框規則（扁平幾何 → box on 標準 3D 框 + daspect 真實比例）
    xlim(ax,xl); ylim(ax,yl); zlim(ax,zl);
    grid(ax,'off'); box(ax,'on');                    % box on 標準 3D 框（粗體框規則；扁平幾何用 B 框法）
    daspect(ax,[1 1 1]);                             % 真實比例（x,y,z 皆 mm）；不拉爆薄板 z
    view(ax,az,el);
    set(ax,'FontSize',15,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    xlabel(ax,'x (mm)','FontWeight','bold'); ylabel(ax,'y (mm)','FontWeight','bold'); zlabel(ax,'z (mm)','FontWeight','bold');
    xt=get(ax,'XTick'); if numel(xt)>3; set(ax,'XTick',xt(1:2:end)); end   % tick 減半
    yt=get(ax,'YTick'); if numel(yt)>3; set(ax,'YTick',yt(1:2:end)); end
    zv = ceil(zl(1)/zstep)*zstep : zstep : floor(zl(2)/zstep)*zstep;       % z 稀疏 tick（真實比例 z 軸短→避免重疊；full 用 1.0/zoom 用 0.5）
    if numel(zv)>=1; set(ax,'ZTick',zv); end
    camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull');
    ax.Toolbar.Visible='off';
end

function export_fig(f, outpng)
    exportgraphics(f, outpng, 'Resolution', 150); close(f);
end

% ============================================================================
% [MODIFIED 2026-07-17] build_pole_profile / arc / draw_plate 已抽到共用 helper
%   code/function/{ntu_pole_profile,ntu_draw_plate}.m（原本本檔與 plot_p1_sensing_disk.m
%   各複製一份、無法跨檔呼叫）。本檔改為呼叫共用版，行為不變。

% [MODIFIED 2026-07-17] sensor 圓柱改畫成「底面在 cz、往 +n 長 H」，與 build_V_matrix 的取樣一致
%   （a = axial_tol*rand → pts = c + a*n̂，即 c 是**底面中心**、不是幾何中心）。
%   原本畫成以 cz 為中心（cz±H/2）→ 與實際取樣位置差 H/2 = 0.05mm。
function draw_cylinder(cx, cy, cz, R, H, col)
    [Xc,Yc,Zc]=cylinder(R,24); Xc=Xc+cx; Yc=Yc+cy; Zc=Zc*H+cz;
    surf(Xc,Yc,Zc,'FaceColor',col,'EdgeColor','none','FaceAlpha',1.0);
    th=linspace(0,2*pi,25);
    for zc=[cz, cz+H]
        patch(cx+R*cos(th),cy+R*sin(th),zc*ones(size(th)),col,'EdgeColor','none');
    end
end
