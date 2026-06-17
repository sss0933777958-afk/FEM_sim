%% plot_pole_sensor_placement.m — B_surface placement on each of the 6 poles
%
%  One side-view figure per pole (P1..P6) showing the B_surface (Hall)
%  sensor placement + its positive normal n+.  Prep for the B_S matrix.
%
%  6 poles = 2 geometric types (notation-glossary.md):
%    lower (milled half-cone): P1(0deg) P3(120deg) P6(240deg)  n+ = +z, out of milled flat
%    upper (natural cone, tilted): P2(180deg) P4(300deg) P5(60deg)  n+ = outward normal to cone face
%
%  Output: magnetic_sim/ANSYS/main/figures/long2016_hexapole_halfcut/Bsurf_placement_P#.png

clear; close all;

out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
if ~exist(out_dir,'dir'); mkdir(out_dir); end

poles = struct( ...
  'name', {'P1','P2','P3','P4','P5','P6'}, ...
  'az',   { 0,   180,  120,  300,  60,   240 }, ...
  'layer',{'lower','upper','lower','upper','upper','lower'});

for k = 1:numel(poles)
    p = poles(k);
    if strcmp(p.layer,'lower')
        fig = figure('Color','w','Position',[100 100 1100 460]);
        draw_lower(fig, p);
    else
        fig = figure('Color','w','Position',[100 100 920 820]);
        draw_upper(fig, p);
    end
    out = fullfile(out_dir, sprintf('Bsurf_placement_%s.png', p.name));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('Saved: %s\n', out);
    close(fig);
end

%% ===== lower pole (milled half-cone) =====
function draw_lower(fig, p)
    POLE_R = 3.175;  POLE_CONE_LEN = 15.876;  cyl_end = 30.0;
    half_angle_deg = atan2d(POLE_R, POLE_CONE_LEN);
    s_along = 4.572;  s_offset = 0.41;
    sx = s_along;  sy = s_offset;

    pole_x = [0, POLE_CONE_LEN, cyl_end, cyl_end, 0];
    pole_y = [0, -POLE_R,       -POLE_R, 0,       0];

    ax = axes(fig); hold(ax,'on');
    fill(ax, pole_x, pole_y, [0.78 0.80 0.84], 'EdgeColor',[0.20 0.20 0.25],'LineWidth',1.8);
    plot(ax, [POLE_CONE_LEN POLE_CONE_LEN], [-POLE_R 0], ':','Color',[0.55 0.55 0.60],'LineWidth',0.8);
    plot(ax, [0 cyl_end], [0 0], '-.','Color',[0.75 0.75 0.80],'LineWidth',0.4);

    % sensor
    plot(ax, sx, sy, 'o','MarkerSize',9,'MarkerFaceColor',[0.85 0.10 0.10],'MarkerEdgeColor',[0.55 0.05 0.05]);
    % n+ arrow (out of milled flat, +y == +z global)
    quiver(ax, sx, sy, 0, 1.6, 0, 'Color',[0 0.5 0],'LineWidth',2.4,'MaxHeadSize',1.0);
    text(ax, sx+0.25, sy+1.75, 'n_+  (\perp milled flat, out of steel)', ...
         'FontSize',10,'FontWeight','bold','Color',[0 0.45 0]);

    % dim 4.572
    dy = 1.35;
    plot(ax,[0 s_along],[dy dy],'-','Color',[0.80 0.10 0.10],'LineWidth',1.4);
    plot(ax,[0 0],dy+[-0.1 0.1],'-','Color',[0.80 0.10 0.10],'LineWidth',1.4);
    plot(ax,[s_along s_along],dy+[-0.1 0.1],'-','Color',[0.80 0.10 0.10],'LineWidth',1.4);
    text(ax,s_along/2,dy+0.15,'4.572 mm','HorizontalAlignment','center', ...
         'VerticalAlignment','bottom','Color',[0.65 0.05 0.05],'FontSize',12,'FontWeight','bold');
    % dim 0.41
    dxr = s_along + 1.8;
    plot(ax,[dxr dxr],[0 s_offset],'-','Color',[0.80 0.10 0.10],'LineWidth',1.4);
    plot(ax,dxr+[-0.12 0.12],[0 0],'-','Color',[0.80 0.10 0.10],'LineWidth',1.4);
    plot(ax,dxr+[-0.12 0.12],[s_offset s_offset],'-','Color',[0.80 0.10 0.10],'LineWidth',1.4);
    plot(ax,[sx dxr],[sy sy],'--','Color',[0.85 0.55 0.55],'LineWidth',0.6);
    plot(ax,[sx dxr],[0 0],'--','Color',[0.85 0.55 0.55],'LineWidth',0.6);
    text(ax,dxr+0.25,s_offset/2,'0.41 mm','HorizontalAlignment','left', ...
         'VerticalAlignment','middle','Color',[0.65 0.05 0.05],'FontSize',12,'FontWeight','bold');

    % annotations
    plot(ax,0,0,'k.','MarkerSize',14);
    text(ax,-0.3,-0.35,'Apex (0, 0)','HorizontalAlignment','right', ...
         'VerticalAlignment','top','FontSize',10,'Color',[0.2 0.2 0.2]);
    text(ax,sx+0.35,sy+0.05,'B_{surface}  (Hall sensor, in air)','HorizontalAlignment','left', ...
         'VerticalAlignment','bottom','FontSize',10,'Color',[0.65 0.05 0.05]);
    text(ax,22,0.55,'Sensing surface — milled flat (y = 0, through pole axis)', ...
         'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9.5,'Color',[0.20 0.20 0.30]);
    plot(ax,[22 22],[0.50 0.05],'-','Color',[0.45 0.45 0.50],'LineWidth',0.5);
    text(ax,7.5,-2.4,sprintf('Lower cone face (%.2f\\circ)',half_angle_deg), ...
         'HorizontalAlignment','center','FontSize',9,'Color',[0.25 0.25 0.30],'Rotation',-half_angle_deg);
    text(ax,(POLE_CONE_LEN+cyl_end)/2,-1.55,'Half-cylinder body (R = 3.175 mm)', ...
         'HorizontalAlignment','center','FontSize',9,'Color',[0.25 0.25 0.30]);

    axis(ax,'equal'); xlim(ax,[-2,cyl_end+2]); ylim(ax,[-POLE_R-1.3,3.7]);
    xlabel(ax,'x [mm]  (along pole axis from apex)'); ylabel(ax,'y [mm]');
    title(ax,sprintf('%s — Lower pole (azimuth %d\\circ) — B_{surface} placement', p.name, p.az));
    set(ax,'FontSize',10,'Box','on','Layer','top'); grid(ax,'off');
    ax.Toolbar.Visible = 'off';
end

%% ===== upper pole (natural cone, hexapole tilt) =====
function draw_upper(fig, p)
    POLE_R = 3.175;  POLE_CONE_LEN = 15.876;  cyl_len = 14.124;
    half_angle_deg = atan2d(POLE_R, POLE_CONE_LEN);
    beta = deg2rad(half_angle_deg);
    tilt_deg = 90 - 54.7356;  tilt = deg2rad(tilt_deg);

    ux = cos(tilt);  uy = sin(tilt);
    px = -sin(tilt); py = cos(tilt);
    sx_sl = cos(tilt+beta);  sy_sl = sin(tilt+beta);
    nx = -sin(tilt+beta);    ny = cos(tilt+beta);

    A = [0;0];
    J = POLE_CONE_LEN*[ux;uy];
    T = (POLE_CONE_LEN+cyl_len)*[ux;uy];
    Ju = J+POLE_R*[px;py];  Jl = J-POLE_R*[px;py];
    Tu = T+POLE_R*[px;py];  Tl = T-POLE_R*[px;py];
    pole_x = [A(1),Ju(1),Tu(1),Tl(1),Jl(1),A(1)];
    pole_y = [A(2),Ju(2),Tu(2),Tl(2),Jl(2),A(2)];

    s_slant = 4.572;  s_offset = 0.41;
    foot = s_slant*[sx_sl;sy_sl];
    sxs = foot(1)+s_offset*nx;  sys = foot(2)+s_offset*ny;

    ax = axes(fig); hold(ax,'on');
    fill(ax, pole_x, pole_y, [0.78 0.80 0.84],'EdgeColor',[0.20 0.20 0.25],'LineWidth',1.8);
    plot(ax,[Ju(1) Jl(1)],[Ju(2) Jl(2)],':','Color',[0.55 0.55 0.60],'LineWidth',0.8);
    axis_end = T + 1.5*[ux;uy];
    plot(ax,[A(1) axis_end(1)],[A(2) axis_end(2)],'-.','Color',[0.78 0.78 0.83],'LineWidth',0.4);

    % sensor
    plot(ax,sxs,sys,'o','MarkerSize',9,'MarkerFaceColor',[0.85 0.10 0.10],'MarkerEdgeColor',[0.55 0.05 0.05]);
    % n+ arrow (outward normal to cone slant) — label pulled out with leader
    quiver(ax,sxs,sys,2.0*nx,2.0*ny,0,'Color',[0 0.5 0],'LineWidth',2.4,'MaxHeadSize',1.0);
    ntip = [sxs+2.0*nx, sys+2.0*ny];
    np_lbl = [6.5, 17.5];
    plot(ax,[np_lbl(1) ntip(1)],[np_lbl(2) ntip(2)],'-','Color',[0.45 0.70 0.45],'LineWidth',0.6);
    text(ax,np_lbl(1),np_lbl(2)+0.5,'n_+  (\perp cone face, out of steel)', ...
         'FontSize',10,'FontWeight','bold','Color',[0 0.45 0],'HorizontalAlignment','center');

    % dim 4.572 along slant
    dim_off = 2.6;
    D1a = A+dim_off*[nx;ny];  D1b = foot+dim_off*[nx;ny];
    plot(ax,[D1a(1) D1b(1)],[D1a(2) D1b(2)],'-','Color',[0.80 0.10 0.10],'LineWidth',1.4);
    tl = 0.18;
    plot(ax,D1a(1)+tl*[-sx_sl sx_sl],D1a(2)+tl*[-sy_sl sy_sl],'-','Color',[0.80 0.10 0.10],'LineWidth',1.4);
    plot(ax,D1b(1)+tl*[-sx_sl sx_sl],D1b(2)+tl*[-sy_sl sy_sl],'-','Color',[0.80 0.10 0.10],'LineWidth',1.4);
    plot(ax,[A(1) D1a(1)],[A(2) D1a(2)],'--','Color',[0.85 0.55 0.55],'LineWidth',0.6);
    plot(ax,[foot(1) D1b(1)],[foot(2) D1b(2)],'--','Color',[0.85 0.55 0.55],'LineWidth',0.6);
    mid = (D1a+D1b)/2+0.35*[nx;ny];
    text(ax,mid(1),mid(2),'4.572 mm','HorizontalAlignment','center','VerticalAlignment','bottom', ...
         'Color',[0.65 0.05 0.05],'FontSize',12,'FontWeight','bold','Rotation',tilt_deg+half_angle_deg);

    % dim 0.41 perpendicular
    plot(ax,foot(1)+tl*[-sx_sl sx_sl],foot(2)+tl*[-sy_sl sy_sl],'-','Color',[0.80 0.10 0.10],'LineWidth',1.4);
    plot(ax,sxs+tl*[-sx_sl sx_sl],sys+tl*[-sy_sl sy_sl],'-','Color',[0.80 0.10 0.10],'LineWidth',1.4);
    plot(ax,[foot(1) sxs],[foot(2) sys],'-','Color',[0.80 0.10 0.10],'LineWidth',1.4);
    lbl2 = [sxs;sys]+0.6*[nx;ny];
    text(ax,lbl2(1)-0.9*sx_sl,lbl2(2)-0.9*sy_sl,'0.41 mm','HorizontalAlignment','center', ...
         'VerticalAlignment','bottom','Color',[0.65 0.05 0.05],'FontSize',12,'FontWeight','bold', ...
         'Rotation',tilt_deg+half_angle_deg);

    % annotations
    plot(ax,A(1),A(2),'k.','MarkerSize',14);
    text(ax,A(1)-0.4*ux-0.7*px,A(2)-0.4*uy-0.7*py,'Apex (0, 0)','HorizontalAlignment','right', ...
         'VerticalAlignment','middle','FontSize',10,'Color',[0.2 0.2 0.2]);
    bs_lbl = [-3.0, 11.5];
    plot(ax,[bs_lbl(1) sxs],[bs_lbl(2) sys],'-','Color',[0.85 0.55 0.55],'LineWidth',0.6);
    text(ax,bs_lbl(1),bs_lbl(2)+0.5,'B_{surface}  (Hall sensor, in air)', ...
         'HorizontalAlignment','left','FontSize',10,'Color',[0.65 0.05 0.05]);
    mid_cone = 0.55*(A+Ju)-0.4*[nx;ny];
    text(ax,mid_cone(1),mid_cone(2),sprintf('Cone face (natural, %.2f\\circ)',half_angle_deg), ...
         'HorizontalAlignment','center','FontSize',9,'Color',[0.20 0.20 0.30],'Rotation',tilt_deg+half_angle_deg);
    mid_cyl = 0.5*(Ju+Tu)-0.5*[px;py];
    text(ax,mid_cyl(1),mid_cyl(2),'Cylinder body (R = 3.175 mm)','HorizontalAlignment','center', ...
         'FontSize',9,'Color',[0.20 0.20 0.30],'Rotation',tilt_deg);
    text(ax,axis_end(1)+0.3*ux,axis_end(2)+0.3*uy,sprintf('Pole axis (tilt %.2f\\circ)',tilt_deg), ...
         'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',9, ...
         'Color',[0.45 0.45 0.55],'Rotation',tilt_deg);
    text(ax,13,-3.4,['Upper pole — NOT milled (natural cone). 4.572 mm along the cone slant; ' ...
         '0.41 mm perpendicular to it.'],'HorizontalAlignment','center','FontSize',8.5, ...
         'Color',[0.35 0.35 0.45],'FontAngle','italic');

    axis(ax,'equal'); xlim(ax,[-4,30]); ylim(ax,[-5.5,24]);
    xlabel(ax,'x [mm]'); ylabel(ax,'y [mm]');
    title(ax,sprintf('%s — Upper pole (azimuth %d\\circ) — B_{surface} placement', p.name, p.az));
    set(ax,'FontSize',10,'Box','on','Layer','top'); grid(ax,'off');
    ax.Toolbar.Visible = 'off';
end
