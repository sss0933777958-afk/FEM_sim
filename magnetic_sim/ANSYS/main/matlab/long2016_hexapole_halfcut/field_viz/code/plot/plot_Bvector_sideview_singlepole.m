function plot_Bvector_sideview_singlepole(out_path)
% PLOT_BVECTOR_SIDEVIEW_SINGLEPOLE
%   Single-pole (filled cone) B field, y=0 side view (x-z, WP-ish frame).
%   RAW data: real FEM nodes (grid-sampled, NO interpolation) + raw physical B
%   (NO source flip). Variable-length arrows (length ~ |B|, clamped), colour = |B|.
%   Full-cone outline (black) + WP marker. NO title (per user).
%
%   out_path (optional): PNG output path. Default = field_viz/figures/singlepole_field.png

    if nargin < 1
        out_path = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
                    'long2016_hexapole_halfcut\field_viz\figures\singlepole_field.png'];
    end

    %% --- Style ---
    FONT_NAME = 'Helvetica';  FONT_LBL = 12;  FONT_CB = 11;  LW_STEEL = 2.0;  DPI = 300;

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');
    cnst = mt_constants();
    res_dir = ansys_path('long2016_hexapole_halfcut', 'data', 'singlepole');

    %% --- geometry constants (mm) ---
    beta      = atan2(3.0, 15.0);
    R_norm_xy = cnst.R_norm_xy * 1e3;
    R_norm_z  = cnst.R_norm_z  * 1e3;
    SPH_OFST  = cnst.SPH_OFST  * 1e3;
    POLE_R    = cnst.POLE_R    * 1e3;
    POLE_LEN  = cnst.POLE_CONE_LEN * 1e3;
    Z_CONE    = SPH_OFST - R_norm_z;

    %% --- load single-pole data; SOURCE convention (flip P1 lower pole -> B out of tip) ---
    fprintf('Loading single-pole data...\n');
    d = import_ansys_data(res_dir, 'all', 'singlepole');
    d.bx = -d.bx;  d.by = -d.by;  d.bz = -d.bz;     % source convention (P1 lower s=-1)
    fprintf('  Matched %d nodes, |B| max %.4g T (source: flipped)\n', numel(d.x), max(sqrt(d.bx.^2+d.by.^2+d.bz.^2)));

    %% --- plot window ---
    xlim_v = [-2, 12];
    zlim_v = [-16, -10];

    %% --- y=0 REAL-NODE grid sampling (steel + air, NO interpolation) ---
    xmm=d.x*1e3; ymm=d.y*1e3; zmm=d.z*1e3;
    inwin = xmm>=xlim_v(1)&xmm<=xlim_v(2)&zmm>=zlim_v(1)&zmm<=zlim_v(2)&abs(ymm)<0.8;
    xs=xmm(inwin); ys=ymm(inwin); zs=zmm(inwin);
    bxs=d.bx(inwin); bys=d.by(inwin); bzs=d.bz(inwin);
    nx=66; nz=40;
    xe=linspace(xlim_v(1),xlim_v(2),nx+1); ze=linspace(zlim_v(1),zlim_v(2),nz+1);
    ix=discretize(xs,xe); iz=discretize(zs,ze); ok=~isnan(ix)&~isnan(iz);
    cellid=(iz(ok)-1)*nx+ix(ok);
    xv=xs(ok);zv=zs(ok);bxv=bxs(ok);byv=bys(ok);bzv=bzs(ok);ayv=abs(ys(ok));
    [gc,~]=findgroups(cellid); gi=(1:numel(cellid))';
    pick=accumarray(gc,gi,[],@(r) r(find(ayv(r)==min(ayv(r)),1)));
    Xg=xv(pick);Zg=zv(pick);Bx_g=bxv(pick);By_g=byv(pick);Bz_g=bzv(pick);
    bsum_g=sqrt(Bx_g.^2+By_g.^2+Bz_g.^2);
    fprintf('  Real-node arrows: %d\n', numel(Xg));

    %% --- UNIFORM-length arrows (all same size); colour = |B| ---
    CAP = 0.045;                                  % T (colour cap)
    arrow_len = 0.9 * (xlim_v(2)-xlim_v(1))/nx;   % uniform arrow length (mm)
    Bip = hypot(Bx_g, Bz_g); Bip(Bip<eps)=eps;    % in-plane magnitude
    bx_q = Bx_g./Bip * arrow_len;                 % unit in-plane direction x len
    bz_q = Bz_g./Bip * arrow_len;
    mag_c = min(bsum_g, CAP);

    %% --- plot ---
    fig = figure('Position', [60 60 1300 620], 'Color', 'w');
    set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LBL);
    n_bins = 28; edges_b = linspace(0, CAP, n_bins+1); cmap_n = turbo(n_bins);
    hold on;
    for k = 1:n_bins
        in = mag_c >= edges_b(k) & mag_c < edges_b(k+1);
        if k==n_bins, in = in | (mag_c >= edges_b(end)); end
        if any(in)
            quiver(Xg(in), Zg(in), bx_q(in), bz_q(in), 0, ...
                   'Color', cmap_n(k,:), 'LineWidth', 1.0, 'MaxHeadSize', 0.6);
        end
    end

    %% --- full-cone outline + WP ---
    T = [R_norm_xy, Z_CONE]; Lsl = sqrt(POLE_LEN^2 + POLE_R^2);
    C_up = T + Lsl*[cos(beta),  sin(beta)];
    C_lo = T + Lsl*[cos(beta), -sin(beta)];
    plot([T(1),C_up(1)],[T(2),C_up(2)],'k-','LineWidth',LW_STEEL);
    plot([T(1),C_lo(1)],[T(2),C_lo(2)],'k-','LineWidth',LW_STEEL);
    plot([C_up(1),C_lo(1)],[C_up(2),C_lo(2)],'k-','LineWidth',LW_STEEL);
    plot(0, SPH_OFST, 'k+', 'MarkerSize', 16, 'LineWidth', 2.4);
    text(0+0.25, SPH_OFST+0.35, 'WP', 'FontSize', 13, 'FontWeight', 'bold');
    hold off;

    %% --- colorbar + axes (NO title) ---
    colormap(turbo); clim([0, CAP]);
    cb = colorbar; ylabel(cb, '|B|  [T]', 'FontSize', FONT_CB); set(cb,'FontSize',FONT_CB);
    axis equal; grid on; xlim(xlim_v); ylim(zlim_v);
    xlabel('x [mm]', 'FontSize', FONT_LBL); ylabel('z [mm]', 'FontSize', FONT_LBL);

    exportgraphics(fig, out_path, 'Resolution', DPI);
    fprintf('Saved: %s\n', out_path);
end
