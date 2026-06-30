function plot_Bvector_sideview_singlepole_interp(out_path)
% PLOT_BVECTOR_SIDEVIEW_SINGLEPOLE_INTERP
%   Single-pole (filled cone) B field, y=0 side view (x-z, WP-ish frame).
%   Same as plot_Bvector_sideview_singlepole, BUT the empty grid cells that
%   fall INSIDE the cone (steel / pole body) are FILLED by 2-D interpolation
%   so the arrow density inside the cone matches the rest of the plot.
%
%   *** INTERPOLATED FIGURE ***  Arrows inside the cone that had no real FEM
%   node in their grid cell are scatteredInterpolant (linear) values built
%   ONLY from cone-interior near-plane nodes (no air smearing). Everything
%   outside the cone stays raw FEM nodes (no interpolation). Plot window is
%   extended to x = 15 mm to show the full cone base.
%
%   out_path (optional): PNG output path.
%       Default = field_viz/figures/singlepole_field_interp.png

    if nargin < 1
        out_path = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
                    'long2016_hexapole_halfcut\field_viz\figures\singlepole_field_interp.png'];
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

    %% --- plot window (x extended to 15 mm to show full cone base) ---
    xlim_v = [-2, 15];
    zlim_v = [-16, -10];

    %% --- cone outline (full cone) + inside-cone test ---
    T = [R_norm_xy, Z_CONE]; Lsl = sqrt(POLE_LEN^2 + POLE_R^2);
    C_up = T + Lsl*[cos(beta),  sin(beta)];
    C_lo = T + Lsl*[cos(beta), -sin(beta)];
    cone_x0 = T(1);  cone_x1 = C_up(1);  cone_zc = T(2);   % vertex x, base x, axis z
    incone  = @(x,z) (x >= cone_x0) & (x <= cone_x1) & ...
                     (abs(z - cone_zc) <= POLE_R .* (x - cone_x0) ./ (cone_x1 - cone_x0));

    %% --- y=0 REAL-NODE grid sampling (steel + air, NO interpolation) ---
    xmm=d.x*1e3; ymm=d.y*1e3; zmm=d.z*1e3;
    inwin = xmm>=xlim_v(1)&xmm<=xlim_v(2)&zmm>=zlim_v(1)&zmm<=zlim_v(2)&abs(ymm)<0.8;
    xs=xmm(inwin); ys=ymm(inwin); zs=zmm(inwin);
    bxs=d.bx(inwin); bys=d.by(inwin); bzs=d.bz(inwin);
    % keep the same cell SIZE as the raw figure (66 cells over old 14 mm width)
    cellw = (12-(-2))/66;                                   % ~0.2121 mm
    nx = round((xlim_v(2)-xlim_v(1))/cellw);  nz = 40;      % nx ~ 80
    xe=linspace(xlim_v(1),xlim_v(2),nx+1); ze=linspace(zlim_v(1),zlim_v(2),nz+1);
    ix=discretize(xs,xe); iz=discretize(zs,ze); ok=~isnan(ix)&~isnan(iz);
    cellid=(iz(ok)-1)*nx+ix(ok);
    xv=xs(ok);zv=zs(ok);bxv=bxs(ok);byv=bys(ok);bzv=bzs(ok);ayv=abs(ys(ok));
    [gc,gid]=findgroups(cellid); gi=(1:numel(cellid))';
    pick=accumarray(gc,gi,[],@(r) r(find(ayv(r)==min(ayv(r)),1)));
    Xg=xv(pick);Zg=zv(pick);Bx_g=bxv(pick);By_g=byv(pick);Bz_g=bzv(pick);
    is_interp = false(numel(Xg),1);                         % raw nodes
    occ_id = gid;                                           % occupied cell ids
    fprintf('  Real-node arrows: %d\n', numel(Xg));

    %% --- INTERPOLATION: fill empty cells whose center is inside the cone ---
    in_cone_node = incone(xs, zs);                          % cone-interior near-plane nodes
    Fbx = scatteredInterpolant(xs(in_cone_node), zs(in_cone_node), bxs(in_cone_node), 'linear','none');
    Fby = scatteredInterpolant(xs(in_cone_node), zs(in_cone_node), bys(in_cone_node), 'linear','none');
    Fbz = scatteredInterpolant(xs(in_cone_node), zs(in_cone_node), bzs(in_cone_node), 'linear','none');
    xc = (xe(1:end-1)+xe(2:end))/2;  zc = (ze(1:end-1)+ze(2:end))/2;
    nfill = 0;
    for jz = 1:nz
        for jx = 1:nx
            id = (jz-1)*nx + jx;
            if ismember(id, occ_id), continue; end          % cell already has a real node
            cx = xc(jx);  cz = zc(jz);
            if ~incone(cx, cz), continue; end               % only inside the cone
            bxi = Fbx(cx,cz);  byi = Fby(cx,cz);  bzi = Fbz(cx,cz);
            if isnan(bxi)||isnan(bzi), continue; end         % outside interpolant hull
            Xg(end+1,1)=cx; Zg(end+1,1)=cz;                  %#ok<AGROW>
            Bx_g(end+1,1)=bxi; By_g(end+1,1)=byi; Bz_g(end+1,1)=bzi; %#ok<AGROW>
            is_interp(end+1,1)=true;                          %#ok<AGROW>
            nfill = nfill + 1;
        end
    end
    fprintf('  Interpolated fill arrows (cone interior): %d\n', nfill);
    bsum_g=sqrt(Bx_g.^2+By_g.^2+Bz_g.^2);

    %% --- UNIFORM-length arrows (all same size); colour = |B| ---
    CAP = 0.045;                                  % T (colour cap)
    arrow_len = 0.9 * (xlim_v(2)-xlim_v(1))/nx;   % uniform arrow length (mm)
    Bip = hypot(Bx_g, Bz_g); Bip(Bip<eps)=eps;    % in-plane magnitude
    bx_q = Bx_g./Bip * arrow_len;                 % unit in-plane direction x len
    bz_q = Bz_g./Bip * arrow_len;
    mag_c = min(bsum_g, CAP);

    %% --- plot ---
    fig = figure('Position', [60 60 1380 620], 'Color', 'w');
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
    plot([T(1),C_up(1)],[T(2),C_up(2)],'k-','LineWidth',LW_STEEL);
    plot([T(1),C_lo(1)],[T(2),C_lo(2)],'k-','LineWidth',LW_STEEL);
    plot([C_up(1),C_lo(1)],[C_up(2),C_lo(2)],'k-','LineWidth',LW_STEEL);
    plot(0, SPH_OFST, 'k+', 'MarkerSize', 16, 'LineWidth', 2.4);
    text(0+0.25, SPH_OFST+0.35, 'WP', 'FontSize', 13, 'FontWeight', 'bold');
    % [INTERP] label so the figure itself flags the interpolated region
    text(xlim_v(2)-0.3, zlim_v(1)+0.45, 'cone interior: interpolated', ...
         'FontSize', 10, 'FontAngle','italic', 'Color',[0.25 0.25 0.25], ...
         'HorizontalAlignment','right');
    hold off;

    %% --- colorbar + axes (NO title) ---
    colormap(turbo); clim([0, CAP]);
    cb = colorbar; ylabel(cb, '|B|  [T]', 'FontSize', FONT_CB); set(cb,'FontSize',FONT_CB);
    axis equal; grid on; xlim(xlim_v); ylim(zlim_v);
    xlabel('x [mm]', 'FontSize', FONT_LBL); ylabel('z [mm]', 'FontSize', FONT_LBL);

    exportgraphics(fig, out_path, 'Resolution', DPI);
    fprintf('Saved: %s\n', out_path);
end
