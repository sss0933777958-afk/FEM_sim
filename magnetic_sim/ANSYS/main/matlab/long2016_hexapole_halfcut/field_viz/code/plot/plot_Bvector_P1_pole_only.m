function plot_Bvector_P1_pole_only()
% PLOT_BVECTOR_P1_POLE_ONLY
%   Side view (xz, y=0) of the ENTIRE P1 magnetic pole — quiver of FEM B field.
%
%   CORRECT geometry per APDL (MT_Sim_P1.txt lines 119-163):
%     P1 = revolved body around the axis (y=0, z=-13) from KP1=(R_norm_xy,0,-13)
%          to KP5=(X0_LOW-10mm, 0, -13) = (37.5, 0, -13).
%     2D generatrix in z=-13 plane (5 lines, closed):
%       KP1 (tip)         = (0.408, 0)
%       KP2 (tip fillet)  ≈ (0.45,  0.04)
%       KP3 (cone base)   = (15.408, 3)   POLE_R=3, POLE_CONE_LEN=15
%       KP4 (cyl end out) = (37.5,   3)
%       KP5 (cyl end ax)  = (37.5,   0)
%     VROTAT around KP1-KP5 axis -> body of revolution.
%     Halfcut: VSBV subtracts BLOCK z in [-13, -7] -> only z <= -13 remains.
%
%   So in xz at y=0, the half-body outline is:
%       tip (0.408, -13)
%        \ cone slope (R(x) = 0 .. 3)
%         \-> (15.408, -16)              (cone base bottom)
%          --- (cylinder bottom, R=3) --> (37.5, -16)
%                                          | (cylinder end face)
%                                          (37.5, -13)
%        <--- halfcut plane z=-13 <------- back to tip
%
%   Total pole length: 37.5 - 0.408 = 37.092 mm.
%
%   No post / yoke / coil / sensor / connector blocks shown.
%   Data: coil1 smrt 4 (P1 excitation, sign-corrected).

    %% --- Style ---
    FONT_NAME = 'Helvetica';
    FONT_LBL  = 12;
    FONT_TTL  = 13;
    FONT_CB   = 11;
    LW_STEEL  = 2.0;
    DPI       = 300;

    %% --- Paths ---
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\hung\analysis\core');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\hexapole-long2016\analysis');

    cnst = mt_constants();

    res_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\coil1\standard';
    out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    %% --- Geometry (APDL mm) ---
    R_norm_xy  = cnst.R_norm_xy * 1e3;          % 0.408
    R_norm_z   = cnst.R_norm_z  * 1e3;          % 0.289
    SPH_OFST   = cnst.SPH_OFST  * 1e3;          % -12.71
    POLE_R     = cnst.POLE_R    * 1e3;          % 3
    POLE_TIP_R = cnst.POLE_TIP_R * 1e3;         % 0.04
    POLE_LEN   = cnst.POLE_CONE_LEN * 1e3;      % 15
    Z_CONE     = SPH_OFST - R_norm_z;           % -13  (revolution axis)

    x_tip   = R_norm_xy;                        %  0.408  KP1
    x_base  = x_tip + POLE_LEN;                 % 15.408  KP3 (cone base, R=POLE_R)
    X0_LOW  = (42 + 53)/2;                      % 47.5    YOKE_MID_R
    x_end   = X0_LOW - 10;                      % 37.5    KP5 (cylinder end on axis)

    pole_total_len = x_end - x_tip;             % 37.092 mm  (matches user's 37.09)

    %% --- Load coil1 (P1 excitation) ---
    fprintf('Loading coil1 (P1 excitation) smrt 4 data...\n');
    d = import_ansys_data(res_dir, 'all', 'coil1');

    %% --- Per-pole sign correction ---
    bbar_mat = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\bs_matrix\Bbar_S_4p572.mat';
    sign_p1 = +1;
    if exist(bbar_mat, 'file')
        Bs = load(bbar_mat, 'col_sign');
        sign_p1 = Bs.col_sign(1);
    end
    fprintf('P1 sign correction: %+d\n', sign_p1);
    if sign_p1 < 0
        d.bx = -d.bx; d.by = -d.by; d.bz = -d.bz;
    end

    %% --- 3D scatteredInterpolant ---
    fprintf('Building 3D scatteredInterpolant from %d nodes...\n', length(d.x));
    F_bx = scatteredInterpolant(d.x, d.y, d.z, d.bx, 'linear', 'none');
    F_by = scatteredInterpolant(d.x, d.y, d.z, d.by, 'linear', 'none');
    F_bz = scatteredInterpolant(d.x, d.y, d.z, d.bz, 'linear', 'none');

    %% --- xz grid (entire pole: cone + cylindrical extension, halfcut lower half) ---
    xlim_v = [-3, 40];
    zlim_v = [-17, -10];
    grid_x = 120;
    grid_z = 28;
    x_centers = linspace(xlim_v(1)+0.15, xlim_v(2)-0.15, grid_x);
    z_centers = linspace(zlim_v(1)+0.15, zlim_v(2)-0.15, grid_z);
    [Xg, Zg] = meshgrid(x_centers, z_centers);
    Xg = Xg(:); Zg = Zg(:);

    x_apdl = Xg * 1e-3;
    y_apdl = zeros(size(Xg));
    z_apdl = Zg * 1e-3;

    Bx_g = F_bx(x_apdl, y_apdl, z_apdl);
    By_g = F_by(x_apdl, y_apdl, z_apdl);
    Bz_g = F_bz(x_apdl, y_apdl, z_apdl);

    keep = ~isnan(Bx_g) & ~isnan(Bz_g);
    Xg = Xg(keep); Zg = Zg(keep);
    Bx_g = Bx_g(keep); By_g = By_g(keep); Bz_g = Bz_g(keep);

    bsum_g = sqrt(Bx_g.^2 + By_g.^2 + Bz_g.^2);
    fprintf('  Grid arrows (after NaN filter): %d\n', length(Xg));
    fprintf('  |B| range:  [%.3g, %.3g] T\n', min(bsum_g), max(bsum_g));

    %% --- Color clipping (so arm-region 0.1 T arrows are visible
    %       without cone peak (0.5 T) saturating the colormap) ---
    b_clip = quantile(bsum_g, 0.90);   % 90th percentile
    fprintf('  Color clip at 90th pct = %.3g T\n', b_clip);

    %% --- Arrow scaling: per-arrow normalized length (uniform visual size)
    %       so weak arrows in arm don't shrink to invisible dots.
    %       Color encodes magnitude; length is uniform. ---
    dx_grid = (xlim_v(2)-xlim_v(1)) / grid_x;
    arrow_max = 0.85 * dx_grid;
    b_norm = bsum_g;
    b_norm(b_norm < eps) = eps;
    scale_per = arrow_max ./ b_norm;
    bx_q = Bx_g .* scale_per;
    bz_q = Bz_g .* scale_per;

    %% --- Plot ---
    fig = figure('Position', [50 50 1500 700], 'Color', 'w');
    set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LBL);

    n_bins  = 24;
    edges_b = linspace(0, b_clip, n_bins + 1);
    cmap_n  = turbo(n_bins);
    lw_rng  = [0.55, 1.90];

    hold on;
    for k = 1:n_bins
        if k < n_bins
            in_bin = bsum_g >= edges_b(k) & bsum_g < edges_b(k+1);
        else
            in_bin = bsum_g >= edges_b(k);   % top bin captures clipped tail
        end
        if any(in_bin)
            lw = lw_rng(1) + (k-1)/(n_bins-1) * (lw_rng(2)-lw_rng(1));
            quiver(Xg(in_bin), Zg(in_bin), bx_q(in_bin), bz_q(in_bin), 0, ...
                   'Color', cmap_n(k,:), 'LineWidth', lw, 'MaxHeadSize', 0.6);
        end
    end

    %% --- P1 pole OUTLINE (revolved cone + cylinder, halfcut lower half) ---
    %   Boundary in xz at y=0 (closed loop):
    %     (x_tip,  Z_CONE)                       tip on axis
    %     -> (x_base, Z_CONE - POLE_R)           cone slope down
    %     -> (x_end,  Z_CONE - POLE_R)           cylinder bottom
    %     -> (x_end,  Z_CONE)                    cylinder right end face
    %     -> back to (x_tip, Z_CONE)             halfcut plane (top)
    plot([x_tip,  x_base], [Z_CONE,            Z_CONE - POLE_R],   'k-', 'LineWidth', LW_STEEL);  % cone slope
    plot([x_base, x_end ], [Z_CONE - POLE_R,   Z_CONE - POLE_R],   'k-', 'LineWidth', LW_STEEL);  % cylinder bottom
    plot([x_end,  x_end ], [Z_CONE - POLE_R,   Z_CONE],            'k-', 'LineWidth', LW_STEEL);  % cylinder end face
    plot([x_end,  x_tip ], [Z_CONE,            Z_CONE],            'k-', 'LineWidth', LW_STEEL);  % halfcut plane

    %% --- WP marker + region labels ---
    plot(0, SPH_OFST, 'k+', 'MarkerSize', 16, 'LineWidth', 2.4);
    text(0 + 0.3, SPH_OFST + 0.35, 'WP', 'FontSize', 13, 'FontWeight', 'bold');

    text(x_tip + 0.4, Z_CONE - 0.55, 'P1 tip', 'FontSize', 12, 'FontWeight', 'bold');
    text((x_tip + x_base)/2, Z_CONE - POLE_R - 0.55, 'cone (R: 0.04 -> 3 mm)', ...
         'FontSize', 11, 'HorizontalAlignment', 'center');
    text((x_base + x_end)/2, Z_CONE - POLE_R - 0.55, 'cylinder (R = 3 mm)', ...
         'FontSize', 11, 'HorizontalAlignment', 'center');
    text(x_end + 0.4, Z_CONE - POLE_R/2, 'pole end', ...
         'FontSize', 11, 'HorizontalAlignment', 'left');

    % Halfcut plane label
    text(38, Z_CONE + 0.25, 'halfcut plane (z = -13)', ...
         'FontSize', 10, 'Color', [0.45 0.45 0.45], 'FontAngle', 'italic', ...
         'HorizontalAlignment', 'right');

    % Total pole length annotation
    y_anno = Z_CONE - POLE_R - 1.4;
    plot([x_tip, x_end], [y_anno, y_anno], 'k-', 'LineWidth', 1.0);
    plot([x_tip, x_tip], [y_anno-0.2, y_anno+0.2], 'k-', 'LineWidth', 1.0);
    plot([x_end, x_end], [y_anno-0.2, y_anno+0.2], 'k-', 'LineWidth', 1.0);
    text((x_tip + x_end)/2, y_anno - 0.55, ...
         sprintf('total pole length = %.2f mm', pole_total_len), ...
         'FontSize', 11, 'HorizontalAlignment', 'center', 'FontAngle', 'italic');

    hold off;

    %% --- Colorbar + axes ---
    colormap(turbo);
    clim([0, b_clip]);
    cb = colorbar;
    ylabel(cb, sprintf('|B|  [Tesla]  (clipped at 90th pct = %.2f T)', b_clip), ...
        'FontSize', FONT_CB);
    set(cb, 'FontSize', FONT_CB);

    axis equal;
    grid on;
    xlim(xlim_v); ylim(zlim_v);
    set(gca, 'GridAlpha', 0.18, 'Layer', 'top');
    xlabel('x [mm]  (APDL frame, +x = P1 cone axis)', 'FontSize', FONT_LBL);
    ylabel('z [mm]  (APDL frame)', 'FontSize', FONT_LBL);
    title(['P1 magnetic pole (revolved cone + cylinder, full 37.09 mm) — ' ...
           'B-field side view xz at y=0 — halfcut lower half (z<=-13), smrt 4 mesh'], ...
          'FontSize', FONT_TTL, 'Interpreter', 'none');

    %% --- Save ---
    out_path = fullfile(out_dir, 'P1_Bvector_pole_only.png');
    exportgraphics(fig, out_path, 'Resolution', DPI);
    fprintf('\nSaved: %s\n', out_path);
end
