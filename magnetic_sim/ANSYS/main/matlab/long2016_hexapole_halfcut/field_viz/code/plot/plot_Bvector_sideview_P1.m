function plot_Bvector_sideview_P1()
% PLOT_BVECTOR_SIDEVIEW_P1
%   Side view of P1-excitation B field in xz plane (y=0).
%   P1 is azimuth 0° lower pole — horizontal, halfcut (D-shape).
%   Cone OUTLINE (black) shows half-cone: tip → lower slant → block → milled flat → tip.
%   B_tip (blue marker) at pole tip, B_surface (red sensor) on milled flat.

    %% --- Style ---
    FONT_NAME = 'Helvetica';
    FONT_LBL  = 12;
    FONT_TTL  = 13;
    FONT_CB   = 11;
    LW_STEEL  = 2.0;
    DPI       = 300;

    %% --- Paths ---
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');

    cnst = mt_constants();

    res_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data\coil1\standard';
    out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    %% --- Geometry constants (mm) ---
    beta       = atan2(3.0, 15.0);              % cone half-angle, 11.31°
    R_norm_xy  = cnst.R_norm_xy * 1e3;          % 0.408 mm
    R_norm_z   = cnst.R_norm_z  * 1e3;          % 0.289 mm
    SPH_OFST   = cnst.SPH_OFST  * 1e3;          % -12.71 mm
    POLE_R     = cnst.POLE_R    * 1e3;          % 3 mm
    POLE_LEN   = cnst.POLE_CONE_LEN * 1e3;      % 15 mm

    % P1 cone axis sits at z = SPH_OFST - R_norm_z = -13 mm
    Z_CONE = SPH_OFST - R_norm_z;               % -13 mm (milled flat z)

    fprintf('Z_CONE (milled flat) = %.3f mm\n', Z_CONE);

    %% --- Load coil1 (= P1) FEM data ---
    fprintf('\nLoading coil1 (P1 excitation) smrt 4 data...\n');
    d = import_ansys_data(res_dir, 'all', 'coil1');

    %% --- Per-pole sign correction (read from Bbar_S mat) ---
    bbar_mat = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\bs_matrix\data\Bbar_S_4p572.mat';
    sign_p1 = +1;
    if exist(bbar_mat, 'file')
        Bs = load(bbar_mat, 'col_sign');
        sign_p1 = Bs.col_sign(1);   % P1 is paper index 1
    end
    fprintf('P1 sign correction factor: %+d\n', sign_p1);
    if sign_p1 < 0
        d.bx = -d.bx; d.by = -d.by; d.bz = -d.bz;
        fprintf('  → Negated B field (display flux yoke→tip→WP)\n');
    end

    %% --- Build 3D scatteredInterpolant ---
    fprintf('Building 3D scatteredInterpolant from %d nodes...\n', length(d.x));
    F_bx = scatteredInterpolant(d.x, d.y, d.z, d.bx, 'linear', 'none');
    F_by = scatteredInterpolant(d.x, d.y, d.z, d.by, 'linear', 'none');
    F_bz = scatteredInterpolant(d.x, d.y, d.z, d.bz, 'linear', 'none');

    %% --- Regular xz grid in y=0 plane ---
    xlim_v = [-3, 12];
    zlim_v = [-15.5, -9];
    grid_x = 60; grid_z = 38;
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

    fprintf('  Regular-grid arrows (after NaN filter): %d\n', length(Xg));

    %% Arrow scaling
    arrow_max = 0.45;
    scale = arrow_max / max(bsum_g);
    bx_q = Bx_g * scale;
    bz_q = Bz_g * scale;

    %% --- Plot ---
    fig = figure('Position', [50 50 1400 800], 'Color', 'w');
    set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LBL);

    n_bins  = 28;
    edges_b = linspace(0, max(bsum_g), n_bins + 1);
    cmap_n  = turbo(n_bins);
    lw_rng  = [0.70, 2.40];

    hold on;
    for k = 1:n_bins
        in_bin = bsum_g >= edges_b(k) & bsum_g < edges_b(k+1);
        if k == n_bins
            in_bin = in_bin | (bsum_g >= edges_b(end));
        end
        if any(in_bin)
            lw = lw_rng(1) + (k-1)/(n_bins-1) * (lw_rng(2)-lw_rng(1));
            quiver(Xg(in_bin), Zg(in_bin), bx_q(in_bin), bz_q(in_bin), 0, ...
                   'Color', cmap_n(k,:), 'LineWidth', lw, 'MaxHeadSize', 0.6);
        end
    end

    %% --- P1 halfcut D-shape OUTLINE ---
    % Tip at (R_norm_xy, Z_CONE) — i.e. (0.408, -13)
    % Lower cone slant: tip down to (R_norm_xy + POLE_LEN, Z_CONE - POLE_R)
    % Block end (right side): vertical line at x = R_norm_xy + POLE_LEN + 8
    T   = [R_norm_xy, Z_CONE];
    C_b = [R_norm_xy + POLE_LEN, Z_CONE - POLE_R];   % cone bottom-right corner
    x_block_end = R_norm_xy + POLE_LEN + 8;
    B_b = [x_block_end, Z_CONE - POLE_R];            % block bottom-right
    B_t = [x_block_end, Z_CONE];                     % block top-right

    % Lower slant
    plot([T(1), C_b(1)], [T(2), C_b(2)], 'k-', 'LineWidth', LW_STEEL);
    % Bottom of cone+block
    plot([C_b(1), B_b(1)], [C_b(2), B_b(2)], 'k-', 'LineWidth', LW_STEEL);
    % Block right edge
    plot([B_b(1), B_t(1)], [B_b(2), B_t(2)], 'k-', 'LineWidth', LW_STEEL);
    % Milled flat (top of D-shape) — from tip to block right top
    plot([T(1), B_t(1)], [T(2), B_t(2)], 'k-', 'LineWidth', LW_STEEL);

    %% --- WP marker ---
    plot(0, SPH_OFST, 'k+', 'MarkerSize', 16, 'LineWidth', 2.4);
    text(0 + 0.3, SPH_OFST + 0.35, 'WP', 'FontSize', 13, 'FontWeight', 'bold');

    %% --- B_tip marker (at pole tip, blue square) ---
    plot(T(1), T(2), 's', 'MarkerSize', 14, ...
         'MarkerFaceColor', [0.10 0.40 0.85], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    text(T(1) - 0.6, T(2) - 0.55, 'B_{tip}', ...
         'FontSize', 14, 'Color', [0.10 0.40 0.85], 'FontWeight', 'bold', ...
         'Interpreter', 'tex');

    %% --- B_surface sensor (red square ON milled flat) ---
    S_ALONG  = 4.572;   % mm along milled flat from tip
    S_OFFSET = 0.41;    % mm perpendicular (up from milled flat)
    sens = [T(1) + S_ALONG, T(2) + S_OFFSET];

    panel_L = 1.0;
    panel_W = 0.20;
    hl = panel_L/2;
    hw = panel_W/2;
    % Sensor panel is HORIZONTAL (aligned with milled flat)
    P1 = sens + [+hl, +hw];
    P2 = sens + [-hl, +hw];
    P3 = sens + [-hl, -hw];
    P4 = sens + [+hl, -hw];
    fill([P1(1) P2(1) P3(1) P4(1)], [P1(2) P2(2) P3(2) P4(2)], ...
         [0.92 0.15 0.15], 'EdgeColor', 'k', 'LineWidth', 1.5);

    % Sensor n+ arrow (perpendicular to milled flat, pointing +z = up)
    arrow_n = 0.6;
    n_tip = sens + [0, arrow_n];
    plot([sens(1), n_tip(1)], [sens(2), n_tip(2)], '-', ...
         'Color', [0.92 0.15 0.15], 'LineWidth', 2.4);
    plot(n_tip(1), n_tip(2), '^', 'MarkerSize', 8, ...
         'MarkerFaceColor', [0.92 0.15 0.15], 'MarkerEdgeColor', [0.92 0.15 0.15]);

    text(sens(1) + 0.4, sens(2) + 0.6, 'B_{surface}', ...
         'FontSize', 14, 'Color', [0.92 0.15 0.15], 'FontWeight', 'bold');

    hold off;

    %% --- Colorbar + axes ---
    colormap(turbo);
    clim([0, max(bsum_g)]);
    cb = colorbar;
    ylabel(cb, '|B|  [Tesla]  (3D magnitude)', 'FontSize', FONT_CB);
    set(cb, 'FontSize', FONT_CB);

    axis equal;
    grid on;
    xlim(xlim_v); ylim(zlim_v);
    xlabel('x [mm]  (APDL frame, +x = P1 cone axis)', 'FontSize', FONT_LBL);
    ylabel('z [mm]  (APDL frame)', 'FontSize', FONT_LBL);
    title(['P1 excitation B-field (sign-corrected: flux yoke→tip→WP) — ' ...
           'lower pole halfcut (D-shape, milled flat at z = -13 mm), smrt 4 mesh'], ...
          'FontSize', FONT_TTL, 'Interpreter', 'none');

    out_path = fullfile(out_dir, 'Bvector_sideview_P1.png');
    exportgraphics(fig, out_path, 'Resolution', DPI);
    fprintf('\nSaved: %s\n', out_path);
end
