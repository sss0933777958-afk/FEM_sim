function plot_Bvector_P1_full_circuit()
% PLOT_BVECTOR_P1_FULL_CIRCUIT
%   Full P1 magnetic structure side view (xz, y=0) — quiver of FEM B field
%   covering the ENTIRE magnetic path:
%     cone (x=0.4-15.4) -> horizontal arm (x=0-37.5) -> post (x=42.5-52.5)
%     -> yoke (z=0-2) -> coil ring (around post).
%
%   Same turbo + line-width scaling style as plot_Bvector_sideview_P1.m,
%   but the x range is now -3 ~ 60 mm to capture the WHOLE magnetic circuit.
%   No sensor markers — just the magnetic structure.
%
%   Data: coil1 smrt 4 baseline (P1 excitation, sign-corrected).

    %% --- Style ---
    FONT_NAME = 'Helvetica';
    FONT_LBL  = 12;
    FONT_TTL  = 13;
    FONT_CB   = 11;
    LW_STEEL  = 1.8;
    LW_HALFCUT = 0.9;
    DPI       = 300;

    COL_STEEL  = [0.78 0.78 0.82];
    COL_COIL   = [0.93 0.60 0.10];
    COL_HALFCUT = [0.45 0.45 0.45];

    %% --- Paths ---
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');

    cnst = mt_constants();

    res_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\coil1\standard';
    out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    %% --- Geometry (APDL mm) ---
    R_norm_xy = cnst.R_norm_xy * 1e3;            % 0.408
    R_norm_z  = cnst.R_norm_z  * 1e3;            % 0.289
    SPH_OFST  = cnst.SPH_OFST  * 1e3;            % -12.71
    POLE_R    = cnst.POLE_R    * 1e3;            % 3
    POLE_TIP_R = cnst.POLE_TIP_R * 1e3;          % 0.04
    POLE_LEN  = cnst.POLE_CONE_LEN * 1e3;        % 15

    YOKE_IN_R  = 42;
    YOKE_OUT_R = 53;
    YOKE_MID_R = (YOKE_IN_R + YOKE_OUT_R)/2;     % 47.5
    YOKE_H     = 2;
    PROT_R     = 5;
    PROT_H     = 7;
    X0_LOW     = YOKE_MID_R;                     % 47.5
    POST_X1    = X0_LOW - PROT_R;                % 42.5
    POST_X2    = X0_LOW + PROT_R;                % 52.5

    ARM_X1 = 0;
    ARM_X2 = X0_LOW - 10;                        % 37.5
    ARM_Z1 = -PROT_H - 6;                        % -13
    ARM_Z2 = -PROT_H;                            % -7
    Z_CONE = ARM_Z1;                             % -13

    x_tip  = R_norm_xy;
    x_base = x_tip + POLE_LEN;                   % 15.408

    COIL_IN_R   = 5;
    COIL_OUT_R  = 8;
    COIL_Z1     = -PROT_H;
    COIL_Z2     = 0;

    %% --- Load coil1 (P1 excitation) ---
    fprintf('Loading coil1 (P1 excitation) smrt 4 data...\n');
    d = import_ansys_data(res_dir, 'all', 'coil1');

    %% --- Per-pole sign correction (read from Bbar_S mat) ---
    bbar_mat = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\bs_matrix\Bbar_S_4p572.mat';
    sign_p1 = +1;
    if exist(bbar_mat, 'file')
        Bs = load(bbar_mat, 'col_sign');
        sign_p1 = Bs.col_sign(1);
    end
    fprintf('P1 sign correction: %+d\n', sign_p1);
    if sign_p1 < 0
        d.bx = -d.bx; d.by = -d.by; d.bz = -d.bz;
        fprintf('  -> B field negated (flux yoke->tip->WP)\n');
    end

    %% --- 3D scatteredInterpolant ---
    fprintf('Building 3D scatteredInterpolant from %d nodes...\n', length(d.x));
    F_bx = scatteredInterpolant(d.x, d.y, d.z, d.bx, 'linear', 'none');
    F_bz = scatteredInterpolant(d.x, d.y, d.z, d.bz, 'linear', 'none');
    F_by = scatteredInterpolant(d.x, d.y, d.z, d.by, 'linear', 'none');

    %% --- xz grid (full circuit view) ---
    xlim_v = [-3, 60];
    zlim_v = [-16.5, 4];
    grid_x = 105;
    grid_z = 34;
    x_centers = linspace(xlim_v(1)+0.3, xlim_v(2)-0.3, grid_x);
    z_centers = linspace(zlim_v(1)+0.3, zlim_v(2)-0.3, grid_z);
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
    bxz_g  = sqrt(Bx_g.^2 + Bz_g.^2);

    fprintf('  Grid arrows (after NaN filter): %d\n', length(Xg));
    fprintf('  |B_xz| range:  [%.3g, %.3g] T\n', min(bxz_g), max(bxz_g));

    %% --- Clip color scale (so faint air-arrows still visible, steel doesn't dominate) ---
    bxz_clip_max = quantile(bxz_g, 0.97);   % 97th percentile to avoid hot-tip outlier
    fprintf('  Color clip at 97th pct = %.3g T\n', bxz_clip_max);

    %% --- Arrow scaling (relative to grid spacing) ---
    dx_grid = x_centers(2) - x_centers(1);
    arrow_max = 0.85 * dx_grid;
    bxz_norm = bxz_g;
    bxz_norm(bxz_norm < eps) = eps;
    scale = arrow_max ./ bxz_norm;
    bx_q = Bx_g .* scale;
    bz_q = Bz_g .* scale;

    %% --- Plot ---
    fig = figure('Position', [30 30 1700 650], 'Color', 'w');
    set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LBL);

    n_bins  = 24;
    edges_b = linspace(0, bxz_clip_max, n_bins + 1);
    cmap_n  = turbo(n_bins);
    lw_rng  = [0.55, 2.10];

    hold on;
    for k = 1:n_bins
        if k < n_bins
            in_bin = bxz_g >= edges_b(k) & bxz_g < edges_b(k+1);
        else
            in_bin = bxz_g >= edges_b(k);   % top bin captures clipped tail
        end
        if any(in_bin)
            lw = lw_rng(1) + (k-1)/(n_bins-1) * (lw_rng(2)-lw_rng(1));
            quiver(Xg(in_bin), Zg(in_bin), bx_q(in_bin), bz_q(in_bin), 0, ...
                   'Color', cmap_n(k,:), 'LineWidth', lw, 'MaxHeadSize', 0.6);
        end
    end

    %% --- Draw STEEL outlines ---
    % Horizontal arm
    arm_x = [ARM_X1, ARM_X2, ARM_X2, ARM_X1, ARM_X1];
    arm_z = [ARM_Z1, ARM_Z1, ARM_Z2, ARM_Z2, ARM_Z1];
    plot(arm_x, arm_z, 'k-', 'LineWidth', LW_STEEL);

    % Cone (halfcut, lower half)
    cone_x = [x_tip, x_base, x_base, x_tip, x_tip];
    cone_z = [Z_CONE, Z_CONE, Z_CONE-POLE_R, Z_CONE-POLE_TIP_R, Z_CONE];
    plot(cone_x, cone_z, 'k-', 'LineWidth', LW_STEEL);

    % Post
    post_x = [POST_X1, POST_X2, POST_X2, POST_X1, POST_X1];
    post_z = [ARM_Z2,  ARM_Z2,  0,       0,       ARM_Z2];
    plot(post_x, post_z, 'k-', 'LineWidth', LW_STEEL);

    % Yoke (right part — from post outward to YOKE_OUT_R)
    yoke_x = [POST_X2-2, YOKE_OUT_R, YOKE_OUT_R, POST_X2-2, POST_X2-2];
    yoke_z = [0, 0, YOKE_H, YOKE_H, 0];
    plot(yoke_x, yoke_z, 'k-', 'LineWidth', LW_STEEL);

    % Yoke (left part — abstracted, continues toward other poles)
    yoke_left_x = [-3, POST_X1+2, POST_X1+2, -3];
    yoke_left_z = [0, 0, YOKE_H, YOKE_H];
    plot(yoke_left_x, yoke_left_z, 'k--', 'LineWidth', LW_STEEL*0.7);

    %% --- Coil ring cross-section (orange filled rectangles) ---
    coil_L_x = [X0_LOW-COIL_OUT_R, X0_LOW-COIL_IN_R, X0_LOW-COIL_IN_R, X0_LOW-COIL_OUT_R, X0_LOW-COIL_OUT_R];
    coil_R_x = [X0_LOW+COIL_IN_R,  X0_LOW+COIL_OUT_R, X0_LOW+COIL_OUT_R, X0_LOW+COIL_IN_R,  X0_LOW+COIL_IN_R];
    coil_z   = [COIL_Z1, COIL_Z1, COIL_Z2, COIL_Z2, COIL_Z1];
    fill(coil_L_x, coil_z, COL_COIL, 'EdgeColor', 'k', 'LineWidth', 1.2, ...
         'FaceAlpha', 0.55);
    fill(coil_R_x, coil_z, COL_COIL, 'EdgeColor', 'k', 'LineWidth', 1.2, ...
         'FaceAlpha', 0.55);
    text(X0_LOW-(COIL_IN_R+COIL_OUT_R)/2, (COIL_Z1+COIL_Z2)/2, char(8857), ...
         'FontSize', 22, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(X0_LOW+(COIL_IN_R+COIL_OUT_R)/2, (COIL_Z1+COIL_Z2)/2, char(8855), ...
         'FontSize', 22, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(X0_LOW, COIL_Z1-1.2, 'SOURC36 coil', ...
         'FontSize', 11, 'HorizontalAlignment', 'center', ...
         'Color', COL_COIL*0.6, 'FontWeight', 'bold');

    %% --- Halfcut plane (dashed) ---
    plot(xlim_v, [Z_CONE Z_CONE], '--', 'Color', COL_HALFCUT, 'LineWidth', LW_HALFCUT);
    text(xlim_v(2)-0.5, Z_CONE+0.4, 'halfcut plane (z = -13)', ...
         'FontSize', 10, 'Color', COL_HALFCUT, ...
         'HorizontalAlignment', 'right', 'FontAngle', 'italic');

    %% --- WP marker ---
    plot(0, SPH_OFST, 'k+', 'MarkerSize', 16, 'LineWidth', 2.4);
    text(-0.4, SPH_OFST+0.7, 'WP', 'FontSize', 13, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'right');

    %% --- Tip + region labels ---
    plot(x_tip, Z_CONE, 'k.', 'MarkerSize', 12);
    text(x_tip+0.6, Z_CONE-1.0, 'P1 tip', 'FontSize', 11, 'FontWeight', 'bold');

    text(8, -14.6, 'Cone (halfcut)', 'FontSize', 11, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center');
    text(20, -10, 'Horizontal arm', 'FontSize', 12, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center');
    text(X0_LOW, -3.5, 'Post', 'FontSize', 12, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center');
    text((POST_X2+YOKE_OUT_R)/2-1, YOKE_H/2, 'Yoke', ...
         'FontSize', 12, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(10, YOKE_H/2+0.45, '(yoke continues -> adjacent poles)', ...
         'FontSize', 10, 'Color', [0.4 0.4 0.4], 'FontAngle', 'italic', ...
         'HorizontalAlignment', 'center');

    hold off;

    %% --- Colorbar + axes ---
    colormap(turbo);
    clim([0, bxz_clip_max]);
    cb = colorbar;
    ylabel(cb, sprintf('|B_{xz}|  [Tesla]  (clipped at 97th pct)'), 'FontSize', FONT_CB);
    set(cb, 'FontSize', FONT_CB);

    axis equal;
    grid on;
    xlim(xlim_v); ylim(zlim_v);
    set(gca, 'GridAlpha', 0.18, 'Layer', 'top');
    xlabel('x [mm]  (APDL frame, +x = P1 pole axis from WP)', 'FontSize', FONT_LBL);
    ylabel('z [mm]', 'FontSize', FONT_LBL);
    title(['P1 magnetic structure (full circuit) — B-field side view, xz at y=0 — ' ...
           'sign-corrected (flux yoke->tip->WP), smrt 4 mesh'], ...
          'FontSize', FONT_TTL, 'Interpreter', 'none');

    %% --- Save ---
    out_path = fullfile(out_dir, 'P1_Bvector_full_circuit.png');
    exportgraphics(fig, out_path, 'Resolution', DPI);
    fprintf('\nSaved: %s\n', out_path);
end
