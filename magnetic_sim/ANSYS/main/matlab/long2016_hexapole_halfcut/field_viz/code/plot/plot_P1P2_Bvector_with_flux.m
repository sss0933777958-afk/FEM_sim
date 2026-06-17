function plot_P1P2_Bvector_with_flux()
% PLOT_P1P2_BVECTOR_WITH_FLUX
%   兩張 2-panel 圖,上下 x 軸對齊:
%     Top:  pole-local frame xz' B-vector quiver(turbo cmap,cone OUTLINE,**沒 sensor**)
%     Bot:  Φ(x) profile (smrt 4)
%   x_local = along pole axis (matches Φ x-axis exactly)

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');

    cnst = mt_constants();

    POLE_R     = cnst.POLE_R    * 1e3;          % 3 mm
    POLE_TIP_R = cnst.POLE_TIP_R * 1e3;         % 0.04 mm
    R_norm     = cnst.R_norm    * 1e3;          % 0.5 mm
    POLE_LEN   = cnst.POLE_CONE_LEN * 1e3;      % 15 mm
    WP_m       = [0; 0; cnst.SPH_OFST];         % m  (APDL frame)

    s_tip       = R_norm;
    s_base      = s_tip + POLE_LEN;
    s_block_end = s_base + 8;
    x_xlim      = [-2, 25];

    data_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\flux_profile';
    fig_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';

    %% ===== P1 =====
    res_dir_p1 = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\coil1\standard';
    pole_axis_p1 = [1; 0; 0];
    up_hat_p1    = [0; 0; 1];
    p1_data = load(fullfile(data_dir, 'P1_flux_profile_smrt4.mat'));
    plot_one(res_dir_p1, 'coil1', pole_axis_p1, up_hat_p1, WP_m, ...
             1, true, ...        % sign_idx=1, is_halfcut=true
             p1_data.x_query_mm, p1_data.Phi_uWb, ...
             'P1', s_tip, s_base, s_block_end, POLE_R, POLE_TIP_R, x_xlim, ...
             fullfile(fig_dir, 'P1_Bvector_with_flux.png'));

    %% ===== P2 (pole-local frame, drawn horizontal) =====
    res_dir_p2 = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\coil5\standard';
    inc      = cnst.upper_incline;
    theta2   = cnst.pole_angles(2) * pi/180;
    pole_axis_p2 = [cos(inc)*cos(theta2); cos(inc)*sin(theta2); sin(inc)];
    up_unnorm    = [0; 0; 1] - sin(inc) * pole_axis_p2;
    up_hat_p2    = up_unnorm / norm(up_unnorm);
    p2_data = load(fullfile(data_dir, 'P2_flux_profile_smrt4.mat'));
    plot_one(res_dir_p2, 'coil5', pole_axis_p2, up_hat_p2, WP_m, ...
             2, false, ...
             p2_data.s_query_mm, p2_data.Phi_uWb, ...
             'P2', s_tip, s_base, s_block_end, POLE_R, POLE_TIP_R, x_xlim, ...
             fullfile(fig_dir, 'P2_Bvector_with_flux.png'));
end


function plot_P2_tilted(cnst, x_phi, phi_arr, POLE_R, POLE_TIP_R, POLE_LEN, ...
                        s_tip, s_base, x_xlim, out_path)
% Top panel: P2 cone in ACTUAL tilted orientation (mirrored APDL xz, y=0)
% Bottom panel: Φ(x) along pole-axis arc-length(NOT linked to top axis)

    LW_STEEL = 2.0;
    inc      = cnst.upper_incline;          % 36.59 deg in rad
    beta     = atan2(3.0, 15.0);            % cone half-angle 11.31 deg
    R_norm_xy = cnst.R_norm_xy * 1e3;        % 0.408 mm
    R_norm_z  = cnst.R_norm_z  * 1e3;        % 0.289 mm
    SPH_OFST  = cnst.SPH_OFST  * 1e3;        % -12.71 mm

    res_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\coil5\standard';
    fprintf('\n=== P2 (tilted +%.2f deg) ===\n', rad2deg(inc));

    d = import_ansys_data(res_dir, 'all', 'coil5');
    bbar_mat = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\bs_matrix\Bbar_S_4p572.mat';
    sign_val = -1;
    if exist(bbar_mat, 'file')
        Bs = load(bbar_mat, 'col_sign');
        sign_val = Bs.col_sign(2);
    end
    fprintf('  Sign correction: %+d\n', sign_val);
    if sign_val < 0
        d.bx = -d.bx; d.by = -d.by; d.bz = -d.bz;
    end

    F_bx = scatteredInterpolant(d.x, d.y, d.z, d.bx, 'linear', 'none');
    F_by = scatteredInterpolant(d.x, d.y, d.z, d.by, 'linear', 'none');
    F_bz = scatteredInterpolant(d.x, d.y, d.z, d.bz, 'linear', 'none');

    % Mirrored xz grid: x' = -x_APDL so cone points right
    xlim_top = [-2, 16];
    zlim_top = [-15, -1];
    grid_x = 60; grid_z = 50;
    xc = linspace(xlim_top(1)+0.18, xlim_top(2)-0.18, grid_x);
    zc = linspace(zlim_top(1)+0.18, zlim_top(2)-0.18, grid_z);
    [Xg, Zg] = meshgrid(xc, zc); Xg = Xg(:); Zg = Zg(:);

    x_apdl = -Xg * 1e-3;
    y_apdl = zeros(size(Xg));
    z_apdl =  Zg * 1e-3;

    Bx = F_bx(x_apdl, y_apdl, z_apdl);
    By = F_by(x_apdl, y_apdl, z_apdl);
    Bz = F_bz(x_apdl, y_apdl, z_apdl);
    keep = ~isnan(Bx);
    Xg = Xg(keep); Zg = Zg(keep);
    Bx = Bx(keep); By = By(keep); Bz = Bz(keep);

    % Mirror Bx
    bx_mir = -Bx;
    B_mag = sqrt(Bx.^2 + By.^2 + Bz.^2);
    fprintf('  Grid arrows: %d\n', length(Xg));

    arrow_max = 0.45;
    scale = arrow_max / max(B_mag);
    bx_q = bx_mir * scale;
    bz_q = Bz     * scale;

    fig = figure('Position', [50 50 1400 900], 'Color', 'w');
    set(fig, 'DefaultAxesFontName', 'Helvetica', 'DefaultAxesFontSize', 11);

    %% --- Panel A: B-vector (actual tilted) ---
    ax1 = subplot(2, 1, 1);
    hold on; box on;

    n_bins = 24;
    edges_b = linspace(0, max(B_mag), n_bins+1);
    cmap_n = turbo(n_bins);
    lw_rng = [0.6, 2.2];
    for k = 1:n_bins
        in_bin = B_mag >= edges_b(k) & B_mag < edges_b(k+1);
        if k == n_bins, in_bin = in_bin | (B_mag >= edges_b(end)); end
        if any(in_bin)
            lw = lw_rng(1) + (k-1)/(n_bins-1)*(lw_rng(2)-lw_rng(1));
            quiver(Xg(in_bin), Zg(in_bin), bx_q(in_bin), bz_q(in_bin), 0, ...
                   'Color', cmap_n(k,:), 'LineWidth', lw, 'MaxHeadSize', 0.6);
        end
    end

    % Cone outline in mirrored xz (tilted +inc with cone half-angle beta)
    T = [R_norm_xy, SPH_OFST + R_norm_z];
    L_slant = sqrt(POLE_LEN^2 + POLE_R^2);
    C_top = T + L_slant * [cos(inc + beta), sin(inc + beta)];
    C_bot = T + L_slant * [cos(inc - beta), sin(inc - beta)];
    plot([T(1), C_top(1)], [T(2), C_top(2)], 'k-', 'LineWidth', LW_STEEL);
    plot([T(1), C_bot(1)], [T(2), C_bot(2)], 'k-', 'LineWidth', LW_STEEL);
    plot([C_top(1), C_bot(1)], [C_top(2), C_bot(2)], 'k-', 'LineWidth', LW_STEEL);

    % WP marker
    plot(0, SPH_OFST, 'k+', 'MarkerSize', 14, 'LineWidth', 2.2);
    text(0.3, SPH_OFST - 0.35, 'WP', 'FontSize', 12, 'FontWeight', 'bold');

    % Tip label
    text(T(1) + 0.2, T(2) + 0.4, 'P2 tip', 'FontSize', 10, ...
         'Color', 'r', 'FontWeight', 'bold');
    % Cone end label (top slant endpoint)
    text(C_top(1) - 0.1, C_top(2) + 0.4, 'cone end', 'FontSize', 10, ...
         'Color', [0.10 0.40 0.85], 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'right');

    colormap(turbo);
    clim([0, max(B_mag)]);
    cb = colorbar;
    ylabel(cb, '|B| [T]', 'FontSize', 11);

    axis equal;
    xlim(xlim_top); ylim(zlim_top);
    xlabel('x'' [mm]  (= -x_{APDL}, mirrored so P2 cone points right)', ...
           'FontSize', 12, 'Interpreter', 'tex');
    ylabel('z [mm]  (APDL frame)', 'FontSize', 12);
    title(sprintf('A. P2 B-vector — actual orientation, upper pole tilted +%.2f°', ...
                  rad2deg(inc)), 'FontSize', 13, 'Interpreter', 'none');
    grid on;

    %% --- Panel B: Φ(x) curve(along pole axis arc-length)---
    ax2 = subplot(2, 1, 2);
    hold on; box on;
    plot(x_phi, phi_arr, 'b-', 'LineWidth', 2.2);
    ymax = max(phi_arr)*1.10; ymin = -ymax*0.05;
    ylim([ymin, ymax]);
    plot([0 0], [ymin ymax], 'k--', 'LineWidth', 1.0);
    plot([s_tip s_tip], [ymin ymax], 'r--', 'LineWidth', 1.0);
    plot([s_base s_base], [ymin ymax], '--', 'Color', [0.10 0.40 0.85], ...
         'LineWidth', 1.0);
    text(0, ymax*0.95, 'WP', 'FontSize', 10);
    text(s_tip, ymax*0.85, 'P2 tip', 'FontSize', 10, 'Color', 'r');
    text(s_base, ymax*0.85, 'cone end', 'FontSize', 10, ...
         'Color', [0.10 0.40 0.85]);
    [phi_max, k_max] = max(phi_arr);
    plot(x_phi(k_max), phi_max, 'bo', 'MarkerSize', 9, 'MarkerFaceColor', 'b');
    text(x_phi(k_max)+0.4, phi_max+0.05, ...
         sprintf('\\Phi_{max} = %.2f \\muWb @ x = %.1f mm', phi_max, x_phi(k_max)), ...
         'FontSize', 11, 'FontWeight', 'bold', 'Color', 'b', 'Interpreter', 'tex');
    xlim(x_xlim);
    xlabel('x [mm]  (pole axis arc-length from WP — bottom panel only)', ...
           'FontSize', 12);
    ylabel('\Phi(x) [\muWb]', 'FontSize', 12, 'Interpreter', 'tex');
    title('B. P2 axial flux Φ(x) = ∫ B_{axial} dA  (along pole axis arc-length)', ...
          'FontSize', 13);
    grid on;

    exportgraphics(fig, out_path, 'Resolution', 200);
    fprintf('Saved: %s\n', out_path);
end


function plot_one(res_dir, coilname, pole_axis, up_hat, WP_m, ...
                  sign_idx, is_halfcut, ...
                  x_phi, phi_arr, ...
                  label, s_tip, s_base, s_block_end, POLE_R, POLE_TIP_R, x_xlim, ...
                  out_path)
    LW_STEEL = 2.0;
    fprintf('\n=== %s ===\n', label);

    %% Load FEM + sign correction
    d = import_ansys_data(res_dir, 'all', coilname);
    bbar_mat = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\bs_matrix\Bbar_S_4p572.mat';
    sign_val = +1;
    if exist(bbar_mat, 'file')
        Bs = load(bbar_mat, 'col_sign');
        sign_val = Bs.col_sign(sign_idx);
    end
    fprintf('  Sign correction: %+d\n', sign_val);
    if sign_val < 0
        d.bx = -d.bx; d.by = -d.by; d.bz = -d.bz;
    end

    F_bx = scatteredInterpolant(d.x, d.y, d.z, d.bx, 'linear', 'none');
    F_by = scatteredInterpolant(d.x, d.y, d.z, d.by, 'linear', 'none');
    F_bz = scatteredInterpolant(d.x, d.y, d.z, d.bz, 'linear', 'none');

    %% Pole-local grid
    x_loc_v = linspace(x_xlim(1), x_xlim(2), 80);
    y_loc_v = linspace(-4, 4, 32);
    [Xg, Yg] = meshgrid(x_loc_v, y_loc_v);
    Xg = Xg(:); Yg = Yg(:);

    % pole-local → APDL: pos_apdl = WP + x_loc·pole_axis + y_loc·up_hat
    Xa_m = WP_m(1) + (Xg*1e-3)*pole_axis(1) + (Yg*1e-3)*up_hat(1);
    Ya_m = WP_m(2) + (Xg*1e-3)*pole_axis(2) + (Yg*1e-3)*up_hat(2);
    Za_m = WP_m(3) + (Xg*1e-3)*pole_axis(3) + (Yg*1e-3)*up_hat(3);

    Bx = F_bx(Xa_m, Ya_m, Za_m);
    By = F_by(Xa_m, Ya_m, Za_m);
    Bz = F_bz(Xa_m, Ya_m, Za_m);

    keep = ~isnan(Bx);
    Xg = Xg(keep); Yg = Yg(keep);
    Bx = Bx(keep); By = By(keep); Bz = Bz(keep);

    % Project to pole-local axes
    B_along = Bx*pole_axis(1) + By*pole_axis(2) + Bz*pole_axis(3);
    B_perp  = Bx*up_hat(1)   + By*up_hat(2)   + Bz*up_hat(3);
    B_mag   = sqrt(Bx.^2 + By.^2 + Bz.^2);

    fprintf('  Grid arrows (after NaN filter): %d\n', length(Xg));

    % Clip at 92nd percentile so moderate-field arrows stay visible
    % (Long2016 verbatim concentrates |B| near cone tip — without clip,
    %  arrow scale is hijacked by the tip and rest becomes invisible)
    b_clip = quantile(B_mag, 0.92);
    fprintf('  |B| range [%.3g, %.3g] T, clip @ 92nd pct = %.3g T\n', ...
            min(B_mag), max(B_mag), b_clip);
    arrow_max = 0.45;
    B_mag_for_scale = min(B_mag, b_clip);
    scale = arrow_max / b_clip;
    bx_q = B_along .* (B_mag_for_scale ./ max(B_mag, eps)) * scale;
    by_q = B_perp  .* (B_mag_for_scale ./ max(B_mag, eps)) * scale;

    %% Plot
    fig = figure('Position', [50 50 1400 760], 'Color', 'w');
    set(fig, 'DefaultAxesFontName', 'Helvetica', 'DefaultAxesFontSize', 11);

    %% --- Panel A: B-vector quiver (top) ---
    ax1 = subplot(2, 1, 1);
    hold on; box on;

    n_bins  = 24;
    edges_b = linspace(0, b_clip, n_bins + 1);
    cmap_n  = turbo(n_bins);
    lw_rng  = [0.6, 2.2];
    B_mag_binned = min(B_mag, b_clip);
    for k = 1:n_bins
        in_bin = B_mag_binned >= edges_b(k) & B_mag_binned < edges_b(k+1);
        if k == n_bins
            in_bin = in_bin | (B_mag_binned >= edges_b(end));
        end
        if any(in_bin)
            lw = lw_rng(1) + (k-1)/(n_bins-1) * (lw_rng(2)-lw_rng(1));
            quiver(Xg(in_bin), Yg(in_bin), bx_q(in_bin), by_q(in_bin), 0, ...
                   'Color', cmap_n(k,:), 'LineWidth', lw, 'MaxHeadSize', 0.6);
        end
    end

    % Cone OUTLINE in pole-local frame
    if is_halfcut
        % D-shape (lower half only)
        % Tip (s_tip, 0) → lower (s_base, -POLE_R) → block (s_block_end, -POLE_R) → (s_block_end, 0) → tip
        plot([s_tip, s_base, s_block_end, s_block_end, s_tip], ...
             [0,    -POLE_R, -POLE_R,    0,           0], ...
             'k-', 'LineWidth', LW_STEEL);
        % milled flat = top edge (already covered by last segment)
    else
        % Full cone: top + bottom slants + block
        plot([s_tip, s_base, s_block_end, s_block_end, s_base, s_tip], ...
             [0,     POLE_R, POLE_R,     -POLE_R,     -POLE_R, 0], ...
             'k-', 'LineWidth', LW_STEEL);
    end

    % WP marker — only for P1 (per user: P2 has WP removed)
    if strcmp(label, 'P1')
        plot(0, 0, 'k+', 'MarkerSize', 14, 'LineWidth', 2.2);
        text(0.2, 0.5, 'WP', 'FontSize', 12, 'FontWeight', 'bold');
    end

    % Cone end / tip annotations (no sensor)
    text(s_tip, -0.6, sprintf('%s tip\n(x=%.2f)', label, s_tip), ...
         'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center');
    text(s_base, POLE_R + 0.4, sprintf('cone end\n(x=%.1f)', s_base), ...
         'FontSize', 10, 'Color', [0.10 0.40 0.85], 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center');

    colormap(turbo);
    clim([0, b_clip]);
    cb = colorbar;
    ylabel(cb, sprintf('|B| [T]  (clip @ 92nd pct = %.2f T)', b_clip), 'FontSize', 11);

    xlim(x_xlim);
    ylim([-POLE_R - 1, POLE_R + 1.2]);
    ylabel('y_{local} [mm]', 'FontSize', 12, 'Interpreter', 'tex');
    if is_halfcut
        ttl = sprintf('A. %s B-vector (pole-local frame, halfcut D-shape, Long2016 verbatim)', label);
    else
        ttl = sprintf('A. %s B-vector (pole-local frame, full cone, Long2016 verbatim)', label);
    end
    title(ttl, 'FontSize', 13);
    grid on;
    set(ax1, 'XTickLabel', []);

    %% --- Panel B: Φ(x) curve (bottom) ---
    ax2 = subplot(2, 1, 2);
    hold on; box on;
    plot(x_phi, phi_arr, 'b-', 'LineWidth', 2.2);
    ymax = max(phi_arr)*1.10;  ymin = -ymax*0.05;
    ylim([ymin, ymax]);
    if strcmp(label, 'P1')
        plot([0 0], [ymin ymax], 'k--', 'LineWidth', 1.0);
        text(0, ymax*0.95, 'WP', 'FontSize', 10);
    end
    plot([s_tip s_tip], [ymin ymax], 'r--', 'LineWidth', 1.0);
    plot([s_base s_base], [ymin ymax], '--', 'Color', [0.10 0.40 0.85], ...
         'LineWidth', 1.0);
    text(s_tip, ymax*0.85, sprintf('%s tip', label), 'FontSize', 10, 'Color', 'r');
    text(s_base, ymax*0.85, 'cone end', 'FontSize', 10, 'Color', [0.10 0.40 0.85]);

    [phi_max, k_max] = max(phi_arr);
    plot(x_phi(k_max), phi_max, 'bo', 'MarkerSize', 9, 'MarkerFaceColor', 'b');
    text(x_phi(k_max)+0.4, phi_max+0.05, ...
         sprintf('\\Phi_{max} = %.2f \\muWb @ x = %.1f mm', phi_max, x_phi(k_max)), ...
         'FontSize', 11, 'FontWeight', 'bold', 'Color', 'b', 'Interpreter', 'tex');

    xlim(x_xlim);
    xlabel('x [mm]   (pole axis arc-length from WP)', 'FontSize', 12);
    ylabel('\Phi(x) [\muWb]', 'FontSize', 12, 'Interpreter', 'tex');
    title(sprintf('B. %s axial flux Φ(x) = ∫ B_{axial} dA', label), 'FontSize', 13);
    grid on;

    linkaxes([ax1, ax2], 'x');

    exportgraphics(fig, out_path, 'Resolution', 200);
    fprintf('Saved: %s\n', out_path);
end
