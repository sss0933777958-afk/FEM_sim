function plot_P2_circuit_with_flux()
% PLOT_P2_CIRCUIT_WITH_FLUX
% Two-panel: top = P2 side view (APDL xz, y=0, upper pole tilted +36.59°)
%             with B-field vectors; bottom = axial flux profile.
%
% Data: coil5 (paper P2) Long2016 verbatim baseline.

    %% --- Paths ---
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    cnst = mt_constants();

    res_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data\coil5\standard';
    out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    %% --- Geometry (mm, APDL) ---
    R_norm_xy_mm  = cnst.R_norm_xy * 1e3;
    POLE_R        = cnst.POLE_R    * 1e3;
    POLE_LEN      = cnst.POLE_CONE_LEN * 1e3;
    POLE_TIP_R_mm = cnst.POLE_TIP_R * 1e3;
    SPH_OFST      = cnst.SPH_OFST  * 1e3;       % -12.711
    YOKE_H        = 2;
    PROT_H        = 7;
    inc           = cnst.upper_incline;          % 36.59° rad
    theta_p2      = pi;                          % azimuth 180°

    pole_axis = [cos(inc)*cos(theta_p2); cos(inc)*sin(theta_p2); sin(inc)];   % (-0.803, 0, 0.596)
    up_unn    = [0;0;1] - sin(inc)*pole_axis;
    up_hat    = up_unn/norm(up_unn);
    pa_xz     = [pole_axis(1); pole_axis(3)];
    uh_xz     = [up_hat(1);   up_hat(3)];
    WP_xz     = [0; SPH_OFST];

    % Actual upper tip in APDL:
    tip_z_upper  = -PROT_H - 6 + cnst.R_norm_z*1e3*2;   % -12.422
    tip_xz       = [R_norm_xy_mm*cos(theta_p2); tip_z_upper];
    cone_end_xz  = tip_xz + POLE_LEN * pa_xz;                  % cone-cylinder boundary
    cone_top     = cone_end_xz + POLE_R * uh_xz;
    cone_bot     = cone_end_xz - POLE_R * uh_xz;

    % Cylinder from cone_end to end_upper (pole geometry end inside block)
    end_upper_xz = [-47.5 + 11.5; YOKE_H + PROT_H + 5];        % = (-36, 14)
    end_upper_xz(1) = -36;                                      % verify
    cyl_top      = end_upper_xz + POLE_R * uh_xz;
    cyl_bot      = end_upper_xz - POLE_R * uh_xz;

    % Upper post / block region (APDL):
    Z0_HIGH      = YOKE_H + PROT_H;            % 9
    X0_HIGH      = 47.5 * cos(theta_p2);       % -47.5
    POST_Z1      = YOKE_H;                     % 2 (yoke top)
    POST_Z2      = Z0_HIGH;                    % 9 (post top = block bottom)
    BLOCK_Z1     = Z0_HIGH;                    % 9
    BLOCK_Z2     = Z0_HIGH + 10;               % 19
    BLOCK_X1     = X0_HIGH - 16.5;             % -64
    BLOCK_X2     = X0_HIGH + 8.5;              % -39
    POST_X1      = X0_HIGH - PROT_H/2;         % -50.5
    POST_X2      = X0_HIGH + PROT_H/2;         % -44.5

    %% --- Load FEM ---
    fprintf('Loading coil5 baseline...\n');
    d = import_ansys_data(res_dir, 'all', 'coil5');

    bbar_mat = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\bs_matrix\Bbar_S_4p572.mat';
    sign_p2 = +1;
    if exist(bbar_mat, 'file')
        Bs = load(bbar_mat, 'col_sign');
        sign_p2 = Bs.col_sign(2);
    end
    fprintf('P2 sign correction: %+d\n', sign_p2);
    if sign_p2 < 0
        d.bx = -d.bx; d.by = -d.by; d.bz = -d.bz;
    end

    F_bx = scatteredInterpolant(d.x, d.y, d.z, d.bx, 'linear', 'none');
    F_by = scatteredInterpolant(d.x, d.y, d.z, d.by, 'linear', 'none');
    F_bz = scatteredInterpolant(d.x, d.y, d.z, d.bz, 'linear', 'none');

    %% --- Figure ---
    fig = figure('Position', [50 30 1500 900], 'Color', 'w');

    %% --- Top: side view + B vectors ---
    ax1 = subplot('Position', [0.07 0.50 0.86 0.45]);
    hold on;

    % P2 cone (tilted triangle, tip → base_top → base_bot → tip)
    fill([tip_xz(1), cone_top(1), cone_bot(1), tip_xz(1)], ...
         [tip_xz(2), cone_top(2), cone_bot(2), tip_xz(2)], ...
         [0.85 0.85 0.85], 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0.7);

    % P2 cylinder (tilted rectangle from cone_end to end_upper)
    fill([cone_top(1), cyl_top(1), cyl_bot(1), cone_bot(1)], ...
         [cone_top(2), cyl_top(2), cyl_bot(2), cone_bot(2)], ...
         [0.85 0.85 0.85], 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0.7);

    % Block (z=9-19, x=BLOCK_X1..BLOCK_X2)
    fill([BLOCK_X1 BLOCK_X2 BLOCK_X2 BLOCK_X1], ...
         [BLOCK_Z1 BLOCK_Z1 BLOCK_Z2 BLOCK_Z2], [0.85 0.85 0.85], ...
         'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.7);

    % Post (z=2-9, x range = ±PROT_R around X0_HIGH)
    fill([POST_X1 POST_X2 POST_X2 POST_X1], ...
         [POST_Z1 POST_Z1 POST_Z2 POST_Z2], [0.85 0.85 0.85], ...
         'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.7);

    % Yoke (rectangle at z=0-2, full width)
    fill([BLOCK_X1-2 BLOCK_X2+2 BLOCK_X2+2 BLOCK_X1-2], ...
         [0 0 YOKE_H YOKE_H], [0.85 0.85 0.85], ...
         'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.7);

    % Coil (orange, z=POST_Z1 to POST_Z2, flanking post)
    COIL_IN_R = 5; COIL_OUT_R = 8;
    COIL_Z1 = POST_Z1; COIL_Z2 = POST_Z2;
    fill([X0_HIGH-COIL_OUT_R X0_HIGH-COIL_IN_R X0_HIGH-COIL_IN_R X0_HIGH-COIL_OUT_R], ...
         [COIL_Z1 COIL_Z1 COIL_Z2 COIL_Z2], [0.93 0.55 0.10], ...
         'EdgeColor','k','LineWidth',1.5,'FaceAlpha',0.55);
    fill([X0_HIGH+COIL_IN_R X0_HIGH+COIL_OUT_R X0_HIGH+COIL_OUT_R X0_HIGH+COIL_IN_R], ...
         [COIL_Z1 COIL_Z1 COIL_Z2 COIL_Z2], [0.93 0.55 0.10], ...
         'EdgeColor','k','LineWidth',1.5,'FaceAlpha',0.55);
    COIL_Z_CTR = (COIL_Z1+COIL_Z2)/2;
    text(X0_HIGH-(COIL_IN_R+COIL_OUT_R)/2, COIL_Z_CTR, char(8857), 'FontSize',16, ...
         'HorizontalAlignment','center','FontWeight','bold');
    text(X0_HIGH+(COIL_IN_R+COIL_OUT_R)/2, COIL_Z_CTR, char(8855), 'FontSize',16, ...
         'HorizontalAlignment','center','FontWeight','bold');

    % B field quiver
    xlim_v = [-68, 5];  zlim_v = [-16, 22];
    grid_x = 120; grid_z = 70;
    xc = linspace(xlim_v(1), xlim_v(2), grid_x);
    zc = linspace(zlim_v(1), zlim_v(2), grid_z);
    [Xg, Zg] = meshgrid(xc, zc); Xg = Xg(:); Zg = Zg(:);
    Bx_q = F_bx(Xg*1e-3, zeros(size(Xg)), Zg*1e-3);
    Bz_q = F_bz(Xg*1e-3, zeros(size(Xg)), Zg*1e-3);
    keep = ~isnan(Bx_q) & ~isnan(Bz_q);
    Xg = Xg(keep); Zg = Zg(keep); Bx_q = Bx_q(keep); Bz_q = Bz_q(keep);
    Bmag = sqrt(Bx_q.^2 + Bz_q.^2);

    b_clip = quantile(Bmag, 0.92);
    dx_grid = (xlim_v(2)-xlim_v(1)) / grid_x;
    arrow_max = 0.85 * dx_grid;
    b_norm = max(Bmag, eps);
    b_eff  = min(Bmag, b_clip);
    scale_per = (arrow_max / b_clip) .* (b_eff ./ b_norm);
    bxq = Bx_q .* scale_per;  bzq = Bz_q .* scale_per;

    n_bins = 24;
    edges_b = linspace(0, b_clip, n_bins+1);
    cmap_n = turbo(n_bins);
    lw_rng = [0.55, 1.90];
    B_binned = min(Bmag, b_clip);
    for k = 1:n_bins
        in_bin = B_binned >= edges_b(k) & B_binned < edges_b(k+1);
        if k == n_bins; in_bin = in_bin | (B_binned >= edges_b(end)); end
        if any(in_bin)
            lw = lw_rng(1) + (k-1)/(n_bins-1)*(lw_rng(2)-lw_rng(1));
            quiver(Xg(in_bin), Zg(in_bin), bxq(in_bin), bzq(in_bin), 0, ...
                   'Color', cmap_n(k,:), 'LineWidth', lw, 'MaxHeadSize', 0.6);
        end
    end

    plot(WP_xz(1), WP_xz(2), 'k+', 'MarkerSize', 14, 'LineWidth', 2);
    text(WP_xz(1)+1, WP_xz(2)-0.3, 'WP', 'FontSize', 12, 'FontWeight','bold');

    colormap(turbo); clim([0 b_clip]);
    cb = colorbar('Position', [0.945 0.50 0.012 0.45]);
    ylabel(cb, sprintf('|B| [T] (clip 92pct = %.2f T)', b_clip), 'FontSize', 10);

    axis equal; grid on; box on;
    xlim(xlim_v); ylim(zlim_v);
    xlabel('x [mm]', 'FontSize', 11);
    ylabel('z [mm]', 'FontSize', 11);
    title('P2 magnetic circuit side view (Long2016 verbatim, upper pole tilted +36.59°)', 'FontSize', 13);

    %% --- Bottom: axial flux profile (full disc, R(s) taper) ---
    ax2 = subplot('Position', [0.07 0.07 0.86 0.32]);
    hold on;

    R_norm_mm = cnst.R_norm * 1e3;   % 0.5
    cone_len_mm = POLE_LEN;
    s_tip = R_norm_mm; s_base = R_norm_mm + cone_len_mm;
    R_func = @(s) clamp_R(s, s_tip, s_base, POLE_TIP_R_mm, POLE_R);

    s_query_mm = linspace(-2, 37.5, 280);
    Phi = zeros(size(s_query_mm));
    N_grid = 81;
    WP_apdl = [0;0;SPH_OFST*1e-3];
    side_hat = cross(pole_axis, up_hat); side_hat = side_hat/norm(side_hat);

    for k = 1:length(s_query_mm)
        s = s_query_mm(k);
        R_mm = R_func(s);
        center_apdl = WP_apdl + (s*1e-3) * pole_axis;
        u_c = linspace(-R_mm, R_mm, N_grid); v_c = linspace(-R_mm, R_mm, N_grid);
        [Ug, Vg] = meshgrid(u_c, v_c);
        R_pts = sqrt(Ug.^2 + Vg.^2);
        in_disc = R_pts <= R_mm;
        du = (u_c(2)-u_c(1))*1e-3; dv = (v_c(2)-v_c(1))*1e-3; dA = du*dv;

        Xq = center_apdl(1) + (Ug*1e-3)*up_hat(1) + (Vg*1e-3)*side_hat(1);
        Yq = center_apdl(2) + (Ug*1e-3)*up_hat(2) + (Vg*1e-3)*side_hat(2);
        Zq = center_apdl(3) + (Ug*1e-3)*up_hat(3) + (Vg*1e-3)*side_hat(3);

        Bx_q = F_bx(Xq(:), Yq(:), Zq(:));
        By_q = F_by(Xq(:), Yq(:), Zq(:));
        Bz_q = F_bz(Xq(:), Yq(:), Zq(:));
        B_axial = Bx_q*pole_axis(1) + By_q*pole_axis(2) + Bz_q*pole_axis(3);
        B_axial = reshape(B_axial, N_grid, N_grid);
        B_axial(~in_disc) = 0; B_axial(isnan(B_axial)) = 0;
        Phi(k) = sum(B_axial(in_disc)) * dA;
    end
    Phi_uWb = abs(Phi) * 1e6;

    plot(s_query_mm, Phi_uWb, 'b-', 'LineWidth', 2.2);

    ymax = max(Phi_uWb)*1.10; ymin = -ymax*0.05;
    plot([0 0], [ymin ymax], 'k--', 'LineWidth', 1.4);
    plot([s_tip s_tip], [ymin ymax], '--', 'Color', [0.92 0.15 0.15], 'LineWidth', 1.4);
    plot([s_base s_base], [ymin ymax], '--', 'Color', [0.15 0.55 0.85], 'LineWidth', 1.4);
    text(0+0.3, ymax*0.93, 'WP', 'FontSize', 11, 'FontWeight','bold');
    text(s_tip+0.3, ymax*0.85, 'P2 tip', 'FontSize', 10, 'Color', [0.92 0.15 0.15], 'FontWeight','bold');
    text(s_base+0.3, ymax*0.85, 'cone end', 'FontSize', 10, 'Color', [0.15 0.55 0.85], 'FontWeight','bold');

    [phi_pk, k_pk] = max(Phi_uWb);
    plot(s_query_mm(k_pk), phi_pk, 'bo', 'MarkerSize', 8, 'MarkerFaceColor','b');
    text(s_query_mm(k_pk)-0.6, phi_pk+0.06, ...
         sprintf('\\Phi_{max} = %.2f \\muWb @ s = %.1f mm', phi_pk, s_query_mm(k_pk)), ...
         'FontSize', 11, 'FontWeight','bold', 'Color','b', 'Interpreter','tex', ...
         'HorizontalAlignment','right');

    grid on; box on;
    xlim([-2 37.5]); ylim([ymin ymax]);
    xlabel('s [mm]  (P2 pole axis arc-length, WP \rightarrow tip \rightarrow cone end \rightarrow cylinder end)', ...
           'FontSize', 11, 'Interpreter','tex');
    ylabel('|\Phi(s)| = |\int B_{axial} dA| [\muWb]', 'FontSize', 11, 'Interpreter','tex');
    title('P2 axial flux profile (full disc with R(s) taper)', 'FontSize', 12, 'Interpreter','tex');

    %% --- Save ---
    out_path = fullfile(out_dir, 'P2_circuit_with_flux.png');
    exportgraphics(fig, out_path, 'Resolution', 250);
    fprintf('Saved: %s\n', out_path);
end


function R = clamp_R(x_mm, x_tip, x_base, R_tip, R_max)
    if x_mm <= x_tip
        R = R_tip;
    elseif x_mm >= x_base
        R = R_max;
    else
        R = R_tip + (x_mm-x_tip)/(x_base-x_tip) * (R_max-R_tip);
    end
end
