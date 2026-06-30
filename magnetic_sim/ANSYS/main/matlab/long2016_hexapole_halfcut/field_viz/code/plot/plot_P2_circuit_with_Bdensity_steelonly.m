function plot_P2_circuit_with_Bdensity_steelonly()
% PLOT_P2_CIRCUIT_WITH_BDENSITY_STEELONLY
% Mirror of plot_P2_circuit_with_flux_steelonly.m, but the bottom panel
% shows the mean axial flux DENSITY <B_axial>(s) = Phi(s) / A_disc(s)
% instead of the integrated flux Phi(s).
%
%   <B_axial>(s) = ( sum_{(u,v) in D(s)} B_axial(s, u, v) * dA ) / A_disc(s)
%   A_disc(s)    = pi * R(s)^2     (FULL disc, P2 is not halfcut)
%
% Steel-only FEM filter for the Phi numerator; top side view uses full
% interpolant so air-side vectors still render.

    %% --- Paths ---
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    cnst = mt_constants();

    res_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data\coil5\standard';
    out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    %% --- Geometry (mm) ---
    R_norm_xy_mm  = cnst.R_norm_xy * 1e3;
    POLE_R        = cnst.POLE_R    * 1e3;
    POLE_LEN      = cnst.POLE_CONE_LEN * 1e3;
    POLE_TIP_R_mm = cnst.POLE_TIP_R * 1e3;
    SPH_OFST      = cnst.SPH_OFST  * 1e3;       % -12.711
    YOKE_H        = 2;
    PROT_H        = 7;
    inc           = cnst.upper_incline;
    theta_p2      = pi;

    pole_axis = [cos(inc)*cos(theta_p2); cos(inc)*sin(theta_p2); sin(inc)];
    up_unn    = [0;0;1] - sin(inc)*pole_axis;
    up_hat    = up_unn/norm(up_unn);
    pa_xz     = [pole_axis(1); pole_axis(3)];
    uh_xz     = [up_hat(1);   up_hat(3)];
    WP_xz     = [0; SPH_OFST];

    tip_z_upper  = -PROT_H - 6 + cnst.R_norm_z*1e3*2;
    tip_xz       = [R_norm_xy_mm*cos(theta_p2); tip_z_upper];
    cone_end_xz  = tip_xz + POLE_LEN * pa_xz;
    cone_top     = cone_end_xz + POLE_R * uh_xz;
    cone_bot     = cone_end_xz - POLE_R * uh_xz;

    end_upper_xz = [-36; YOKE_H + PROT_H + 5];
    cyl_top      = end_upper_xz + POLE_R * uh_xz;
    cyl_bot      = end_upper_xz - POLE_R * uh_xz;

    Z0_HIGH      = YOKE_H + PROT_H;
    X0_HIGH      = 47.5 * cos(theta_p2);
    POST_Z1      = YOKE_H;
    POST_Z2      = Z0_HIGH;
    BLOCK_Z1     = Z0_HIGH;
    BLOCK_Z2     = Z0_HIGH + 10;
    BLOCK_X1     = X0_HIGH - 16.5;
    BLOCK_X2     = X0_HIGH + 8.5;
    POST_X1      = X0_HIGH - PROT_H/2;
    POST_X2      = X0_HIGH + PROT_H/2;
    POST_R_mm    = PROT_H/2;

    %% Cylinder length along pole axis
    cyl_vec_xz = end_upper_xz - cone_end_xz;
    L_cyl_mm   = norm(cyl_vec_xz);
    L_pole_total_mm = POLE_LEN + L_cyl_mm;

    %% --- Load FEM ---
    fprintf('Loading coil5 baseline...\n');
    d = import_ansys_data(res_dir, 'all', 'coil5');
    n_all = length(d.x);
    fprintf('  Total FEM nodes: %d\n', n_all);

    %% --- P2 sign correction ---
    bbar_mat = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\bs_matrix\data\Bbar_S_4p572.mat';
    sign_p2 = +1;
    if exist(bbar_mat, 'file')
        Bs = load(bbar_mat, 'col_sign');
        sign_p2 = Bs.col_sign(2);
    end
    fprintf('P2 sign correction: %+d\n', sign_p2);
    if sign_p2 < 0
        d.bx = -d.bx; d.by = -d.by; d.bz = -d.bz;
    end

    %% --- Filter to STEEL only ---
    tip_apdl = [R_norm_xy_mm*cos(theta_p2)*1e-3; 0; tip_z_upper*1e-3];
    is_steel = mask_p2_steel(d.x, d.y, d.z, tip_apdl, pole_axis, ...
                              POLE_TIP_R_mm*1e-3, POLE_R*1e-3, ...
                              POLE_LEN*1e-3, L_pole_total_mm*1e-3, ...
                              [BLOCK_X1-2, BLOCK_X2+2]*1e-3, [0, YOKE_H]*1e-3, ...
                              [BLOCK_X1, BLOCK_X2]*1e-3, [BLOCK_Z1, BLOCK_Z2]*1e-3, ...
                              [X0_HIGH; 0]*1e-3, POST_R_mm*1e-3, [POST_Z1, POST_Z2]*1e-3);
    n_steel = sum(is_steel);
    fprintf('  Steel nodes (filtered): %d (=%.1f%%)\n', n_steel, 100*n_steel/n_all);

    d_st.x  = d.x(is_steel);   d_st.y  = d.y(is_steel);   d_st.z  = d.z(is_steel);
    d_st.bx = d.bx(is_steel);  d_st.by = d.by(is_steel);  d_st.bz = d.bz(is_steel);

    %% --- Build TWO interpolants ---
    fprintf('Building scatteredInterpolant (full, for side view)...\n');
    F_bx_full = scatteredInterpolant(d.x, d.y, d.z, d.bx, 'linear', 'none');
    F_bz_full = scatteredInterpolant(d.x, d.y, d.z, d.bz, 'linear', 'none');
    fprintf('Building scatteredInterpolant (steel-only, for <B_axial> integration)...\n');
    F_bx = scatteredInterpolant(d_st.x, d_st.y, d_st.z, d_st.bx, 'linear', 'none');
    F_by = scatteredInterpolant(d_st.x, d_st.y, d_st.z, d_st.by, 'linear', 'none');
    F_bz = scatteredInterpolant(d_st.x, d_st.y, d_st.z, d_st.bz, 'linear', 'none');

    %% --- Figure ---
    fig = figure('Position', [50 30 1500 900], 'Color', 'w');

    %% --- Top: side view + B vectors (identical to flux script) ---
    ax1 = subplot('Position', [0.07 0.50 0.86 0.45]);
    hold on;

    fill([tip_xz(1), cone_top(1), cone_bot(1), tip_xz(1)], ...
         [tip_xz(2), cone_top(2), cone_bot(2), tip_xz(2)], ...
         [0.85 0.85 0.85], 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0.7);
    fill([cone_top(1), cyl_top(1), cyl_bot(1), cone_bot(1)], ...
         [cone_top(2), cyl_top(2), cyl_bot(2), cone_bot(2)], ...
         [0.85 0.85 0.85], 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0.7);
    fill([BLOCK_X1 BLOCK_X2 BLOCK_X2 BLOCK_X1], ...
         [BLOCK_Z1 BLOCK_Z1 BLOCK_Z2 BLOCK_Z2], [0.85 0.85 0.85], ...
         'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.7);
    fill([POST_X1 POST_X2 POST_X2 POST_X1], ...
         [POST_Z1 POST_Z1 POST_Z2 POST_Z2], [0.85 0.85 0.85], ...
         'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.7);
    fill([BLOCK_X1-2 BLOCK_X2+2 BLOCK_X2+2 BLOCK_X1-2], ...
         [0 0 YOKE_H YOKE_H], [0.85 0.85 0.85], ...
         'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.7);

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

    xlim_v = [-68, 5];  zlim_v = [-16, 22];
    grid_x = 120; grid_z = 70;
    xc = linspace(xlim_v(1), xlim_v(2), grid_x);
    zc = linspace(zlim_v(1), zlim_v(2), grid_z);
    [Xg, Zg] = meshgrid(xc, zc); Xg = Xg(:); Zg = Zg(:);
    Bx_q = F_bx_full(Xg*1e-3, zeros(size(Xg)), Zg*1e-3);
    Bz_q = F_bz_full(Xg*1e-3, zeros(size(Xg)), Zg*1e-3);
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
    title('P2 magnetic circuit side view  (steel-only interp, upper pole tilted +36.59°)', 'FontSize', 13);

    %% --- Bottom: mean axial flux DENSITY <B_axial>(s) ---
    ax2 = subplot('Position', [0.07 0.07 0.86 0.32]);
    hold on;

    R_norm_mm = cnst.R_norm * 1e3;
    s_tip = R_norm_mm; s_base = R_norm_mm + POLE_LEN;
    R_func = @(s) clamp_R(s, s_tip, s_base, POLE_TIP_R_mm, POLE_R);

    s_query_mm = linspace(-2, 37.5, 280);
    Phi   = zeros(size(s_query_mm));
    Adisc = zeros(size(s_query_mm));   % FULL disc area [m^2]
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

        Bx_q_i = F_bx(Xq(:), Yq(:), Zq(:));
        By_q_i = F_by(Xq(:), Yq(:), Zq(:));
        Bz_q_i = F_bz(Xq(:), Yq(:), Zq(:));
        B_axial = Bx_q_i*pole_axis(1) + By_q_i*pole_axis(2) + Bz_q_i*pole_axis(3);
        B_axial = reshape(B_axial, N_grid, N_grid);
        B_axial(~in_disc) = 0; B_axial(isnan(B_axial)) = 0;
        Phi(k)   = sum(B_axial(in_disc)) * dA;
        Adisc(k) = pi * (R_mm*1e-3)^2;   % FULL disc analytic area
    end

    % Mean B_axial = Phi / A_disc
    Bmean = zeros(size(s_query_mm));
    valid = Adisc > 0;
    Bmean(valid) = Phi(valid) ./ Adisc(valid);
    Bmean_mT = abs(Bmean) * 1e3;     % T -> mT

    % Plot from PEAK out to block end (drop the short rising stub before peak).
    valid_s = s_query_mm >= s_tip;
    s_full     = s_query_mm(valid_s);
    Bmean_full = Bmean_mT(valid_s);
    [B_pk, k_pk] = max(Bmean_full);
    s_plot     = s_full(k_pk:end);
    Bmean_plot = Bmean_full(k_pk:end);

    plot(s_plot, Bmean_plot, 'b-', 'LineWidth', 2.2);
    plot(s_plot(1), B_pk, 'bo', 'MarkerSize', 8, 'MarkerFaceColor','b');

    ymax = B_pk * 1.18; ymin = -ymax*0.03;
    plot([0 0],           [ymin ymax], 'k--',                          'LineWidth', 1.4);
    plot([s_tip s_tip],   [ymin ymax], '--', 'Color', [0.92 0.15 0.15], 'LineWidth', 1.4);
    plot([s_base s_base], [ymin ymax], '--', 'Color', [0.15 0.55 0.85], 'LineWidth', 1.4);

    text(0+0.3,      ymax*0.95, 'WP',      'FontSize', 11, 'FontWeight','bold');
    text(s_tip+0.3,  ymax*0.83, 'P2 tip',  'FontSize', 10, ...
         'Color', [0.92 0.15 0.15], 'FontWeight','bold');
    text(s_base+0.3, ymax*0.83, 'cone end','FontSize', 10, ...
         'Color', [0.15 0.55 0.85], 'FontWeight','bold');

    text(s_plot(1)+1.0, B_pk+ymax*0.04, ...
         sprintf('\\langle B_{axial} \\rangle_{max} = %.1f mT @ s = %.2f mm', B_pk, s_plot(1)), ...
         'FontSize', 11, 'FontWeight','bold', 'Color','b', 'Interpreter','tex', ...
         'HorizontalAlignment','left');

    grid on; box on;
    xlim([-2, s_query_mm(end)]); ylim([ymin ymax]);
    xlabel('s [mm]  (P2 pole axis arc-length, WP \rightarrow tip \rightarrow cone end \rightarrow cylinder end)', ...
           'FontSize', 11, 'Interpreter','tex');
    ylabel('\langle B_{axial} \rangle (s) = \Phi(s) / A_{disc}(s)  [mT]', ...
           'FontSize', 11, 'Interpreter','tex');
    title('P2 mean axial flux density  (steel-only interp, full disc with R(s) taper)', ...
          'FontSize', 12, 'Interpreter','tex');

    %% --- Save ---
    out_path = fullfile(out_dir, 'P2_circuit_with_Bdensity_steelonly.png');
    exportgraphics(fig, out_path, 'Resolution', 250);
    fprintf('Saved: %s\n', out_path);
    fprintf('  <B_axial>_max = %.2f mT @ s = %.2f mm\n', B_pk, s_plot(k_pk));
end


function is_steel = mask_p2_steel(x_m, y_m, z_m, tip_apdl, pole_axis, ...
                                   R_tip_m, R_base_m, L_cone_m, L_pole_total_m, ...
                                   yoke_x_m, yoke_z_m, block_x_m, block_z_m, ...
                                   post_xy_m, post_R_m, post_z_m)
    rel_x = x_m - tip_apdl(1);
    rel_y = y_m - tip_apdl(2);
    rel_z = z_m - tip_apdl(3);
    proj = rel_x*pole_axis(1) + rel_y*pole_axis(2) + rel_z*pole_axis(3);

    perp_x = rel_x - proj*pole_axis(1);
    perp_y = rel_y - proj*pole_axis(2);
    perp_z = rel_z - proj*pole_axis(3);
    perp_dist = sqrt(perp_x.^2 + perp_y.^2 + perp_z.^2);

    in_x_cone = (proj >= 0) & (proj <= L_cone_m);
    R_at_proj = R_tip_m + (proj/L_cone_m) * (R_base_m - R_tip_m);
    in_cone = in_x_cone & (perp_dist <= R_at_proj);

    in_x_cyl = (proj > L_cone_m) & (proj <= L_pole_total_m);
    in_cyl = in_x_cyl & (perp_dist <= R_base_m);

    in_yoke = (x_m >= yoke_x_m(1)) & (x_m <= yoke_x_m(2)) & ...
              (z_m >= yoke_z_m(1)) & (z_m <= yoke_z_m(2)) & ...
              (abs(y_m) <= 11e-3);

    in_block = (x_m >= block_x_m(1)) & (x_m <= block_x_m(2)) & ...
               (z_m >= block_z_m(1)) & (z_m <= block_z_m(2)) & ...
               (abs(y_m) <= 11e-3);

    r_from_post = sqrt((x_m - post_xy_m(1)).^2 + (y_m - post_xy_m(2)).^2);
    in_post = (r_from_post <= post_R_m) & (z_m >= post_z_m(1)) & (z_m <= post_z_m(2));

    is_steel = in_cone | in_cyl | in_yoke | in_block | in_post;
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
