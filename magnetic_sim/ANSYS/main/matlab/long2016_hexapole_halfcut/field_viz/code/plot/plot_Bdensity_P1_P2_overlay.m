function plot_Bdensity_P1_P2_overlay()
% PLOT_BDENSITY_P1_P2_OVERLAY
% Overlay P1 mean axial flux density <B_x>(x) and P2 mean axial flux
% density <B_axial>(s) on a single panel for direct comparison.
%
%   <B>(x or s) = Phi(x or s) / A_disc(x or s)
%       P1:  A_disc = pi R(x)^2 / 2   (half-disc, halfcut)
%       P2:  A_disc = pi R(s)^2       (full disc)
%
% Each curve is plotted from PEAK out to block end (drop short rising stub).
% Both abscissae start at WP (=0) so they share a single x-axis labelled
% "x or s [mm]".

    %% --- Paths ---
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    cnst = mt_constants();

    res_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut';
    out_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    %% --- Shared geometry ---
    R_norm_xy = cnst.R_norm_xy * 1e3;
    R_norm_mm = cnst.R_norm    * 1e3;
    POLE_R    = cnst.POLE_R    * 1e3;
    POLE_LEN  = cnst.POLE_CONE_LEN * 1e3;
    POLE_TIP_R_mm = cnst.POLE_TIP_R * 1e3;
    SPH_OFST  = cnst.SPH_OFST  * 1e3;
    Z_CONE    = -13.0;
    X0_LOW    = 47.5;
    x_end_p1  = X0_LOW - 10;
    PROT_R    = 5; PROT_H = 7;

    %% =================== P1 <B_x>(x) ===================
    fprintf('--- P1 (coil1) ---\n');
    d1 = import_ansys_data(fullfile(res_root, 'coil1', 'standard'), 'all', 'coil1');

    x_tip_p1  = R_norm_xy;
    x_base_p1 = x_tip_p1 + POLE_LEN;
    is_steel_p1 = mask_p1_steel(d1.x, d1.y, d1.z, x_tip_p1, x_base_p1, x_end_p1, ...
                                 POLE_R, POLE_TIP_R_mm, Z_CONE, ...
                                 X0_LOW, PROT_R, PROT_H);
    fprintf('  Steel nodes: %d (%.1f%%)\n', sum(is_steel_p1), 100*sum(is_steel_p1)/length(d1.x));
    F_bx_p1 = scatteredInterpolant(d1.x(is_steel_p1), d1.y(is_steel_p1), ...
                                    d1.z(is_steel_p1), d1.bx(is_steel_p1), 'linear', 'none');

    R_func_p1 = @(x) clamp_R(x, x_tip_p1, x_base_p1, POLE_TIP_R_mm, POLE_R);
    x_query = linspace(-2, 37.5, 280);
    Phi_p1   = zeros(size(x_query));
    Adisc_p1 = zeros(size(x_query));
    N_grid   = 81;
    for k = 1:length(x_query)
        xq = x_query(k); R_mm = R_func_p1(xq);
        y_c = linspace(-R_mm, R_mm, N_grid);
        z_c = linspace(Z_CONE-R_mm, Z_CONE, N_grid);
        [Yg, Zg2] = meshgrid(y_c, z_c);
        in_disc = (Yg.^2 + (Zg2-Z_CONE).^2) <= R_mm^2 & Zg2 <= Z_CONE;
        dy = (y_c(2)-y_c(1))*1e-3; dz = (z_c(2)-z_c(1))*1e-3; dA = dy*dz;
        Bx_q = F_bx_p1(xq*1e-3*ones(size(Yg)), Yg*1e-3, Zg2*1e-3);
        Bx_q(~in_disc) = 0; Bx_q(isnan(Bx_q)) = 0;
        Phi_p1(k)   = sum(Bx_q(:)) * dA;
        Adisc_p1(k) = pi * (R_mm*1e-3)^2 / 2;
    end
    Bmean_p1 = zeros(size(x_query));
    valid = Adisc_p1 > 0;
    Bmean_p1(valid) = abs(Phi_p1(valid)) ./ Adisc_p1(valid) * 1e3;   % mT

    valid_x = x_query >= x_tip_p1;
    x_full_p1 = x_query(valid_x);
    B_full_p1 = Bmean_p1(valid_x);
    [B_pk_p1, k_pk_p1] = max(B_full_p1);
    x_p1_plot = x_full_p1(k_pk_p1:end);
    B_p1_plot = B_full_p1(k_pk_p1:end);
    fprintf('  <B_x>_max_P1 = %.2f mT @ x = %.2f mm\n', B_pk_p1, x_p1_plot(1));

    %% =================== P2 <B_axial>(s) ===================
    fprintf('--- P2 (coil5) ---\n');
    d2 = import_ansys_data(fullfile(res_root, 'coil5', 'standard'), 'all', 'coil5');

    bbar_mat = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\bs_matrix\Bbar_S_4p572.mat';
    sign_p2 = +1;
    if exist(bbar_mat, 'file')
        Bs = load(bbar_mat, 'col_sign');
        sign_p2 = Bs.col_sign(2);
    end
    fprintf('  P2 sign: %+d\n', sign_p2);
    if sign_p2 < 0
        d2.bx = -d2.bx; d2.by = -d2.by; d2.bz = -d2.bz;
    end

    inc      = cnst.upper_incline;
    theta_p2 = pi;
    pole_axis = [cos(inc)*cos(theta_p2); cos(inc)*sin(theta_p2); sin(inc)];
    up_unn    = [0;0;1] - sin(inc)*pole_axis;
    up_hat    = up_unn/norm(up_unn);
    side_hat  = cross(pole_axis, up_hat); side_hat = side_hat/norm(side_hat);

    YOKE_H = 2;
    tip_z_upper = -PROT_H - 6 + cnst.R_norm_z*1e3*2;
    tip_apdl = [R_norm_xy*cos(theta_p2)*1e-3; 0; tip_z_upper*1e-3];

    Z0_HIGH = YOKE_H + PROT_H;
    X0_HIGH = 47.5 * cos(theta_p2);
    POST_Z1 = YOKE_H; POST_Z2 = Z0_HIGH;
    BLOCK_Z1 = Z0_HIGH; BLOCK_Z2 = Z0_HIGH + 10;
    BLOCK_X1 = X0_HIGH - 16.5; BLOCK_X2 = X0_HIGH + 8.5;
    POST_R_mm = PROT_H/2;

    end_upper_xz = [-36; YOKE_H + PROT_H + 5];
    pa_xz        = [pole_axis(1); pole_axis(3)];
    tip_xz       = [R_norm_xy*cos(theta_p2); tip_z_upper];
    cone_end_xz  = tip_xz + POLE_LEN * pa_xz;
    L_cyl_mm    = norm(end_upper_xz - cone_end_xz);
    L_pole_tot  = POLE_LEN + L_cyl_mm;

    is_steel_p2 = mask_p2_steel(d2.x, d2.y, d2.z, tip_apdl, pole_axis, ...
                                 POLE_TIP_R_mm*1e-3, POLE_R*1e-3, ...
                                 POLE_LEN*1e-3, L_pole_tot*1e-3, ...
                                 [BLOCK_X1-2, BLOCK_X2+2]*1e-3, [0, YOKE_H]*1e-3, ...
                                 [BLOCK_X1, BLOCK_X2]*1e-3, [BLOCK_Z1, BLOCK_Z2]*1e-3, ...
                                 [X0_HIGH; 0]*1e-3, POST_R_mm*1e-3, [POST_Z1, POST_Z2]*1e-3);
    fprintf('  Steel nodes: %d (%.1f%%)\n', sum(is_steel_p2), 100*sum(is_steel_p2)/length(d2.x));

    F_bx_p2 = scatteredInterpolant(d2.x(is_steel_p2), d2.y(is_steel_p2), ...
                                    d2.z(is_steel_p2), d2.bx(is_steel_p2), 'linear','none');
    F_by_p2 = scatteredInterpolant(d2.x(is_steel_p2), d2.y(is_steel_p2), ...
                                    d2.z(is_steel_p2), d2.by(is_steel_p2), 'linear','none');
    F_bz_p2 = scatteredInterpolant(d2.x(is_steel_p2), d2.y(is_steel_p2), ...
                                    d2.z(is_steel_p2), d2.bz(is_steel_p2), 'linear','none');

    s_tip_p2 = R_norm_mm; s_base_p2 = R_norm_mm + POLE_LEN;
    R_func_p2 = @(s) clamp_R(s, s_tip_p2, s_base_p2, POLE_TIP_R_mm, POLE_R);
    WP_apdl = [0;0;SPH_OFST*1e-3];

    s_query = x_query;
    Phi_p2   = zeros(size(s_query));
    Adisc_p2 = zeros(size(s_query));
    for k = 1:length(s_query)
        s = s_query(k); R_mm = R_func_p2(s);
        center_apdl = WP_apdl + (s*1e-3) * pole_axis;
        u_c = linspace(-R_mm, R_mm, N_grid); v_c = linspace(-R_mm, R_mm, N_grid);
        [Ug, Vg] = meshgrid(u_c, v_c);
        in_disc = (Ug.^2 + Vg.^2) <= R_mm^2;
        du = (u_c(2)-u_c(1))*1e-3; dv = (v_c(2)-v_c(1))*1e-3; dA = du*dv;
        Xq = center_apdl(1) + (Ug*1e-3)*up_hat(1) + (Vg*1e-3)*side_hat(1);
        Yq = center_apdl(2) + (Ug*1e-3)*up_hat(2) + (Vg*1e-3)*side_hat(2);
        Zq = center_apdl(3) + (Ug*1e-3)*up_hat(3) + (Vg*1e-3)*side_hat(3);
        Bxi = F_bx_p2(Xq(:), Yq(:), Zq(:));
        Byi = F_by_p2(Xq(:), Yq(:), Zq(:));
        Bzi = F_bz_p2(Xq(:), Yq(:), Zq(:));
        B_ax = Bxi*pole_axis(1) + Byi*pole_axis(2) + Bzi*pole_axis(3);
        B_ax = reshape(B_ax, N_grid, N_grid);
        B_ax(~in_disc) = 0; B_ax(isnan(B_ax)) = 0;
        Phi_p2(k)   = sum(B_ax(in_disc)) * dA;
        Adisc_p2(k) = pi * (R_mm*1e-3)^2;
    end
    Bmean_p2 = zeros(size(s_query));
    valid = Adisc_p2 > 0;
    Bmean_p2(valid) = abs(Phi_p2(valid)) ./ Adisc_p2(valid) * 1e3;   % mT

    valid_s = s_query >= s_tip_p2;
    s_full_p2 = s_query(valid_s);
    B_full_p2 = Bmean_p2(valid_s);
    [B_pk_p2, k_pk_p2] = max(B_full_p2);
    s_p2_plot = s_full_p2(k_pk_p2:end);
    B_p2_plot = B_full_p2(k_pk_p2:end);
    fprintf('  <B_axial>_max_P2 = %.2f mT @ s = %.2f mm\n', B_pk_p2, s_p2_plot(1));

    %% =================== Plot overlay ===================
    fig = figure('Position', [80 80 1500 700], 'Color', 'w');
    ax  = axes('Position', [0.08 0.13 0.88 0.78]);
    hold on;

    h1 = plot(x_p1_plot, B_p1_plot, '-',  'Color', [0.10 0.40 0.85], 'LineWidth', 2.4);
    h2 = plot(s_p2_plot, B_p2_plot, '--', 'Color', [0.85 0.20 0.20], 'LineWidth', 2.4);

    plot(x_p1_plot(1), B_pk_p1, 'o', 'Color',[0.10 0.40 0.85], 'MarkerFaceColor',[0.10 0.40 0.85], 'MarkerSize', 8);
    plot(s_p2_plot(1), B_pk_p2, 's', 'Color',[0.85 0.20 0.20], 'MarkerFaceColor',[0.85 0.20 0.20], 'MarkerSize', 9);

    ymax_all = max(B_pk_p1, B_pk_p2) * 1.18;
    ymin_all = -ymax_all*0.03;

    % Reference verticals (use P1 cone geometry — close to P2 too)
    plot([0 0],                 [ymin_all ymax_all], 'k--', 'LineWidth', 1.3);
    plot([x_tip_p1 x_tip_p1],   [ymin_all ymax_all], '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0);
    plot([x_base_p1 x_base_p1], [ymin_all ymax_all], '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0);

    text(0+0.3,         ymax_all*0.96, 'WP',       'FontSize', 11, 'FontWeight','bold');
    text(x_tip_p1+0.3,  ymax_all*0.88, 'tip',      'FontSize', 10, 'Color', [0.40 0.40 0.40], 'FontWeight','bold');
    text(x_base_p1+0.3, ymax_all*0.82, 'cone end', 'FontSize', 10, 'Color', [0.40 0.40 0.40], 'FontWeight','bold');

    % Peak labels in upper-left empty zone (curves both fall fast after peak,
    % so any x ~ 5-10 mm has plenty of vertical space above)
    text(5.0, ymax_all*0.95, ...
         sprintf('\\langle B_x \\rangle_{max}^{P1} = %.1f mT @ x = %.2f mm', B_pk_p1, x_p1_plot(1)), ...
         'FontSize', 11, 'FontWeight','bold', 'Color', [0.10 0.40 0.85], 'Interpreter','tex', ...
         'HorizontalAlignment','left');
    text(5.0, ymax_all*0.87, ...
         sprintf('\\langle B_{axial} \\rangle_{max}^{P2} = %.1f mT @ s = %.2f mm', B_pk_p2, s_p2_plot(1)), ...
         'FontSize', 11, 'FontWeight','bold', 'Color', [0.85 0.20 0.20], 'Interpreter','tex', ...
         'HorizontalAlignment','left');

    legend([h1 h2], ...
        {'P1: lower halfcut, half-disc, axis along +x', ...
         'P2: upper full cone tilted +36.59°, full disc, axis along s'}, ...
        'Location', 'northeast', 'FontSize', 11, 'Box', 'on');

    grid on; box on;
    xlim([-2, x_query(end)]);
    ylim([ymin_all ymax_all]);
    xlabel('distance along pole axis [mm]  (WP \rightarrow tip \rightarrow cone end \rightarrow block end)', ...
           'FontSize', 12, 'Interpreter','tex');
    ylabel('\langle B \rangle = \Phi / A_{disc}  [mT]', 'FontSize', 12, 'Interpreter','tex');
    title('P1 vs P2 mean axial flux density  (steel-only interp, baseline, peak \rightarrow block end)', ...
          'FontSize', 13, 'Interpreter','tex');

    %% --- Save ---
    out_path = fullfile(out_dir, 'Bdensity_P1_P2_overlay.png');
    exportgraphics(fig, out_path, 'Resolution', 250);
    fprintf('Saved: %s\n', out_path);
end


function is_steel = mask_p1_steel(x_m, y_m, z_m, x_tip_mm, x_base_mm, x_end_mm, ...
                                   POLE_R_mm, R_tip_mm, Z_CONE_mm, ...
                                   X0_LOW_mm, PROT_R_mm, PROT_H_mm)
    x = x_m*1e3; y = y_m*1e3; z = z_m*1e3;
    r_axis = sqrt(y.^2 + (z - Z_CONE_mm).^2);
    in_x_cone = (x >= x_tip_mm) & (x <= x_base_mm);
    R_at_x = R_tip_mm + (x - x_tip_mm)/(x_base_mm - x_tip_mm) * (POLE_R_mm - R_tip_mm);
    in_cone = in_x_cone & (r_axis <= R_at_x) & (z <= Z_CONE_mm);
    in_x_cyl = (x >= x_base_mm) & (x <= x_end_mm);
    in_cyl = in_x_cyl & (r_axis <= POLE_R_mm) & (z <= Z_CONE_mm);
    in_yoke_lower = (x >= x_end_mm) & (x <= X0_LOW_mm + 11) & ...
                    (z >= -(PROT_H_mm + 10)) & (z <= -PROT_H_mm) & (abs(y) <= 11);
    POST_X1_mm = X0_LOW_mm - PROT_R_mm;
    in_yoke_step = (x >= x_end_mm) & (x <= POST_X1_mm) & ...
                   (z >= -(PROT_H_mm + 2)) & (z <= -PROT_H_mm) & (abs(y) <= 11);
    r_post = sqrt((x - X0_LOW_mm).^2 + y.^2);
    in_post = (r_post <= PROT_R_mm) & (z >= -PROT_H_mm) & (z <= 0);
    is_steel = in_cone | in_cyl | in_yoke_lower | in_yoke_step | in_post;
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
              (z_m >= yoke_z_m(1)) & (z_m <= yoke_z_m(2)) & (abs(y_m) <= 11e-3);
    in_block = (x_m >= block_x_m(1)) & (x_m <= block_x_m(2)) & ...
               (z_m >= block_z_m(1)) & (z_m <= block_z_m(2)) & (abs(y_m) <= 11e-3);
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
