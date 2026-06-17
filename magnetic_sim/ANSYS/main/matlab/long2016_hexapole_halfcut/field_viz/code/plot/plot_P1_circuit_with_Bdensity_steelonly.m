function plot_P1_circuit_with_Bdensity_steelonly()
% PLOT_P1_CIRCUIT_WITH_BDENSITY_STEELONLY
% Same layout as plot_P1_circuit_with_flux_steelonly.m, BUT the bottom
% subplot shows the mean axial flux DENSITY <B_x>(x) = Phi(x) / A_disc(x)
% instead of the integrated flux Phi(x).
%
% Definition (parallel to the flux script):
%   <B_x>(x) = ( sum_{(i,j) in D(x)} B_x(x,y_i,z_j) * dA ) / A_disc(x)
%   A_disc(x) = pi * R(x)^2 / 2     (half-disc area)
%
% Steel-only FEM filter for the Phi numerator (same as flux script).
% Top side view uses the FULL interpolant (so the field outside steel still
% renders correctly for visualization).

    %% --- Paths ---
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\hung\analysis\core');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\hexapole-long2016\analysis');
    cnst = mt_constants();

    res_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\coil1\standard';
    out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    %% --- Geometry (mm, APDL) ---
    R_norm_xy = cnst.R_norm_xy * 1e3;
    POLE_R    = cnst.POLE_R    * 1e3;
    POLE_LEN  = cnst.POLE_CONE_LEN * 1e3;
    SPH_OFST  = cnst.SPH_OFST  * 1e3;
    Z_CONE    = -13.0;
    x_tip     = R_norm_xy;
    x_base    = x_tip + POLE_LEN;
    X0_LOW    = 47.5;
    x_end     = X0_LOW - 10;
    PROT_R    = 5; PROT_H = 7;
    POST_X1   = X0_LOW - PROT_R;
    POST_X2   = X0_LOW + PROT_R;
    COIL_IN_R = 5; COIL_OUT_R = 8;

    %% --- Load FEM ---
    fprintf('Loading coil1 baseline...\n');
    d = import_ansys_data(res_dir, 'all', 'coil1');
    n_all = length(d.x);
    fprintf('  Total FEM nodes: %d\n', n_all);

    %% --- Filter to STEEL only ---
    is_steel = mask_p1_steel(d.x, d.y, d.z, x_tip, x_base, x_end, ...
                              POLE_R, cnst.POLE_TIP_R*1e3, Z_CONE, ...
                              X0_LOW, PROT_R, PROT_H);
    n_steel = sum(is_steel);
    fprintf('  Steel nodes (filtered): %d (=%.1f%%)\n', n_steel, 100*n_steel/n_all);

    d_st.x  = d.x(is_steel);   d_st.y  = d.y(is_steel);   d_st.z  = d.z(is_steel);
    d_st.bx = d.bx(is_steel);  d_st.bz = d.bz(is_steel);

    %% --- Build TWO interpolants ---
    fprintf('Building scatteredInterpolant (full, for side view)...\n');
    F_bx_full = scatteredInterpolant(d.x, d.y, d.z, d.bx, 'linear', 'none');
    F_bz_full = scatteredInterpolant(d.x, d.y, d.z, d.bz, 'linear', 'none');
    fprintf('Building scatteredInterpolant (steel-only, for <B_x> integration)...\n');
    F_bx = scatteredInterpolant(d_st.x, d_st.y, d_st.z, d_st.bx, 'linear', 'none');

    %% --- Figure ---
    fig = figure('Position', [50 30 1500 900], 'Color', 'w');

    %% --- Top: side view + B vectors (identical to flux script) ---
    ax1 = subplot('Position', [0.07 0.55 0.86 0.40]);
    hold on;

    plot([x_tip, x_base], [Z_CONE, Z_CONE-POLE_R], 'k-', 'LineWidth', 2);
    plot([x_base, x_end], [Z_CONE-POLE_R, Z_CONE-POLE_R], 'k-', 'LineWidth', 2);
    plot([x_end, x_tip], [Z_CONE, Z_CONE], 'k-', 'LineWidth', 2);
    fill([x_tip, x_base, x_end, x_tip], [Z_CONE, Z_CONE-POLE_R, Z_CONE-POLE_R, Z_CONE], ...
         [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.7);

    fill([x_end, X0_LOW+11, X0_LOW+11, x_end], ...
         [-PROT_H-10, -PROT_H-10, -PROT_H, -PROT_H], ...
         [0.85 0.85 0.85], 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.7);
    fill([x_end, POST_X1, POST_X1, x_end], ...
         [-PROT_H, -PROT_H, -PROT_H-2, -PROT_H-2], ...
         [0.85 0.85 0.85], 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.7);

    fill([POST_X1, POST_X2, POST_X2, POST_X1], ...
         [-PROT_H, -PROT_H, 0, 0], [0.85 0.85 0.85], 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.7);

    coilL_x = [X0_LOW-COIL_OUT_R, X0_LOW-COIL_IN_R];
    coilR_x = [X0_LOW+COIL_IN_R, X0_LOW+COIL_OUT_R];
    fill([coilL_x(1) coilL_x(2) coilL_x(2) coilL_x(1)], ...
         [-PROT_H -PROT_H 0 0], [0.93 0.55 0.10], 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.55);
    fill([coilR_x(1) coilR_x(2) coilR_x(2) coilR_x(1)], ...
         [-PROT_H -PROT_H 0 0], [0.93 0.55 0.10], 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.55);
    text(mean(coilL_x), -PROT_H/2, char(8857), 'FontSize', 16, 'HorizontalAlignment', 'center', 'FontWeight','bold');
    text(mean(coilR_x), -PROT_H/2, char(8855), 'FontSize', 16, 'HorizontalAlignment', 'center', 'FontWeight','bold');

    grid_x = 130; grid_z = 38;
    xc = linspace(-3, 60, grid_x);  zc = linspace(-18, 2, grid_z);
    [Xg, Zg] = meshgrid(xc, zc); Xg = Xg(:); Zg = Zg(:);
    Bx_q = F_bx_full(Xg*1e-3, zeros(size(Xg)), Zg*1e-3);
    Bz_q = F_bz_full(Xg*1e-3, zeros(size(Xg)), Zg*1e-3);
    keep = ~isnan(Bx_q) & ~isnan(Bz_q);
    Xg = Xg(keep); Zg = Zg(keep); Bx_q = Bx_q(keep); Bz_q = Bz_q(keep);
    Bmag = sqrt(Bx_q.^2 + Bz_q.^2);

    b_clip = quantile(Bmag, 0.92);
    dx_grid = (60 - (-3)) / grid_x;
    arrow_max = 0.85 * dx_grid;
    b_norm = max(Bmag, eps);
    b_eff  = min(Bmag, b_clip);
    scale_per = (arrow_max / b_clip) .* (b_eff ./ b_norm);
    bxq = Bx_q .* scale_per;
    bzq = Bz_q .* scale_per;

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

    plot(0, SPH_OFST, 'k+', 'MarkerSize', 14, 'LineWidth', 2);
    text(-0.5, SPH_OFST+0.5, 'WP', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','right');
    plot([-2, x_end], [Z_CONE, Z_CONE], '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);

    colormap(turbo); clim([0 b_clip]);
    cb = colorbar('Position', [0.945 0.55 0.012 0.40]);
    ylabel(cb, sprintf('|B| [T] (clip 92pct = %.2f T)', b_clip), 'FontSize', 10);

    axis equal; grid on; box on;
    xlim([-3, 60]); ylim([-18, 2]);
    xlabel('x [mm]', 'FontSize', 11);
    ylabel('z [mm]', 'FontSize', 11);
    title('P1 magnetic flux side view  (steel-only interpolation)', 'FontSize', 13);

    %% --- Bottom: mean axial flux DENSITY <B_x>(x) ---
    ax2 = subplot('Position', [0.07 0.08 0.86 0.36]);
    hold on;

    POLE_TIP_R_mm = cnst.POLE_TIP_R * 1e3;
    R_func = @(x) clamp_R(x, x_tip, x_base, POLE_TIP_R_mm, POLE_R);

    x_query_mm = linspace(-2, 37.5, 280);
    Phi   = zeros(size(x_query_mm));
    Adisc = zeros(size(x_query_mm));   % half-disc area [m^2]
    N_grid = 81;
    for k = 1:length(x_query_mm)
        x_q  = x_query_mm(k);
        R_mm = R_func(x_q);
        y_c  = linspace(-R_mm, R_mm, N_grid);
        z_c  = linspace(Z_CONE-R_mm, Z_CONE, N_grid);
        [Yg, Zg2] = meshgrid(y_c, z_c);
        in_disc = (Yg.^2 + (Zg2-Z_CONE).^2) <= R_mm^2 & Zg2 <= Z_CONE;
        dy = (y_c(2)-y_c(1))*1e-3; dz = (z_c(2)-z_c(1))*1e-3; dA = dy*dz;
        Bx_q = F_bx(x_q*1e-3*ones(size(Yg)), Yg*1e-3, Zg2*1e-3);
        Bx_q(~in_disc) = 0; Bx_q(isnan(Bx_q)) = 0;
        Phi(k)   = sum(Bx_q(:)) * dA;
        Adisc(k) = pi * (R_mm*1e-3)^2 / 2;  % half-disc analytic area
    end

    % Mean B_x = Phi / A_disc.  Guard against A_disc -> 0 outside cone (WP side).
    Bmean = zeros(size(x_query_mm));
    valid = Adisc > 0;
    Bmean(valid) = Phi(valid) ./ Adisc(valid);
    Bmean_mT = abs(Bmean) * 1e3;     % convert T -> mT

    % Plot from PEAK out to block end (drop the short rising stub before peak).
    % "從 block 端出發，畫到最大值即可" — keep the long tail [x_pk, x_end].
    valid_x = x_query_mm >= x_tip;
    x_full     = x_query_mm(valid_x);
    Bmean_full = Bmean_mT(valid_x);
    [B_pk, k_pk] = max(Bmean_full);
    x_plot     = x_full(k_pk:end);
    Bmean_plot = Bmean_full(k_pk:end);

    plot(x_plot, Bmean_plot, 'b-', 'LineWidth', 2.2);
    plot(x_plot(1), B_pk, 'bo', 'MarkerSize', 8, 'MarkerFaceColor','b');

    % Vertical reference lines
    ymax = B_pk * 1.18; ymin = -ymax*0.03;
    plot([0 0],           [ymin ymax], 'k--',                          'LineWidth', 1.4);
    plot([x_tip x_tip],   [ymin ymax], '--', 'Color', [0.92 0.15 0.15], 'LineWidth', 1.4);
    plot([x_base x_base], [ymin ymax], '--', 'Color', [0.15 0.55 0.85], 'LineWidth', 1.4);

    text(0+0.3,      ymax*0.95, 'WP',      'FontSize', 11, 'FontWeight','bold');
    text(x_tip+0.3,  ymax*0.83, 'P1 tip',  'FontSize', 10, ...
         'Color', [0.92 0.15 0.15], 'FontWeight','bold');
    text(x_base+0.3, ymax*0.83, 'cone end','FontSize', 10, ...
         'Color', [0.15 0.55 0.85], 'FontWeight','bold');

    % Peak annotation to the RIGHT of peak marker (peak now at LEFT end of curve)
    text(x_plot(1)+1.0, B_pk+ymax*0.04, ...
         sprintf('\\langle B_x \\rangle_{max} = %.1f mT @ x = %.2f mm', B_pk, x_plot(1)), ...
         'FontSize', 11, 'FontWeight','bold', 'Color','b', 'Interpreter','tex', ...
         'HorizontalAlignment','left');

    grid on; box on;
    xlim([-2, x_query_mm(end)]); ylim([ymin ymax]);
    xlabel('x [mm]  (P1 pole axis, WP \rightarrow tip \rightarrow cone end \rightarrow cylinder end)', ...
           'FontSize', 11, 'Interpreter','tex');
    ylabel('\langle B_x \rangle (x) = \Phi(x) / A_{disc}(x)  [mT]', ...
           'FontSize', 11, 'Interpreter','tex');
    title('P1 mean axial flux density  (steel-only interp,  half-disc with R(x))', ...
          'FontSize', 12, 'Interpreter','tex');

    %% --- Save ---
    out_path = fullfile(out_dir, 'P1_circuit_with_Bdensity_steelonly.png');
    exportgraphics(fig, out_path, 'Resolution', 250);
    fprintf('Saved: %s\n', out_path);
    fprintf('  <B_x>_max = %.2f mT @ x = %.2f mm\n', B_pk, x_plot(k_pk));
end


function is_steel = mask_p1_steel(x_m, y_m, z_m, x_tip_mm, x_base_mm, x_end_mm, ...
                                   POLE_R_mm, R_tip_mm, Z_CONE_mm, ...
                                   X0_LOW_mm, PROT_R_mm, PROT_H_mm)
    x  = x_m * 1e3;  y  = y_m * 1e3;  z  = z_m * 1e3;
    r_axis = sqrt(y.^2 + (z - Z_CONE_mm).^2);

    in_x_cone = (x >= x_tip_mm) & (x <= x_base_mm);
    R_at_x = R_tip_mm + (x - x_tip_mm)/(x_base_mm - x_tip_mm) * (POLE_R_mm - R_tip_mm);
    in_cone = in_x_cone & (r_axis <= R_at_x) & (z <= Z_CONE_mm);

    in_x_cyl = (x >= x_base_mm) & (x <= x_end_mm);
    in_cyl = in_x_cyl & (r_axis <= POLE_R_mm) & (z <= Z_CONE_mm);

    in_yoke_lower = (x >= x_end_mm) & (x <= X0_LOW_mm + 11) & ...
                    (z >= -(PROT_H_mm + 10)) & (z <= -PROT_H_mm) & ...
                    (abs(y) <= 11);

    POST_X1_mm = X0_LOW_mm - PROT_R_mm;
    in_yoke_step = (x >= x_end_mm) & (x <= POST_X1_mm) & ...
                   (z >= -(PROT_H_mm + 2)) & (z <= -PROT_H_mm) & ...
                   (abs(y) <= 11);

    r_post = sqrt((x - X0_LOW_mm).^2 + y.^2);
    in_post = (r_post <= PROT_R_mm) & (z >= -PROT_H_mm) & (z <= 0);

    is_steel = in_cone | in_cyl | in_yoke_lower | in_yoke_step | in_post;
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
