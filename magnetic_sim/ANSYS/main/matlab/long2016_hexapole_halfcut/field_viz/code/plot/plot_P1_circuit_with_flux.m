function plot_P1_circuit_with_flux()
% PLOT_P1_CIRCUIT_WITH_FLUX
% Two-panel: top = P1 side view (APDL xz, y=0) with COIL+BLOCK+POLE+POST
%             and B-field vectors; bottom = axial flux profile.

    %% --- Paths ---
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    cnst = mt_constants();

    res_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data\coil1\standard';
    out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    %% --- Geometry (mm, APDL) ---
    R_norm_xy = cnst.R_norm_xy * 1e3;            % 0.408
    POLE_R    = cnst.POLE_R    * 1e3;            % 3
    POLE_LEN  = cnst.POLE_CONE_LEN * 1e3;        % 15
    SPH_OFST  = cnst.SPH_OFST  * 1e3;            % -12.711
    Z_CONE    = -13.0;                            % halfcut plane
    x_tip     = R_norm_xy;                        % 0.408
    x_base    = x_tip + POLE_LEN;                 % 15.408
    X0_LOW    = 47.5;
    x_end     = X0_LOW - 10;                      % 37.5 (cylinder end)
    PROT_R    = 5; PROT_H = 7;
    POST_X1   = X0_LOW - PROT_R;                  % 42.5
    POST_X2   = X0_LOW + PROT_R;                  % 52.5
    COIL_IN_R = 5; COIL_OUT_R = 8;

    %% --- Load FEM ---
    fprintf('Loading coil1 baseline...\n');
    d = import_ansys_data(res_dir, 'all', 'coil1');
    F_bx = scatteredInterpolant(d.x, d.y, d.z, d.bx, 'linear', 'none');
    F_bz = scatteredInterpolant(d.x, d.y, d.z, d.bz, 'linear', 'none');

    %% --- Figure ---
    fig = figure('Position', [50 30 1500 900], 'Color', 'w');

    %% --- Top: side view + B vectors ---
    ax1 = subplot('Position', [0.07 0.55 0.86 0.40]);
    hold on;

    % Pole (halfcut D-shape outline: tip → cone slant → cylinder bottom → end → halfcut top)
    plot([x_tip, x_base], [Z_CONE, Z_CONE-POLE_R], 'k-', 'LineWidth', 2);     % cone slant
    plot([x_base, x_end], [Z_CONE-POLE_R, Z_CONE-POLE_R], 'k-', 'LineWidth', 2); % cylinder bottom
    plot([x_end, x_tip], [Z_CONE, Z_CONE], 'k-', 'LineWidth', 2);              % halfcut top
    fill([x_tip, x_base, x_end, x_tip], [Z_CONE, Z_CONE-POLE_R, Z_CONE-POLE_R, Z_CONE], ...
         [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.7);

    % Block (combined union of 4 sub-blocks)
    blk_x = [x_end, x_end, POST_X1, POST_X1, POST_X2, X0_LOW+11, X0_LOW+11, x_end];
    blk_z = [Z_CONE-POLE_R, -PROT_H-2, -PROT_H-2, -PROT_H, -PROT_H, -PROT_H, Z_CONE-POLE_R-1, Z_CONE-POLE_R];
    % Simplify: use rectangle for the entire block bounding region
    fill([x_end, X0_LOW+11, X0_LOW+11, x_end], ...
         [-PROT_H-10, -PROT_H-10, -PROT_H, -PROT_H], ...
         [0.85 0.85 0.85], 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.7);
    % Block top step (37.5-42.5 region 8mm; 42.5-58.5 region 10mm)
    fill([x_end, POST_X1, POST_X1, x_end], ...
         [-PROT_H, -PROT_H, -PROT_H-2, -PROT_H-2], ...
         [0.85 0.85 0.85], 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.7);

    % Post (vertical)
    fill([POST_X1, POST_X2, POST_X2, POST_X1], ...
         [-PROT_H, -PROT_H, 0, 0], [0.85 0.85 0.85], 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.7);

    % Coil (orange, 2 sides of post)
    coilL_x = [X0_LOW-COIL_OUT_R, X0_LOW-COIL_IN_R];
    coilR_x = [X0_LOW+COIL_IN_R, X0_LOW+COIL_OUT_R];
    fill([coilL_x(1) coilL_x(2) coilL_x(2) coilL_x(1)], ...
         [-PROT_H -PROT_H 0 0], [0.93 0.55 0.10], 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.55);
    fill([coilR_x(1) coilR_x(2) coilR_x(2) coilR_x(1)], ...
         [-PROT_H -PROT_H 0 0], [0.93 0.55 0.10], 'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceAlpha', 0.55);
    text(mean(coilL_x), -PROT_H/2, char(8857), 'FontSize', 16, 'HorizontalAlignment', 'center', 'FontWeight','bold');
    text(mean(coilR_x), -PROT_H/2, char(8855), 'FontSize', 16, 'HorizontalAlignment', 'center', 'FontWeight','bold');

    % B field quiver
    grid_x = 130; grid_z = 38;
    xc = linspace(-3, 60, grid_x);  zc = linspace(-18, 2, grid_z);
    [Xg, Zg] = meshgrid(xc, zc); Xg = Xg(:); Zg = Zg(:);
    Bx_q = F_bx(Xg*1e-3, zeros(size(Xg)), Zg*1e-3);
    Bz_q = F_bz(Xg*1e-3, zeros(size(Xg)), Zg*1e-3);
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
    title('P1 magnetic flux side view', 'FontSize', 13);

    %% --- Bottom: axial flux profile ---
    ax2 = subplot('Position', [0.07 0.08 0.86 0.36]);
    hold on;

    % Compute flux profile (same as plot_P1_flux_profile.m)
    POLE_TIP_R_mm = cnst.POLE_TIP_R * 1e3;
    cone_len_mm   = POLE_LEN;
    R_func = @(x) clamp_R(x, x_tip, x_base, POLE_TIP_R_mm, POLE_R);

    x_query_mm = linspace(-2, 37.5, 280);
    Phi = zeros(size(x_query_mm));
    N_grid = 81;
    for k = 1:length(x_query_mm)
        x_q = x_query_mm(k);
        R_mm = R_func(x_q);
        y_c = linspace(-R_mm, R_mm, N_grid);
        z_c = linspace(Z_CONE-R_mm, Z_CONE, N_grid);
        [Yg, Zg2] = meshgrid(y_c, z_c);
        in_disc = (Yg.^2 + (Zg2-Z_CONE).^2) <= R_mm^2 & Zg2 <= Z_CONE;
        dy = (y_c(2)-y_c(1))*1e-3; dz = (z_c(2)-z_c(1))*1e-3; dA = dy*dz;
        Bx_q = F_bx(x_q*1e-3*ones(size(Yg)), Yg*1e-3, Zg2*1e-3);
        Bx_q(~in_disc) = 0; Bx_q(isnan(Bx_q)) = 0;
        Phi(k) = sum(Bx_q(:)) * dA;
    end
    Phi_uWb = abs(Phi) * 1e6;

    plot(x_query_mm, Phi_uWb, 'b-', 'LineWidth', 2.2);

    ymax = max(Phi_uWb)*1.10; ymin = -ymax*0.05;
    plot([0 0], [ymin ymax], 'k--', 'LineWidth', 1.4);
    plot([x_tip x_tip], [ymin ymax], '--', 'Color', [0.92 0.15 0.15], 'LineWidth', 1.4);
    plot([x_base x_base], [ymin ymax], '--', 'Color', [0.15 0.55 0.85], 'LineWidth', 1.4);

    text(0+0.3, ymax*0.93, 'WP', 'FontSize', 11, 'FontWeight','bold');
    text(x_tip+0.3, ymax*0.85, 'P1 tip', 'FontSize', 10, 'Color', [0.92 0.15 0.15], 'FontWeight','bold');
    text(x_base+0.3, ymax*0.85, 'cone end', 'FontSize', 10, 'Color', [0.15 0.55 0.85], 'FontWeight','bold');

    [phi_pk, k_pk] = max(Phi_uWb);
    plot(x_query_mm(k_pk), phi_pk, 'bo', 'MarkerSize', 8, 'MarkerFaceColor','b');
    text(x_query_mm(k_pk)-0.6, phi_pk+0.06, ...
         sprintf('\\Phi_{max} = %.2f \\muWb @ x = %.1f mm', phi_pk, x_query_mm(k_pk)), ...
         'FontSize', 11, 'FontWeight','bold', 'Color','b', 'Interpreter','tex', ...
         'HorizontalAlignment','right');

    grid on; box on;
    xlim([-2 37.5]); ylim([ymin ymax]);
    xlabel('x [mm]  (P1 pole axis, WP \rightarrow tip \rightarrow cone end \rightarrow cylinder end)', ...
           'FontSize', 11, 'Interpreter','tex');
    ylabel('|\Phi(x)| = |\int B_x dA| [\muWb]', 'FontSize', 11, 'Interpreter','tex');
    title('P1 axial flux profile (half-disc with R(x), z \leq -13 mm)', 'FontSize', 12, 'Interpreter','tex');

    %% --- Save ---
    out_path = fullfile(out_dir, 'P1_circuit_with_flux.png');
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
