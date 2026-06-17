function plot_P1_fig25c_charge(save_final)
% PLOT_P1_FIG25C_CHARGE
%   P1 magnetic circuit arrow figure, in the style of Long 2016 dissertation
%   Fig. 2.5(c): FEM B-field arrows concentrated near the sharp tip, with
%   the equivalent point magnetic charge ([J] fit, cube +/-40 um) marked.
%
%   View: TOP VIEW - xy plane just below the halfcut surface (z = -13.02 mm,
%         20 um inside the steel). The cone is cut along its full length so
%         the flux convergence toward the tip is visible. The [J] charge
%         (0.16 mm deeper) is projected onto this view.
%   Zoom: concentrated near the tip so flux convergence into the tip is visible.
%
%   Color: paper style blue -> red (weak air field = blue, strong pole flux = red).
%
%   Charge: magnetic_sim/ANSYS/main/MATLAB_data/long2016_hexapole_halfcut/charge_fit/joint_6coil_40um_fit.mat
%           best.pos(:,1) in WP frame -> APDL frame via z + SPH_OFST.
%
%   Usage:
%     plot_P1_fig25c_charge()        % preview (figures dir, *_preview.png)
%     plot_P1_fig25c_charge(true)    % save final (300 dpi)
%
%   Data: coil1 (P1 excitation) baseline = Long2016 verbatim, 0.6 A.

    if nargin < 1, save_final = false; end

    %% --- Style ---
    FONT_NAME  = 'Helvetica';
    FONT_LBL   = 12;
    FONT_TTL   = 13;
    FONT_ANNOT = 12;
    LW_STEEL   = 1.6;
    DPI        = 300;

    COL_LOW    = [0.15 0.25 0.85];   % blue (weak field, air)
    COL_HIGH   = [0.85 0.08 0.08];   % red  (strong flux in pole)
    COL_CHARGE = [0 0 0];            % black dot
    COL_FLUX   = [0.05 0.55 0.10];   % green dashed ellipse
    COL_ARROW  = [0.85 0.08 0.08];   % red annotation arrows

    %% --- Paths ---
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\hung\analysis\core');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\hexapole-long2016\analysis');

    cnst = mt_constants();

    res_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\coil1\standard';
    data_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\MATLAB_data';
    out_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    %% --- Geometry (APDL mm) ---
    R_norm_xy  = cnst.R_norm_xy * 1e3;          % 0.408
    R_norm_z   = cnst.R_norm_z  * 1e3;          % 0.289
    SPH_OFST   = cnst.SPH_OFST  * 1e3;          % -12.71
    POLE_R     = cnst.POLE_R    * 1e3;          % 3
    POLE_TIP_R = cnst.POLE_TIP_R * 1e3;         % 0.04
    POLE_LEN   = cnst.POLE_CONE_LEN * 1e3;      % 15
    Z_CONE     = SPH_OFST - R_norm_z;           % -13  (cone axis / halfcut plane)

    x_tip   = R_norm_xy;                        %  0.408
    x_base  = x_tip + POLE_LEN;                 % 15.408 (cone base)

    %% --- [J] fit charge position (WP frame -> APDL frame) ---
    fit_mat = fullfile(data_dir, 'charge_fit', 'joint_6coil_40um_fit.mat');
    S = load(fit_mat, 'best', 'cube_half');
    chg_wp = S.best.pos(:, 1);                            % P1, metres, WP frame
    chg_x  = chg_wp(1) * 1e3;                             % mm
    chg_y  = chg_wp(2) * 1e3;                             % mm
    chg_z  = chg_wp(3) * 1e3 + SPH_OFST;                  % mm, APDL frame
    fprintf('[J] P1 charge (APDL frame): x = %.4f, y = %.4f, z = %.4f mm\n', chg_x, chg_y, chg_z);
    fprintf('[J] fit mean error: %.2f%% (cube +/-%.0f um)\n', ...
        S.best.mean_err*100, S.cube_half*1e6);

    % Sampling plane: just below the halfcut surface (inside the steel along
    % the full pole length, so tip flux convergence is visible)
    z_plane = Z_CONE - 0.02;                              % -13.02 mm

    %% --- Load coil1 (P1 excitation) ---
    fprintf('Loading coil1 (P1 excitation) data...\n');
    d = import_ansys_data(res_dir, 'all', 'coil1');

    %% --- Per-pole sign correction ---
    bbar_mat = fullfile(data_dir, 'bs_matrix', 'Bbar_S_4p572.mat');
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

    %% --- xy grid (top view, tight zoom on tip to show flux convergence) ---
    xlim_v = [-1.5, 3.5];
    ylim_v = [-1.5, 1.5];
    grid_x = 100;
    grid_y = 60;
    x_centers = linspace(xlim_v(1)+0.1, xlim_v(2)-0.1, grid_x);
    y_centers = linspace(ylim_v(1)+0.1, ylim_v(2)-0.1, grid_y);
    [Xg, Yg] = meshgrid(x_centers, y_centers);
    Xg = Xg(:); Yg = Yg(:);

    Bx_g = F_bx(Xg*1e-3, Yg*1e-3, z_plane*1e-3 * ones(size(Xg)));
    By_g = F_by(Xg*1e-3, Yg*1e-3, z_plane*1e-3 * ones(size(Xg)));

    keep = ~isnan(Bx_g) & ~isnan(By_g);
    Xg = Xg(keep); Yg = Yg(keep);
    Bx_g = Bx_g(keep); By_g = By_g(keep);

    bsum_g = sqrt(Bx_g.^2 + By_g.^2);
    fprintf('  Grid arrows: %d | in-plane |B| range [%.3g, %.3g] T\n', ...
        length(Xg), min(bsum_g), max(bsum_g));

    %% --- Color: paper-style blue -> red by |B| (clipped at 90th pct) ---
    b_clip = quantile(bsum_g, 0.90);

    %% --- Arrow scaling: uniform visual length, color encodes |B| ---
    dx_grid = (xlim_v(2)-xlim_v(1)) / grid_x;
    arrow_max = 0.95 * dx_grid;
    b_norm = max(bsum_g, eps);
    bx_q = Bx_g ./ b_norm * arrow_max;
    by_q = By_g ./ b_norm * arrow_max;

    %% --- Figure ---
    fig = figure('Position', [60 60 1500 900], 'Color', 'w');
    set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LBL);
    ax = axes;
    hold(ax, 'on');

    % Binned blue -> red quiver
    n_bins  = 16;
    edges_b = linspace(0, b_clip, n_bins + 1);
    t_bins  = ((1:n_bins) - 0.5)' / n_bins;
    cmap_rb = (1 - t_bins) * COL_LOW + t_bins * COL_HIGH;
    lw_rng  = [0.55, 1.9];

    for k = 1:n_bins
        if k < n_bins
            in_bin = bsum_g >= edges_b(k) & bsum_g < edges_b(k+1);
        else
            in_bin = bsum_g >= edges_b(k);     % top bin includes clipped tail
        end
        if any(in_bin)
            lw = lw_rng(1) + (k-1)/(n_bins-1) * (lw_rng(2)-lw_rng(1));
            quiver(ax, Xg(in_bin), Yg(in_bin), bx_q(in_bin), by_q(in_bin), 0, ...
                   'Color', cmap_rb(k,:), 'LineWidth', lw, 'MaxHeadSize', 0.6);
        end
    end

    %% --- P1 cone outline (top view projection: symmetric triangle) ---
    %   tip (0.408, 0) -> base (15.408, +/-3); cone extends beyond x-range
    plot(ax, [x_tip, x_base], [+POLE_TIP_R, +POLE_R], 'k-', 'LineWidth', LW_STEEL);
    plot(ax, [x_tip, x_base], [-POLE_TIP_R, -POLE_R], 'k-', 'LineWidth', LW_STEEL);
    plot(ax, [x_tip, x_tip],  [-POLE_TIP_R, +POLE_TIP_R], 'k-', 'LineWidth', LW_STEEL);

    %% --- WP marker ---
    plot(ax, 0, 0, 'k+', 'MarkerSize', 16, 'LineWidth', 2.4);
    text(ax, 0, 0.12, 'WP', 'FontSize', 13, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

    %% --- Magnetic charge marker ([J] fit, projected onto this plane) ---
    plot(ax, chg_x, chg_y, 'o', 'MarkerSize', 16, ...
        'MarkerFaceColor', COL_CHARGE, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    text(ax, chg_x + 0.08, chg_y + 0.12, sprintf('q  (depth %.2f mm)', Z_CONE - chg_z), ...
        'FontSize', 14, 'FontWeight', 'bold', 'FontAngle', 'italic');

    % "Magnetic charge" annotation (black box + red arrow, paper style)
    annot_chg = [-0.5, -1.0];
    plot(ax, [annot_chg(1), chg_x - 0.03], [annot_chg(2), chg_y - 0.06], '-', ...
        'Color', COL_ARROW, 'LineWidth', 2.2);
    fill_arrowhead(ax, chg_x - 0.03, chg_y - 0.06, ...
        atan2(chg_y - annot_chg(2), chg_x - annot_chg(1)), 0.10, COL_ARROW);
    text(ax, annot_chg(1), annot_chg(2) - 0.10, {'Magnetic', 'charge'}, ...
        'FontSize', FONT_ANNOT + 1, 'FontWeight', 'bold', ...
        'EdgeColor', 'k', 'BackgroundColor', 'w', 'Margin', 4, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

    %% --- Labels ---
    text(ax, x_tip - 0.05, -0.15, 'P1 tip', 'FontSize', 11, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'right');

    hold(ax, 'off');

    %% --- Colorbar + axes ---
    colormap(ax, cmap_rb);
    clim(ax, [0, b_clip]);
    cb = colorbar(ax);
    ylabel(cb, sprintf('in-plane |B|  [T]  (clipped at 90th pct = %.2f T)', b_clip), ...
        'FontSize', 11);

    axis(ax, 'equal');
    xlim(ax, xlim_v); ylim(ax, ylim_v);
    set(ax, 'Layer', 'top');
    xlabel(ax, 'x [mm]  (APDL frame, +x = P1 pole axis)', 'FontSize', FONT_LBL);
    ylabel(ax, 'y [mm]', 'FontSize', FONT_LBL);
    title(ax, sprintf(['P1 pole magnetic circuit - top view (xy at z = %.2f mm, ' ...
        '20 um below halfcut surface) - [J] fit cube +/-40 um, mean err %.2f%%'], ...
        z_plane, S.best.mean_err*100), ...
        'FontSize', FONT_TTL, 'Interpreter', 'none');

    %% --- Save ---
    if save_final
        out_path = fullfile(out_dir, 'P1_fig25c_charge.png');
    else
        out_path = fullfile(out_dir, 'P1_fig25c_charge_preview.png');
    end
    exportgraphics(fig, out_path, 'Resolution', DPI);
    fprintf('Saved: %s\n', out_path);
end


function fill_arrowhead(ax, x, y, angle, sz, color)
% Filled triangular arrowhead at (x, y) pointing along `angle` (rad).
    half_w = 0.4 * sz;
    pts = [0, 0; -sz, half_w; -sz, -half_w]';
    R = [cos(angle), -sin(angle); sin(angle), cos(angle)];
    pts = R * pts;
    fill(ax, x + pts(1,:), y + pts(2,:), color, 'EdgeColor', color);
end
