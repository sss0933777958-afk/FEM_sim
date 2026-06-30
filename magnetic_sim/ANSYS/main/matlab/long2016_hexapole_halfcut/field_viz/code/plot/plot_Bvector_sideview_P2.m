function plot_Bvector_sideview_P2()
% PLOT_BVECTOR_SIDEVIEW_P2
%   Side view of P2-excitation B field in mirrored APDL xz frame (y=0).
%   P2 is azimuth 180° upper pole, tilted up by upper_incline ≈ 36.59°.
%   Mirror x' = -x so the cone points to the right (matches reference style).
%
%   Style: turbo cmap, denser/thicker arrows for clearer 磁路;
%   cone OUTLINE ONLY (no fill); sensor as tilted GREEN rectangle aligned
%   with cone slant; no tip marker.

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

    apdl_coil_for_P2 = find(strcmp(cnst.apdl_to_paper, 'P2'));   % = 5
    assert(apdl_coil_for_P2 == 5, 'apdl_to_paper map changed?');

    res_dir = sprintf('G:\\my_workspace\\code\\FEM_sim\\kuo\\results\\long2016_hexapole_halfcut\\coil%d', apdl_coil_for_P2);
    out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    %% --- Geometry constants (mm) ---
    inc        = cnst.upper_incline;            % 36.59 deg in rad
    beta       = atan2(3.0, 15.0);              % cone half-angle, 11.31°
    R_norm_xy  = cnst.R_norm_xy * 1e3;          % 0.408 mm
    R_norm_z   = cnst.R_norm_z  * 1e3;          % 0.289 mm
    SPH_OFST   = cnst.SPH_OFST  * 1e3;          % -12.71 mm
    POLE_R     = cnst.POLE_R    * 1e3;          % 3 mm
    POLE_LEN   = cnst.POLE_CONE_LEN * 1e3;      % 15 mm

    fprintf('P2 upper-pole tilt (above horizontal):  %.2f deg\n', rad2deg(inc));
    fprintf('Cone half-angle β:                     %.2f deg\n', rad2deg(beta));

    %% --- Load coil5 (= P2) FEM data ---
    fprintf('\nLoading coil%d (P2 excitation) field grid...\n', apdl_coil_for_P2);
    d = import_ansys_data(res_dir, 'all', sprintf('coil%d', apdl_coil_for_P2));

    %% --- Per-pole sign correction (read from Bbar_S mat) ---
    %  P2 is an upper pole; raw FEM has B_bar_S(2,2) < 0 (per memory).
    %  Sign-correction flips the entire coil5 B field so flux flows yoke→tip→WP.
    bbar_mat = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\bs_matrix\data\Bbar_S_4p572.mat';
    if exist(bbar_mat, 'file')
        Bs = load(bbar_mat, 'col_sign');
        sign_p2 = Bs.col_sign(2);   % P2 is paper index 2
    else
        warning('Bbar_S_4p572.mat not found; using sign_p2 = -1 (P2 default)');
        sign_p2 = -1;
    end
    fprintf('P2 sign correction factor: %+d\n', sign_p2);
    if sign_p2 < 0
        d.bx = -d.bx;
        d.by = -d.by;
        d.bz = -d.bz;
        fprintf('  → Negated B field (P2 sim has flux tip→yoke, displaying flux yoke→tip→WP)\n');
    end

    %% --- Build 3D scatteredInterpolant (use ALL nodes, not slice) ---
    % This handles sparse outer-air mesh by interpolating from nearby 3D nodes.
    fprintf('Building 3D scatteredInterpolant from %d nodes...\n', length(d.x));
    F_bx = scatteredInterpolant(d.x, d.y, d.z, d.bx, 'linear', 'none');
    F_by = scatteredInterpolant(d.x, d.y, d.z, d.by, 'linear', 'none');
    F_bz = scatteredInterpolant(d.x, d.y, d.z, d.bz, 'linear', 'none');

    %% --- Regular xz grid in y=0 plane (mirrored x') ---
    xlim_v = [-2, 16];
    zlim_v = [-15, 0];
    grid_x = 50; grid_z = 42;       % regular xz grid
    x_centers = linspace(xlim_v(1)+0.18, xlim_v(2)-0.18, grid_x);
    z_centers = linspace(zlim_v(1)+0.18, zlim_v(2)-0.18, grid_z);
    [Xg, Zg] = meshgrid(x_centers, z_centers);
    Xg = Xg(:); Zg = Zg(:);

    % Convert mirrored x' back to APDL x for interpolation: x_apdl = -x'
    x_apdl = -Xg * 1e-3;            % m
    y_apdl = zeros(size(Xg));       % y = 0 plane
    z_apdl =  Zg * 1e-3;            % m

    Bx_g = F_bx(x_apdl, y_apdl, z_apdl);
    By_g = F_by(x_apdl, y_apdl, z_apdl);
    Bz_g = F_bz(x_apdl, y_apdl, z_apdl);

    % Drop NaN cells (outside FEM domain)
    keep = ~isnan(Bx_g) & ~isnan(Bz_g);
    Xg = Xg(keep); Zg = Zg(keep);
    Bx_g = Bx_g(keep); By_g = By_g(keep); Bz_g = Bz_g(keep);

    % Mirror Bx for mirrored-x' frame
    bx_mir_g = -Bx_g;
    bsum_g   = sqrt(Bx_g.^2 + By_g.^2 + Bz_g.^2);

    fprintf('  Regular-grid arrows (after NaN filter): %d\n', length(Xg));

    %% Linear arrow scaling
    arrow_max = 0.45;
    bmax_arr = max(bsum_g);
    scale = arrow_max / bmax_arr;
    bx_q = bx_mir_g * scale;
    bz_q = Bz_g     * scale;
    x_q  = Xg;
    z_q  = Zg;
    bmag_q = bsum_g;

    %% --- Plot ---
    fig = figure('Position', [50 50 1400 1000], 'Color', 'w');
    set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LBL);

    n_bins  = 28;
    edges_b = linspace(0, max(bmag_q), n_bins + 1);
    cmap_n  = turbo(n_bins);
    lw_rng  = [0.70, 2.40];                      % thicker (was 0.50-1.85)

    hold on;
    for k = 1:n_bins
        in_bin = bmag_q >= edges_b(k) & bmag_q < edges_b(k+1);
        if k == n_bins
            in_bin = in_bin | (bmag_q >= edges_b(end));
        end
        if any(in_bin)
            lw = lw_rng(1) + (k-1)/(n_bins-1) * (lw_rng(2)-lw_rng(1));
            quiver(x_q(in_bin), z_q(in_bin), bx_q(in_bin), bz_q(in_bin), 0, ...
                   'Color', cmap_n(k,:), 'LineWidth', lw, 'MaxHeadSize', 0.6);
        end
    end

    %% --- P2 cone OUTLINE only (no fill) ---
    % Tip in mirrored xz: (R_norm_xy, SPH_OFST + R_norm_z)
    T   = [R_norm_xy, SPH_OFST + R_norm_z];

    % Cone slant length = sqrt(15² + 3²) = 15.30 mm
    L_slant = sqrt(POLE_LEN^2 + POLE_R^2);

    % Cone top edge direction (β above pole axis, pole axis at +inc above horiz)
    C_top = T + L_slant * [cos(inc + beta), sin(inc + beta)];
    C_bot = T + L_slant * [cos(inc - beta), sin(inc - beta)];

    plot([T(1), C_top(1)], [T(2), C_top(2)], 'k-', 'LineWidth', LW_STEEL);
    plot([T(1), C_bot(1)], [T(2), C_bot(2)], 'k-', 'LineWidth', LW_STEEL);
    plot([C_top(1), C_bot(1)], [C_top(2), C_bot(2)], 'k-', 'LineWidth', LW_STEEL);

    %% --- WP marker ---
    plot(0, SPH_OFST, 'k+', 'MarkerSize', 16, 'LineWidth', 2.4);
    text(0 + 0.4, SPH_OFST - 0.4, 'WP', 'FontSize', 13, 'FontWeight', 'bold');

    %% --- Surface sensor (GREEN tilted panel aligned with cone slant) ---
    % Sensor position in pole-local (per gen_Vout_Vin_4p572.m, upper-pole 4.572 mm):
    S_ALONG  = 4.572;   % mm along cone surface from tip
    S_OFFSET = 0.41;    % mm perpendicular to cone surface

    %  s_axial_pl = S_ALONG·cos(β) − S_OFFSET·sin(β);
    %  s_up_pl    = S_ALONG·sin(β) + S_OFFSET·cos(β);
    %  Map pole-local → mirrored xz: x' = axial·cos(inc) + up·(-sin(inc)),
    %                                z' = axial·sin(inc) + up· cos(inc).
    s_axial = S_ALONG*cos(beta) - S_OFFSET*sin(beta);
    s_up    = S_ALONG*sin(beta) + S_OFFSET*cos(beta);
    sens    = T + s_axial * [cos(inc), sin(inc)] + s_up * [-sin(inc), cos(inc)];

    % Sensor panel: tangent direction (aligned with cone slant, top edge direction)
    %               normal  direction (perpendicular outward, away from cone)
    slant_tangent = [cos(inc + beta), sin(inc + beta)];
    slant_normal  = [-sin(inc + beta), cos(inc + beta)];   % perp +90°

    panel_L = 1.0;     % mm panel length (along tangent)
    panel_W = 0.20;    % mm panel thickness (along normal)
    hw = panel_W/2;
    hl = panel_L/2;
    P1 = sens + (+hl)*slant_tangent + (+hw)*slant_normal;
    P2 = sens + (-hl)*slant_tangent + (+hw)*slant_normal;
    P3 = sens + (-hl)*slant_tangent + (-hw)*slant_normal;
    P4 = sens + (+hl)*slant_tangent + (-hw)*slant_normal;
    fill([P1(1) P2(1) P3(1) P4(1)], [P1(2) P2(2) P3(2) P4(2)], ...
         [0.20 0.75 0.25], 'EdgeColor', 'k', 'LineWidth', 1.5);

    % Small arrow showing sensor n+ direction (perpendicular outward)
    arrow_n = 0.6;   % mm
    n_tip = sens + arrow_n * slant_normal;
    plot([sens(1), n_tip(1)], [sens(2), n_tip(2)], '-', ...
         'Color', [0.92 0.15 0.15], 'LineWidth', 2.4);
    plot(n_tip(1), n_tip(2), '^', 'MarkerSize', 8, ...
         'MarkerFaceColor', [0.92 0.15 0.15], 'MarkerEdgeColor', [0.92 0.15 0.15]);

    text(sens(1) + 0.8, sens(2) + 0.6, 'B_{surface}', ...
         'FontSize', 13, 'Color', [0.92 0.15 0.15], 'FontWeight', 'bold');
    text(n_tip(1) + 0.2, n_tip(2), 'n_+', ...
         'FontSize', 12, 'Color', [0.92 0.15 0.15], 'FontWeight', 'bold');

    hold off;

    %% --- Colorbar + axes ---
    colormap(turbo);
    clim([0, max(bmag_q)]);
    cb = colorbar;
    ylabel(cb, '|B|  [Tesla]  (3D magnitude)', 'FontSize', FONT_CB);
    set(cb, 'FontSize', FONT_CB);

    axis equal;
    grid on;
    xlim(xlim_v); ylim(zlim_v);
    xlabel('x'' [mm]  (= -x_{APDL}, mirrored so P2 cone points right)', 'FontSize', FONT_LBL);
    ylabel('z [mm]  (APDL frame)', 'FontSize', FONT_LBL);
    title(sprintf(['P2 excitation B-field (sign-corrected: flux yoke→tip→WP) — '...
                   'upper pole tilted +%.2f° (upper incline)'], rad2deg(inc)), ...
          'FontSize', FONT_TTL, 'Interpreter', 'none');

    out_path = fullfile(out_dir, 'Bvector_sideview_P2.png');
    exportgraphics(fig, out_path, 'Resolution', DPI);
    fprintf('\nSaved: %s\n', out_path);
end
