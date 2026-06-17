function plot_Bvector_P1_coil_block_pole()
% PLOT_BVECTOR_P1_COIL_BLOCK_POLE
%   Side view (xz at y=0) of the entire P1 magnetic structure:
%     COIL (SOURC36 ring cross-section, both sides of post)
%   + BLOCK (4 connector blocks between pole-end and post bottom)
%   + POLE (revolved cone + cylinder, halfcut lower half)
%   + POST (vertical cylindrical cross-section, supports coil)
%
%   All geometry traced from MT_Sim_P1.txt verbatim:
%
%   Pole (revolved around z=-13 axis, halfcut keeps z <= -13):
%       tip   (0.408, -13)
%       cone  -> base (15.408, -16)
%       cyl   -> end  (37.5,  -16)
%
%   4 connector blocks (NOT halfcut — span full z range):
%       Block 1: x=[42.5, 58.5], z=[-17, -7]   (biggest, under post)
%       Block 2: x=[42.5, 45.5], z=[-17, -12]
%       Block 3: x=[37.5, 45.5], z=[-17,  -9]
%       Block 4: x=[37.5, 42.5], z=[-17, -12]
%       y=0 union outline (stepped):
%         x=[37.5, 42.5]: z=[-17, -9]   (height 8 mm)
%         x=[42.5, 58.5]: z=[-17, -7]   (height 10 mm)
%
%   Post (CYL4 vertical cylinder at x=47.5, R=5 mm):
%       x=[42.5, 52.5], z=[-7, 0]
%
%   Coil ring (SOURC36, R_in=5, R_out=8, H=7, centered at (47.5, 0, -3.5)):
%       Left  cross-section:  x=[39.5, 42.5], z=[-7, 0]
%       Right cross-section:  x=[52.5, 55.5], z=[-7, 0]
%
%   Data: coil1 smrt 4 (P1 excitation, sign-corrected).

    %% --- Style ---
    FONT_NAME = 'Helvetica';
    FONT_LBL  = 12;
    FONT_TTL  = 13;
    FONT_CB   = 11;
    LW_STEEL  = 2.0;
    LW_COIL   = 1.5;
    LW_HALFCUT = 1.0;
    DPI       = 300;

    COL_COIL    = [0.93 0.55 0.10];
    COL_HALFCUT = [0.45 0.45 0.45];

    %% --- Paths ---
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\hung\analysis\core');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\hexapole-long2016\analysis');

    cnst = mt_constants();

    res_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\coil1\standard';
    out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    %% --- Geometry (APDL mm) — verbatim from MT_Sim_P1.txt ---
    R_norm_xy = cnst.R_norm_xy * 1e3;            % 0.408
    R_norm_z  = cnst.R_norm_z  * 1e3;            % 0.289
    SPH_OFST  = cnst.SPH_OFST  * 1e3;            % -12.71
    POLE_R    = cnst.POLE_R    * 1e3;            % 3
    POLE_LEN  = cnst.POLE_CONE_LEN * 1e3;        % 15
    Z_CONE    = SPH_OFST - R_norm_z;             % -13

    % Pole (revolved cone + cylinder)
    x_tip   = R_norm_xy;                         % 0.408
    x_base  = x_tip + POLE_LEN;                  % 15.408
    X0_LOW  = (42 + 53)/2;                       % 47.5    (YOKE_MID_R)
    x_end   = X0_LOW - 10;                       % 37.5

    % Post
    PROT_R  = 5;
    PROT_H  = 7;
    POST_X1 = X0_LOW - PROT_R;                   % 42.5
    POST_X2 = X0_LOW + PROT_R;                   % 52.5
    POST_Z1 = -PROT_H;                           % -7
    POST_Z2 = 0;

    % 4 connector blocks (line 107-110 of MT_Sim_P1.txt)
    B1 = struct('x1', X0_LOW-5,  'x2', X0_LOW+11,  'z1', -PROT_H-10, 'z2', POST_Z1);   % [42.5, 58.5] x [-17, -7]
    B2 = struct('x1', X0_LOW-5,  'x2', X0_LOW-2,   'z1', -PROT_H-10, 'z2', -PROT_H-5); % [42.5, 45.5] x [-17, -12]
    B3 = struct('x1', X0_LOW-10, 'x2', X0_LOW-2,   'z1', -PROT_H-10, 'z2', -PROT_H-2); % [37.5, 45.5] x [-17, -9]
    B4 = struct('x1', X0_LOW-10, 'x2', X0_LOW-5,   'z1', -PROT_H-10, 'z2', -PROT_H-5); % [37.5, 42.5] x [-17, -12]

    % Coil ring (SOURC36, both sides at y=0 cross-section)
    COIL_IN_R  = 5;
    COIL_OUT_R = 8;
    COIL_Z1    = -PROT_H;
    COIL_Z2    = 0;
    coilL_x = [X0_LOW-COIL_OUT_R, X0_LOW-COIL_IN_R];   % [39.5, 42.5]
    coilR_x = [X0_LOW+COIL_IN_R,  X0_LOW+COIL_OUT_R];  % [52.5, 55.5]

    %% --- Load coil1 (P1 excitation) ---
    fprintf('Loading coil1 (P1 excitation) Long2016 verbatim data...\n');
    d = import_ansys_data(res_dir, 'all', 'coil1');

    bbar_mat = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\bs_matrix\Bbar_S_4p572.mat';
    sign_p1 = +1;
    if exist(bbar_mat, 'file')
        Bs = load(bbar_mat, 'col_sign');
        sign_p1 = Bs.col_sign(1);
    end
    fprintf('P1 sign correction: %+d\n', sign_p1);
    if sign_p1 < 0
        d.bx = -d.bx; d.by = -d.by; d.bz = -d.bz;
    end

    F_bx = scatteredInterpolant(d.x, d.y, d.z, d.bx, 'linear', 'none');
    F_by = scatteredInterpolant(d.x, d.y, d.z, d.by, 'linear', 'none');
    F_bz = scatteredInterpolant(d.x, d.y, d.z, d.bz, 'linear', 'none');

    %% --- xz grid (covers coil + block + pole) ---
    xlim_v = [-3, 60];
    zlim_v = [-18, 2];
    grid_x = 130;
    grid_z = 40;
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
    fprintf('  Grid arrows: %d  |B| range: [%.3g, %.3g] T\n', length(Xg), min(bsum_g), max(bsum_g));

    %% --- Color clip + uniform arrow length ---
    b_clip = quantile(bsum_g, 0.92);
    fprintf('  Color clip (92nd pct) = %.3g T\n', b_clip);

    dx_grid = (xlim_v(2)-xlim_v(1)) / grid_x;
    arrow_max = 0.85 * dx_grid;
    b_norm = bsum_g; b_norm(b_norm < eps) = eps;
    scale_per = arrow_max ./ b_norm;
    bx_q = Bx_g .* scale_per;
    bz_q = Bz_g .* scale_per;

    %% --- Plot ---
    fig = figure('Position', [30 30 1700 700], 'Color', 'w');
    set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LBL);

    n_bins  = 24;
    edges_b = linspace(0, b_clip, n_bins + 1);
    cmap_n  = turbo(n_bins);
    lw_rng  = [0.55, 1.90];

    hold on;
    for k = 1:n_bins
        if k < n_bins
            in_bin = bsum_g >= edges_b(k) & bsum_g < edges_b(k+1);
        else
            in_bin = bsum_g >= edges_b(k);
        end
        if any(in_bin)
            lw = lw_rng(1) + (k-1)/(n_bins-1) * (lw_rng(2)-lw_rng(1));
            quiver(Xg(in_bin), Zg(in_bin), bx_q(in_bin), bz_q(in_bin), 0, ...
                   'Color', cmap_n(k,:), 'LineWidth', lw, 'MaxHeadSize', 0.6);
        end
    end

    %% --- Draw POLE outline (revolved cone + cylinder, halfcut z<=-13) ---
    plot([x_tip,  x_base], [Z_CONE,           Z_CONE - POLE_R],  'k-', 'LineWidth', LW_STEEL); % cone slope
    plot([x_base, x_end ], [Z_CONE - POLE_R,  Z_CONE - POLE_R],  'k-', 'LineWidth', LW_STEEL); % cylinder bottom
    % cylinder end face: from cylinder bottom up — but at x=37.5 the block takes over,
    % so cylinder right face is INTERIOR (no air boundary). Skip plotting.
    plot([x_end,  x_tip ], [Z_CONE,           Z_CONE],           'k-', 'LineWidth', LW_STEEL); % halfcut top

    %% --- Draw BLOCK union outline (4 blocks combined at y=0) ---
    %  Stepped polygon from cylinder-end up to post bottom:
    %    (37.5, -13)   cylinder top meets block 3 left at halfcut plane
    %    -> (37.5, -9) block 3 left top
    %    -> (42.5, -9) block 3 top right
    %    -> (42.5, -7) step up to block 1 top
    %  Then block 1 top continues to x=POST_X1 where post takes over
    %  Block bottom edge: from (37.5, -17) to (58.5, -17)
    %  Block right edge:  (58.5, -17) up to (58.5, -7)
    %  Below pole+block junction: (37.5, -16) down to (37.5, -17) — step

    % Block top contour (left to right, above cylinder)
    plot([x_end,   x_end  ], [Z_CONE,          -PROT_H - 2],   'k-', 'LineWidth', LW_STEEL); % (37.5,-13)→(37.5,-9)
    plot([x_end,   POST_X1], [-PROT_H - 2,     -PROT_H - 2],   'k-', 'LineWidth', LW_STEEL); % (37.5,-9)→(42.5,-9)
    plot([POST_X1, POST_X1], [-PROT_H - 2,     POST_Z1],       'k-', 'LineWidth', LW_STEEL); % (42.5,-9)→(42.5,-7)
    % (POST_X1, -7) → (POST_X2, -7) is interior (post sits on top)
    % After post: (POST_X2, -7) → (B1.x2, -7)
    plot([POST_X2, B1.x2  ], [POST_Z1,         POST_Z1],       'k-', 'LineWidth', LW_STEEL); % (52.5,-7)→(58.5,-7)
    plot([B1.x2,   B1.x2  ], [POST_Z1,         B1.z1],         'k-', 'LineWidth', LW_STEEL); % (58.5,-7)→(58.5,-17)
    % Block bottom (continuous from x=37.5 to x=58.5)
    plot([x_end,   B1.x2  ], [B1.z1,           B1.z1],         'k-', 'LineWidth', LW_STEEL); % (37.5,-17)→(58.5,-17)
    % Step at x=37.5 between block bottom (-17) and cylinder bottom (-16)
    plot([x_end,   x_end  ], [B1.z1,           Z_CONE-POLE_R], 'k-', 'LineWidth', LW_STEEL); % (37.5,-17)→(37.5,-16)

    %% --- Draw POST outline (vertical cylinder cross-section) ---
    plot([POST_X1, POST_X1], [POST_Z1, POST_Z2], 'k-', 'LineWidth', LW_STEEL); % left edge
    plot([POST_X1, POST_X2], [POST_Z2, POST_Z2], 'k-', 'LineWidth', LW_STEEL); % top
    plot([POST_X2, POST_X2], [POST_Z2, POST_Z1], 'k-', 'LineWidth', LW_STEEL); % right edge

    %% --- Draw COIL ring (orange filled rectangles, left + right of post) ---
    fill([coilL_x(1) coilL_x(2) coilL_x(2) coilL_x(1)], ...
         [COIL_Z1   COIL_Z1   COIL_Z2   COIL_Z2  ], ...
         COL_COIL, 'EdgeColor', 'k', 'LineWidth', LW_COIL, 'FaceAlpha', 0.55);
    fill([coilR_x(1) coilR_x(2) coilR_x(2) coilR_x(1)], ...
         [COIL_Z1   COIL_Z1   COIL_Z2   COIL_Z2  ], ...
         COL_COIL, 'EdgeColor', 'k', 'LineWidth', LW_COIL, 'FaceAlpha', 0.55);
    text(mean(coilL_x), (COIL_Z1+COIL_Z2)/2, char(8857), ...
         'FontSize', 20, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(mean(coilR_x), (COIL_Z1+COIL_Z2)/2, char(8855), ...
         'FontSize', 20, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

    %% --- Halfcut plane (extends across the pole region only) ---
    plot([xlim_v(1)+1, x_end], [Z_CONE, Z_CONE], '--', 'Color', COL_HALFCUT, ...
         'LineWidth', LW_HALFCUT);
    text(x_end-0.3, Z_CONE-0.6, 'halfcut plane (z=-13, only POLE cut)', ...
         'FontSize', 10, 'Color', COL_HALFCUT, ...
         'HorizontalAlignment', 'right', 'FontAngle', 'italic');

    %% --- WP marker + labels ---
    plot(0, SPH_OFST, 'k+', 'MarkerSize', 16, 'LineWidth', 2.4);
    text(-0.3, SPH_OFST + 0.8, 'WP', 'FontSize', 13, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'right');

    text(x_tip + 0.5, Z_CONE - 1.2, 'P1 tip', 'FontSize', 11, 'FontWeight', 'bold');
    text((x_tip + x_base)/2, Z_CONE - POLE_R - 0.6, 'CONE', ...
         'FontSize', 11, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text((x_base + x_end)/2, Z_CONE - POLE_R - 0.6, 'CYLINDER', ...
         'FontSize', 11, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text((x_end + POST_X1)/2, -PROT_H - 5, 'BLOCK', ...
         'FontSize', 11, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', ...
         'Rotation', 0);
    text((POST_X1 + POST_X2)/2, (POST_Z1 + POST_Z2)/2, 'POST', ...
         'FontSize', 11, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(mean(coilR_x) + 4, (COIL_Z1 + COIL_Z2)/2, 'COIL', ...
         'FontSize', 11, 'HorizontalAlignment', 'left', 'FontWeight', 'bold', ...
         'Color', COL_COIL*0.7);

    % Pole length annotation
    y_anno = -17.3;
    plot([x_tip, x_end], [y_anno, y_anno], 'k-', 'LineWidth', 0.8);
    plot([x_tip, x_tip], [y_anno-0.2, y_anno+0.2], 'k-', 'LineWidth', 0.8);
    plot([x_end, x_end], [y_anno-0.2, y_anno+0.2], 'k-', 'LineWidth', 0.8);
    text((x_tip + x_end)/2, y_anno - 0.4, ...
         sprintf('pole length = %.2f mm', x_end - x_tip), ...
         'FontSize', 10, 'HorizontalAlignment', 'center', 'FontAngle', 'italic');

    hold off;

    %% --- Colorbar + axes ---
    colormap(turbo);
    clim([0, b_clip]);
    cb = colorbar;
    ylabel(cb, sprintf('|B|  [Tesla]  (clipped at 92nd pct = %.2f T)', b_clip), ...
        'FontSize', FONT_CB);

    axis equal;
    grid on;
    xlim(xlim_v); ylim(zlim_v);
    set(gca, 'GridAlpha', 0.18, 'Layer', 'top');
    xlabel('x [mm]  (APDL frame, +x = P1 pole axis from WP)', 'FontSize', FONT_LBL);
    ylabel('z [mm]', 'FontSize', FONT_LBL);
    title(['P1 magnetic structure: COIL + BLOCK + POLE + POST — ' ...
           'B-field side view xz at y=0, halfcut, Long2016 verbatim (VADD fix)'], ...
          'FontSize', FONT_TTL, 'Interpreter', 'none');

    %% --- Save ---
    out_path = fullfile(out_dir, 'P1_Bvector_coil_block_pole.png');
    exportgraphics(fig, out_path, 'Resolution', DPI);
    fprintf('\nSaved: %s\n', out_path);
end
