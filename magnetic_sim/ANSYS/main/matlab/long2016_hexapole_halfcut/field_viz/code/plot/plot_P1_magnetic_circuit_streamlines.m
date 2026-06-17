function plot_P1_magnetic_circuit_streamlines()
% PLOT_P1_MAGNETIC_CIRCUIT_STREAMLINES
%   Magnetic circuit diagram of P1 (lower pole) drawn as streamlines (flux
%   lines) traced from FEM B field — raw data, no sign flip, no Bbar_S sign
%   correction (col_sign auto-detect untouched).
%
%   Geometry overlaid: COIL + BLOCK + POLE (cone+cylinder, halfcut) + POST.
%   Streamlines colored by local |B|.  Side view xz at y=0, APDL frame.
%
%   Closed flux loop (CCW current from above, moment +z) should show:
%     - In post:    flux UP (+z)
%     - In yoke:    flux radially out / azimuthal (to other poles)
%     - In block:   flux from below into post
%     - In cone:    flux from cylinder → tip → out into WP air
%     - Outer air:  returning loop from yoke top around back to block bottom

    %% --- Style ---
    FONT_LBL  = 12;
    FONT_TTL  = 13;
    FONT_CB   = 11;
    LW_STEEL  = 2.0;
    LW_STREAM = 1.4;
    DPI       = 300;

    COL_COIL    = [0.93 0.55 0.10];
    COL_HALFCUT = [0.45 0.45 0.45];

    %% --- Paths ---
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');

    cnst = mt_constants();

    res_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\coil1\standard';
    out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    %% --- Geometry (APDL mm) ---
    R_norm_xy = cnst.R_norm_xy * 1e3;
    R_norm_z  = cnst.R_norm_z  * 1e3;
    SPH_OFST  = cnst.SPH_OFST  * 1e3;
    POLE_R    = cnst.POLE_R    * 1e3;
    POLE_LEN  = cnst.POLE_CONE_LEN * 1e3;
    Z_CONE    = SPH_OFST - R_norm_z;          % -13

    x_tip  = R_norm_xy;
    x_base = x_tip + POLE_LEN;
    X0_LOW = (42 + 53)/2;                     % 47.5
    x_end  = X0_LOW - 10;                     % 37.5

    PROT_R  = 5;
    PROT_H  = 7;
    POST_X1 = X0_LOW - PROT_R;
    POST_X2 = X0_LOW + PROT_R;
    POST_Z1 = -PROT_H;
    POST_Z2 = 0;

    YOKE_IN_R  = 42;
    YOKE_OUT_R = 53;
    YOKE_H     = 2;

    % 4 connector blocks (union at y=0): stepped shape
    %   x=[37.5, 42.5]: z=[-17, -9]
    %   x=[42.5, 58.5]: z=[-17, -7]
    B_x1 = x_end;          % 37.5
    B_x2 = X0_LOW + 11;    % 58.5
    B_z_lo = -PROT_H - 10; % -17

    COIL_IN_R  = 5;
    COIL_OUT_R = 8;

    %% --- Load coil1 (P1 excitation) — RAW FEM, no sign flip ---
    fprintf('Loading coil1 (P1 excitation) smrt 4 raw FEM data...\n');
    d = import_ansys_data(res_dir, 'all', 'coil1');
    fprintf('  Nodes loaded: %d\n', length(d.x));

    F_bx = scatteredInterpolant(d.x, d.y, d.z, d.bx, 'linear', 'none');
    F_bz = scatteredInterpolant(d.x, d.y, d.z, d.bz, 'linear', 'none');

    %% --- Uniform 2D grid in xz at y=0 for streamline computation ---
    xlim_v = [-3, 60];
    zlim_v = [-18, 4];
    nx = 260;
    nz = 90;
    x_g = linspace(xlim_v(1), xlim_v(2), nx);
    z_g = linspace(zlim_v(1), zlim_v(2), nz);
    [Xg, Zg] = meshgrid(x_g, z_g);

    fprintf('Interpolating B field onto %dx%d grid...\n', nx, nz);
    BX = F_bx(Xg*1e-3, zeros(size(Xg)), Zg*1e-3);
    BZ = F_bz(Xg*1e-3, zeros(size(Xg)), Zg*1e-3);
    BX(isnan(BX)) = 0;
    BZ(isnan(BZ)) = 0;
    Bmag = sqrt(BX.^2 + BZ.^2);

    %% --- Seed points for streamline tracing ---
    %  Strategy: line-seeds across key cross-sections so the integrator captures
    %  the dominant flux paths but doesn't drown the figure in lines.
    seed_x = []; seed_z = [];

    % Inside POST vertical line (capture flux UP through post)
    seed_x = [seed_x, ones(1,7)*X0_LOW];
    seed_z = [seed_z, linspace(POST_Z1+0.5, POST_Z2-0.5, 7)];

    % Inside CYLINDER horizontal line (capture flux toward tip)
    seed_x = [seed_x, linspace(x_base+1, x_end-1, 8)];
    seed_z = [seed_z, ones(1,8)*(Z_CONE - POLE_R/2)];

    % Inside CONE diagonal line (capture flux through cone)
    seed_x = [seed_x, linspace(x_tip+0.5, x_base-0.5, 6)];
    seed_z = [seed_z, Z_CONE - (linspace(x_tip+0.5, x_base-0.5, 6) - x_tip)/(x_base - x_tip)*POLE_R*0.6];

    % Inside BLOCK lower region (capture flux entering from below)
    seed_x = [seed_x, linspace(B_x1+1, B_x2-1, 7)];
    seed_z = [seed_z, ones(1,7)*(B_z_lo + 1.5)];

    % Inside YOKE (capture flux leaving post top)
    seed_x = [seed_x, linspace(POST_X1+0.5, YOKE_OUT_R-0.5, 5)];
    seed_z = [seed_z, ones(1,5)*(YOKE_H/2)];

    % Air around tip + WP (capture flux exiting tip)
    seed_x = [seed_x, -2.5, -1.5, -0.5, 0.5];
    seed_z = [seed_z, Z_CONE,  Z_CONE,  Z_CONE+0.3,  Z_CONE+0.6];

    % Air above the pole (capture upward leakage from halfcut plane)
    seed_x = [seed_x, linspace(2, 30, 8)];
    seed_z = [seed_z, ones(1,8)*(Z_CONE + 2)];

    fprintf('Total seed points: %d\n', length(seed_x));

    %% --- Compute streamlines (BOTH forward and reverse for closed loops) ---
    fprintf('Tracing streamlines...\n');
    opts = [0.1, 2000];   % step size 0.1, max 2000 steps

    XY_fwd = stream2(Xg, Zg, BX, BZ, seed_x, seed_z, opts);
    XY_rev = stream2(Xg, Zg, -BX, -BZ, seed_x, seed_z, opts);

    %% --- Plot ---
    fig = figure('Position', [30 30 1700 700], 'Color', 'w');
    set(fig, 'DefaultAxesFontSize', FONT_LBL);
    hold on;

    % Background |B| as filled contour (subtle)
    Bclip = quantile(Bmag(:), 0.96);
    contourf(Xg, Zg, min(Bmag, Bclip), 24, 'LineColor', 'none');
    colormap(turbo);
    clim([0, Bclip]);
    cb = colorbar;
    ylabel(cb, sprintf('|B|  [Tesla]  (background, clipped at 96th pct = %.2f T)', Bclip), ...
        'FontSize', FONT_CB);
    set(cb, 'FontSize', FONT_CB);

    % Streamlines on top
    sl_fwd = streamline(XY_fwd);
    sl_rev = streamline(XY_rev);
    set([sl_fwd; sl_rev], 'Color', [0.1 0.1 0.1], 'LineWidth', LW_STREAM);

    % Add arrowheads at midpoint of each streamline
    add_arrows_to_streams(XY_fwd, [0 0 0]);

    %% --- Draw STEEL outlines on top of streamlines ---
    %  POLE (cone + cylinder, halfcut lower half)
    plot([x_tip,  x_base], [Z_CONE,           Z_CONE - POLE_R],  'k-', 'LineWidth', LW_STEEL);
    plot([x_base, x_end ], [Z_CONE - POLE_R,  Z_CONE - POLE_R],  'k-', 'LineWidth', LW_STEEL);
    plot([x_end,  x_tip ], [Z_CONE,           Z_CONE],           'k-', 'LineWidth', LW_STEEL);

    %  BLOCK union stepped polygon
    plot([x_end,   x_end  ], [Z_CONE,          -PROT_H - 2],   'k-', 'LineWidth', LW_STEEL);
    plot([x_end,   POST_X1], [-PROT_H - 2,     -PROT_H - 2],   'k-', 'LineWidth', LW_STEEL);
    plot([POST_X1, POST_X1], [-PROT_H - 2,     POST_Z1],       'k-', 'LineWidth', LW_STEEL);
    plot([POST_X2, B_x2   ], [POST_Z1,         POST_Z1],       'k-', 'LineWidth', LW_STEEL);
    plot([B_x2,    B_x2   ], [POST_Z1,         B_z_lo],        'k-', 'LineWidth', LW_STEEL);
    plot([x_end,   B_x2   ], [B_z_lo,          B_z_lo],        'k-', 'LineWidth', LW_STEEL);
    plot([x_end,   x_end  ], [B_z_lo,          Z_CONE-POLE_R], 'k-', 'LineWidth', LW_STEEL);

    %  POST
    plot([POST_X1, POST_X1], [POST_Z1, POST_Z2], 'k-', 'LineWidth', LW_STEEL);
    plot([POST_X1, POST_X2], [POST_Z2, POST_Z2], 'k-', 'LineWidth', LW_STEEL);
    plot([POST_X2, POST_X2], [POST_Z2, POST_Z1], 'k-', 'LineWidth', LW_STEEL);

    %  YOKE (right portion in view)
    plot([POST_X2-2, YOKE_OUT_R, YOKE_OUT_R, POST_X2-2, POST_X2-2], ...
         [0,         0,          YOKE_H,     YOKE_H,    0], 'k-', 'LineWidth', LW_STEEL);

    %% --- COIL ---
    coilL_x = [X0_LOW-COIL_OUT_R, X0_LOW-COIL_IN_R, X0_LOW-COIL_IN_R, X0_LOW-COIL_OUT_R];
    coilR_x = [X0_LOW+COIL_IN_R,  X0_LOW+COIL_OUT_R, X0_LOW+COIL_OUT_R, X0_LOW+COIL_IN_R];
    coil_z  = [POST_Z1, POST_Z1, POST_Z2, POST_Z2];
    fill(coilL_x, coil_z, COL_COIL, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.65);
    fill(coilR_x, coil_z, COL_COIL, 'EdgeColor', 'k', 'LineWidth', 1.2, 'FaceAlpha', 0.65);
    text(mean([X0_LOW-COIL_OUT_R, X0_LOW-COIL_IN_R]), (POST_Z1+POST_Z2)/2, char(8857), ...
         'FontSize', 20, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(mean([X0_LOW+COIL_IN_R, X0_LOW+COIL_OUT_R]), (POST_Z1+POST_Z2)/2, char(8855), ...
         'FontSize', 20, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

    %% --- Halfcut plane (annotation, dashed) ---
    plot([xlim_v(1)+1, x_end], [Z_CONE Z_CONE], '--', 'Color', COL_HALFCUT, 'LineWidth', 1.0);

    %% --- WP marker + labels ---
    plot(0, SPH_OFST, 'w+', 'MarkerSize', 18, 'LineWidth', 3);
    plot(0, SPH_OFST, 'k+', 'MarkerSize', 16, 'LineWidth', 2);
    text(-0.4, SPH_OFST + 1.2, 'WP', 'FontSize', 13, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'right', 'BackgroundColor', [1 1 1 0.6]);

    text(x_tip + 0.6, Z_CONE - 1.2, 'P1 tip', 'FontSize', 11, 'FontWeight', 'bold', ...
         'BackgroundColor', [1 1 1 0.6]);
    text((x_tip + x_base)/2, Z_CONE - POLE_R - 0.8, 'CONE', 'FontSize', 11, ...
         'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 1 0.6]);
    text((x_base + x_end)/2, Z_CONE - POLE_R - 0.8, 'CYLINDER', 'FontSize', 11, ...
         'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 1 0.6]);
    text((x_end + POST_X1)/2, -PROT_H - 5.5, 'BLOCK', 'FontSize', 11, ...
         'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 1 0.6]);
    text((POST_X1 + POST_X2)/2, -3.5, 'POST', 'FontSize', 11, ...
         'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 1 0.6]);
    text((POST_X2+YOKE_OUT_R)/2-1, YOKE_H/2, 'YOKE', 'FontSize', 11, ...
         'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 1 0.6]);
    text(X0_LOW, POST_Z1-1.0, 'COIL', 'FontSize', 11, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'Color', COL_COIL*0.6, 'BackgroundColor', [1 1 1 0.8]);

    hold off;

    %% --- Axes ---
    axis equal;
    xlim(xlim_v); ylim(zlim_v);
    box on;
    set(gca, 'Layer', 'top');
    xlabel('x [mm]  (APDL frame, +x = P1 pole axis from WP)', 'FontSize', FONT_LBL);
    ylabel('z [mm]', 'FontSize', FONT_LBL);
    title(['P1 magnetic circuit (streamlines from raw FEM B field) — ' ...
           'side view xz at y=0, halfcut, smrt 4 mesh, CCW coil (moment +z)'], ...
          'FontSize', FONT_TTL, 'Interpreter', 'none');

    %% --- Save ---
    out_path = fullfile(out_dir, 'P1_magnetic_circuit_streamlines.png');
    exportgraphics(fig, out_path, 'Resolution', DPI);
    fprintf('\nSaved: %s\n', out_path);
end


function add_arrows_to_streams(XY_list, color)
% Add a single arrowhead at the midpoint of each streamline.
    for k = 1:length(XY_list)
        pts = XY_list{k};
        if isempty(pts) || size(pts, 1) < 3
            continue
        end
        mid = round(size(pts,1)/2);
        if mid <= 1, continue; end
        p0 = pts(mid-1, :);
        p1 = pts(mid+1, :);
        dx = p1(1) - p0(1);
        dz = p1(2) - p0(2);
        nrm = sqrt(dx^2 + dz^2);
        if nrm < 1e-6, continue; end
        % Short arrow at midpoint
        dx = dx/nrm * 0.6;
        dz = dz/nrm * 0.6;
        quiver(pts(mid,1)-dx/2, pts(mid,2)-dz/2, dx, dz, 0, ...
               'Color', color, 'LineWidth', 1.1, 'MaxHeadSize', 4, 'AutoScale','off');
    end
end
