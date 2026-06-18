function plot_P1_geom_vectors_neg()
% PLOT_P1_GEOM_VECTORS_NEG
%
% Variant of PLOT_P1_GEOM_VECTORS_NHAT_AHAT with:
%   - FEM B data multiplied by -1 (all three arrows reversed)
%   - gray dashed a-hat guide lines removed (and dropped from legend)
%   - everything else (geometry, annotations, frame) unchanged
%   - output to a SEPARATE file (original figure/script left intact)
%
% Side-view xz plane (y=0) of P1 halfcut cone, with magnetic flux vectors
% sampled along a sweep path (anchor → air, 2 mm) and decomposed into
% oblique pair B = α·n̂ + β·â, where:
%   n̂ = sweep direction (outward normal at sensor face, "blue")
%   â = direction toward P1 tip (axis-pointing,        "red")
%
% Data: Long2016 verbatim baseline P1 (magnetic_sim/ANSYS/main/ANSYS_data/long2016_hexapole_halfcut/coil1/)
%       I_in = 1.0 A, halfcut geometry.

    %% --- Paths ---
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');

    res_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data\coil1\standard';
    out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    %% --- Geometry (mm, halfcut: lower D-shape only, z ≤ Z_CONE) ---
    R_norm_xy     = 0.500*sqrt(2/3);     % 0.408
    POLE_R        = 3.0;
    POLE_CONE_LEN = 15.0;
    Z_CONE        = -13.0;
    X0_LOW        = 47.5;
    X_SHOULDER    = R_norm_xy + POLE_CONE_LEN;

    fig = figure('Position', [80 80 1200 650], 'Color', 'w');
    hold on;

    % HALFCUT D-shape cone (lower half only)
    cone_x = [R_norm_xy, X_SHOULDER, X_SHOULDER, R_norm_xy];
    cone_z = [Z_CONE,    Z_CONE-POLE_R, Z_CONE, Z_CONE];
    fill(cone_x, cone_z, [0.85 0.85 0.85], 'EdgeColor', 'k', 'LineWidth', 1.6);

    % Holder BLOCKs (rough side view)
    rectangle('Position', [X0_LOW-5, Z_CONE-10, 16, 10], ...
        'FaceColor', [0.92 0.92 0.92], 'EdgeColor', 'k', 'LineWidth', 1.0);
    rectangle('Position', [X0_LOW-10, Z_CONE-10, 8, 8], ...
        'FaceColor', [0.95 0.95 0.95], 'EdgeColor', 'k', 'LineWidth', 1.0);

    % WP marker
    MT_CTR_ZOFST = -7.0 - 6.0 + 0.500/sqrt(3);   % -12.711
    plot(0, MT_CTR_ZOFST, 'k+', 'MarkerSize', 14, 'LineWidth', 2);
    text(-0.05, MT_CTR_ZOFST + 0.25, 'WP', 'HorizontalAlignment', 'right', ...
         'FontSize', 11, 'FontWeight', 'bold');

    %% --- Sensor: 4.572 mm along axis + 0.41 mm perpendicular to milled flat ---
    % P1 lower halfcut: milled flat is HORIZONTAL (z = Z_CONE = -13 mm)
    % → n̂ = perpendicular to milled flat = pure +z (no tilt)
    nhx = 0;  nhz = 1;
    anchor   = [R_norm_xy + 4.572, Z_CONE];   % on milled flat
    d_sensor = 0.41;                            % mm, canonical Hall sensor placement
    d_max    = 2.0;                             % mm, end of sweep (for profile visualization)
    epx      = anchor(1) + d_sensor*nhx;        % sensor position
    epz      = anchor(2) + d_sensor*nhz;
    sweep_end_x = anchor(1) + d_max*nhx;
    sweep_end_z = anchor(2) + d_max*nhz;
    tip_x_p1 = R_norm_xy;
    tip_z_p1 = Z_CONE;

    % Axis marker at d=0
    plot(anchor(1), anchor(2), 'ko', 'MarkerSize', 10, ...
         'MarkerFaceColor', [1 0.7 0], 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
    % [MODIFIED] 'axis (r=0)' text label removed

    % Green sensor plate at canonical 0.41 mm (parallel to milled flat)
    face_len = 0.45;  face_thk = 0.08;
    v1x = nhz;  v1z = -nhx;
    sq_x = epx + face_len/2 * [+v1x, +v1x, -v1x, -v1x, +v1x] ...
               + face_thk/2 * [+nhx, -nhx, -nhx, +nhx, +nhx];
    sq_z = epz + face_len/2 * [+v1z, +v1z, -v1z, -v1z, +v1z] ...
               + face_thk/2 * [+nhz, -nhz, -nhz, +nhz, +nhz];
    fill(sq_x, sq_z, [0.4 0.85 0.4], 'EdgeColor', [0 0.5 0], 'LineWidth', 1.5);
    text(epx + 0.4, epz + 0.05, sprintf('sensor @ %.2f mm', d_sensor), ...
         'FontSize', 10, 'Color', [0 0.5 0], 'FontWeight', 'bold');

    %% --- Load FEM (Long2016 verbatim P1 baseline) ---
    fprintf('Loading coil1 baseline (Long2016 verbatim, halfcut)...\n');
    d = import_ansys_data(res_dir, 'all', 'coil1');
    F_bx = scatteredInterpolant(d.x, d.y, d.z, d.bx, 'linear', 'none');
    F_by = scatteredInterpolant(d.x, d.y, d.z, d.by, 'linear', 'none');
    F_bz = scatteredInterpolant(d.x, d.y, d.z, d.bz, 'linear', 'none');

    %% --- Sample 11 points along sweep (every 0.2 mm) ---
    d_arrows = 0:0.2:2.0;   % mm
    N = length(d_arrows);

    B_all   = zeros(N, 2);   % [Bx Bz] in T at sample point
    a_all   = zeros(N, 1);   % α coefficient (T)  — along n̂
    b_all   = zeros(N, 1);   % β coefficient (T)  — along â
    nh_all  = zeros(N, 2);   % n̂ vector (constant per sample)
    ah_all  = zeros(N, 2);   % â vector (toward P1 tip, varies per sample)
    pos_all = zeros(N, 2);
    for j = 1:N
        pos_x = anchor(1) + d_arrows(j) * nhx;
        pos_z = anchor(2) + d_arrows(j) * nhz;

        % â = unit vector from sample → P1 tip
        dx_t  = tip_x_p1 - pos_x;
        dz_t  = tip_z_p1 - pos_z;
        dist_t = sqrt(dx_t^2 + dz_t^2);
        ahx = dx_t/dist_t;  ahz = dz_t/dist_t;

        % FEM B at sample (in APDL meters, y=0 slice)
        Bx = -F_bx(pos_x*1e-3, 0, pos_z*1e-3);   % [MODIFIED] negate data
        Bz = -F_bz(pos_x*1e-3, 0, pos_z*1e-3);   % [MODIFIED] negate data
        if isnan(Bx) || isnan(Bz)
            Bx = 0; Bz = 0;
        end

        % Oblique decomposition B = α·n̂ + β·â  (E·λ = B)
        E = [nhx ahx; nhz ahz];
        lambda = E \ [Bx; Bz];

        B_all(j,:)   = [Bx, Bz];
        a_all(j)     = lambda(1);     % along n̂
        b_all(j)     = lambda(2);     % along â
        nh_all(j,:)  = [nhx, nhz];
        ah_all(j,:)  = [ahx, ahz];
        pos_all(j,:) = [pos_x, pos_z];
    end

    %% --- Scale visualisation ---
    Bmag = sqrt(B_all(:,1).^2 + B_all(:,2).^2);
    scale_B = 1.5 / max(Bmag) * 0.45;     % mm per Tesla, shorter for B
    scale_c = 1.5 / max(Bmag);            % mm per Tesla, full for components
    fprintf('Vector scale: 1 mm = %.0f mT (component), 1 mm = %.0f mT (B)\n', ...
            1/scale_c*1e3, 1/scale_B*1e3);

    %% --- Draw arrows (dashed â guide lines removed) ---
    h_B = []; h_n = []; h_a = [];
    for j = 1:N
        x0 = pos_all(j,1);  z0 = pos_all(j,2);

        % [MODIFIED] dashed gray â guide line removed

        % β·â vector (RED, toward tip)
        vbx = b_all(j) * ah_all(j,1) * scale_c;
        vbz = b_all(j) * ah_all(j,2) * scale_c;
        h_beta = quiver(x0, z0, vbx, vbz, 0, 'Color', [0.9 0.2 0.2], ...
                        'LineWidth', 1.6, 'MaxHeadSize', 0.6, 'AutoScale', 'off');

        % α·n̂ vector (BLUE, sweep direction)
        vax = a_all(j) * nh_all(j,1) * scale_c;
        vaz = a_all(j) * nh_all(j,2) * scale_c;
        h_alpha = quiver(x0, z0, vax, vaz, 0, 'Color', [0 0.45 0.95], ...
                         'LineWidth', 1.4, 'MaxHeadSize', 1.5, 'AutoScale', 'off');

        % B vector (BLACK, shorter, on top)
        vBx = B_all(j,1) * scale_B;  vBz = B_all(j,2) * scale_B;
        h_full = quiver(x0, z0, vBx, vBz, 0, 'Color', 'k', 'LineWidth', 1.8, ...
                        'MaxHeadSize', 0.9, 'AutoScale', 'off');

        if j == 1
            h_B = h_full;  h_n = h_alpha;  h_a = h_beta;
        end
    end

    %% --- Disc-averaged B·n̂ at canonical sensor position (Ø0.3 mm) ---
    R_disc_mm = 0.15;
    nhat_3d = [nhx; 0; nhz];
    v1_3d = [0; 1; 0];
    v2_3d = cross(nhat_3d, v1_3d);  v2_3d = v2_3d/norm(v2_3d);
    N_grid = 41;
    [Ug, Vg] = meshgrid(linspace(-R_disc_mm, R_disc_mm, N_grid), ...
                        linspace(-R_disc_mm, R_disc_mm, N_grid));
    mask = sqrt(Ug.^2 + Vg.^2) <= R_disc_mm;
    Px = epx*1e-3 + Ug*1e-3*v1_3d(1) + Vg*1e-3*v2_3d(1);
    Py =       0  + Ug*1e-3*v1_3d(2) + Vg*1e-3*v2_3d(2);
    Pz = epz*1e-3 + Ug*1e-3*v1_3d(3) + Vg*1e-3*v2_3d(3);
    Bx_d = F_bx(Px(:), Py(:), Pz(:));
    By_d = F_by(Px(:), Py(:), Pz(:));
    Bz_d = F_bz(Px(:), Py(:), Pz(:));
    valid = ~isnan(Bx_d) & mask(:);
    Bn_disc = Bx_d(valid)*nhat_3d(1) + By_d(valid)*nhat_3d(2) + Bz_d(valid)*nhat_3d(3);
    Bn_avg = -mean(Bn_disc);   % [MODIFIED] negate data
    fprintf('P1 disc-avg B·n̂ @ %.2f mm = %.4f T (%.2f mT), N=%d\n', ...
            d_sensor, Bn_avg, Bn_avg*1e3, nnz(valid));
    % [MODIFIED] disc-avg ⟨B·n̂⟩ text box removed

    %% --- Legend (dashed â direction entry removed) ---
    legend([h_B, h_n, h_a], ...
        {'$\mathbf{B}$ (FEM)', ...               % [MODIFIED] scale text removed
         '$\alpha\,\hat{\mathbf{n}}$', ...
         '$\beta\,\hat{\mathbf{a}}$'}, ...
        'Interpreter', 'latex', 'FontSize', 11, 'Location', 'northeast');

    %% --- Axes ---
    axis equal; grid on; box on;
    xlim([3.5, 6.5]); ylim([-13.5, -11.5]);
    set(gca, 'FontSize', 11);
    xlabel('x [mm]', 'FontSize', 12);
    ylabel('z [mm]', 'FontSize', 12);
    % [MODIFIED] title removed

    %% --- Save ---
    out_path = fullfile(out_dir, 'P1_geom_vectors_neg.png');   % [MODIFIED] separate output file
    exportgraphics(fig, out_path, 'Resolution', 300);
    fprintf('Saved: %s\n', out_path);
end
