function plot_P2_geom_vectors_nhat_ahat()
% PLOT_P2_GEOM_VECTORS_NHAT_AHAT  (APDL xz frame, original tilt)
%
% P2 = upper pole at azimuth 180°, cone tilted +36.59° from horizontal.
% Drawn in APDL (lab) frame so the tilt is preserved.
%
% Sensor mounted on TOP slant of cone at s = 4.572 mm from tip.
% n̂ = outward normal of top slant (in APDL xz).
% â = unit vector toward P2 tip.
%
% Sensor value: B·n̂ averaged over Ø0.3 mm Hall sensor disc (matches the
% averaging used by gen_Bbar_S — disc oriented perpendicular to n̂).
%
% Data: Long2016 verbatim baseline, coil5 (P2 paper).

    %% --- Paths ---
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');

    res_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data\coil5\standard';
    out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    %% --- Geometry constants (APDL mm) ---
    cnst = mt_constants();
    R_norm_mm    = cnst.R_norm     * 1e3;   % 0.5
    POLE_R       = cnst.POLE_R     * 1e3;
    POLE_LEN     = cnst.POLE_CONE_LEN * 1e3;
    SPH_OFST_mm  = cnst.SPH_OFST   * 1e3;   % -12.711
    inc          = cnst.upper_incline;       % 36.59° in rad
    theta_p2     = pi;                        % azimuth 180°
    beta_cone    = atan2(POLE_R, POLE_LEN);  % 11.31°
    R_disc_mm    = 0.15;                      % Hall sensor radius (Ø0.3 mm)

    % Pole axis & up_hat in APDL (y=0 slice contains both since theta=180°)
    pa_xz = [cos(inc)*cos(theta_p2); sin(inc)];      % (-0.803, +0.596)
    up_unn = [0; 1] - sin(inc)*pa_xz;
    uh_xz  = up_unn / norm(up_unn);                   % (+0.596, +0.803)
    fprintf('P2 pole_axis (xz) = (%.3f, %.3f)\n', pa_xz);
    fprintf('P2 up_hat   (xz) = (%.3f, %.3f)\n', uh_xz);

    % WP in APDL xz
    WP_xz = [0; SPH_OFST_mm];   % (0, -12.711)

    %% --- Build cone outline in APDL xz (TRIANGLE tilted) ---
    fig = figure('Position', [80 80 1300 750], 'Color', 'w');
    hold on;

    tip_xz       = WP_xz + R_norm_mm * pa_xz;
    base_center  = WP_xz + (R_norm_mm + POLE_LEN) * pa_xz;
    base_top     = base_center + POLE_R * uh_xz;
    base_bottom  = base_center - POLE_R * uh_xz;

    cone_x = [tip_xz(1), base_top(1), base_bottom(1), tip_xz(1)];
    cone_z = [tip_xz(2), base_top(2), base_bottom(2), tip_xz(2)];
    fill(cone_x, cone_z, [0.85 0.85 0.85], 'EdgeColor', 'k', 'LineWidth', 1.6);

    % Holder block — rough rectangle behind cone base (in APDL xz)
    holder_w = 8; holder_h = 2*POLE_R + 4;
    hold_center = base_center + 4 * pa_xz;
    hxs = hold_center(1) + holder_w/2 * [-1 +1 +1 -1] .* pa_xz(1) ...
                          + holder_h/2 * [-1 -1 +1 +1] .* uh_xz(1);
    hzs = hold_center(2) + holder_w/2 * [-1 +1 +1 -1] .* pa_xz(2) ...
                          + holder_h/2 * [-1 -1 +1 +1] .* uh_xz(2);
    fill(hxs, hzs, [0.92 0.92 0.92], 'EdgeColor', 'k', 'LineWidth', 1.0);

    % WP marker
    plot(WP_xz(1), WP_xz(2), 'k+', 'MarkerSize', 14, 'LineWidth', 2);
    text(WP_xz(1) - 0.2, WP_xz(2) + 0.3, 'WP', 'HorizontalAlignment', 'right', ...
         'FontSize', 11, 'FontWeight', 'bold');

    %% --- Sensor: anchor on TOP slant at s=4.572, n̂ = top-slant outward normal ---
    s_anchor = R_norm_mm + 4.572;
    R_at_anchor = (s_anchor - R_norm_mm) / POLE_LEN * POLE_R;
    anchor_xz = WP_xz + s_anchor * pa_xz + R_at_anchor * uh_xz;

    % n̂ in APDL xz = top-slant outward normal
    nhat_xz = -sin(beta_cone) * pa_xz + cos(beta_cone) * uh_xz;
    nhat_xz = nhat_xz / norm(nhat_xz);
    fprintf('n̂ (xz) = (%.3f, %.3f)\n', nhat_xz);

    d_sensor = 0.41;       % canonical Hall sensor placement (perpendicular distance)
    d_max    = 2.0;        % end of sweep (for profile visualization)
    epx = anchor_xz(1) + d_sensor * nhat_xz(1);   % sensor disc center
    epz = anchor_xz(2) + d_sensor * nhat_xz(2);
    sweep_end_x = anchor_xz(1) + d_max * nhat_xz(1);
    sweep_end_z = anchor_xz(2) + d_max * nhat_xz(2);

    % Axis marker
    plot(anchor_xz(1), anchor_xz(2), 'ko', 'MarkerSize', 10, ...
         'MarkerFaceColor', [1 0.7 0], 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
    text(anchor_xz(1) + 0.3, anchor_xz(2) - 0.3, 'axis (r=0)', ...
         'FontSize', 10, 'Color', [0.5 0.35 0]);

    % Green sensor plate at d_max
    face_len = 0.45; face_thk = 0.08;
    v1x = nhat_xz(2);  v1z = -nhat_xz(1);
    sq_x = epx + face_len/2*[+v1x +v1x -v1x -v1x +v1x] ...
               + face_thk/2*[+nhat_xz(1) -nhat_xz(1) -nhat_xz(1) +nhat_xz(1) +nhat_xz(1)];
    sq_z = epz + face_len/2*[+v1z +v1z -v1z -v1z +v1z] ...
               + face_thk/2*[+nhat_xz(2) -nhat_xz(2) -nhat_xz(2) +nhat_xz(2) +nhat_xz(2)];
    fill(sq_x, sq_z, [0.4 0.85 0.4], 'EdgeColor', [0 0.5 0], 'LineWidth', 1.5);
    text(epx + 0.3, epz - 0.05, sprintf('sensor @ %.2f mm', d_sensor), ...
         'FontSize', 10, 'Color', [0 0.5 0], 'FontWeight', 'bold');

    %% --- Load FEM (Long2016 verbatim baseline P2 = coil5) ---
    fprintf('Loading coil5 baseline...\n');
    d = import_ansys_data(res_dir, 'all', 'coil5');

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

    F_bx = scatteredInterpolant(d.x, d.y, d.z, d.bx, 'linear', 'none');
    F_by = scatteredInterpolant(d.x, d.y, d.z, d.by, 'linear', 'none');
    F_bz = scatteredInterpolant(d.x, d.y, d.z, d.bz, 'linear', 'none');

    %% --- Sample 11 points along sweep (centerline) ---
    d_arrows = 0:0.2:2.0;   % mm
    N = length(d_arrows);

    B_xz   = zeros(N, 2);
    a_all  = zeros(N, 1);
    b_all  = zeros(N, 1);
    nh_all = zeros(N, 2);
    ah_all = zeros(N, 2);
    pos_all = zeros(N, 2);

    for j = 1:N
        px = anchor_xz(1) + d_arrows(j)*nhat_xz(1);
        pz = anchor_xz(2) + d_arrows(j)*nhat_xz(2);

        Bx = F_bx(px*1e-3, 0, pz*1e-3);
        Bz = F_bz(px*1e-3, 0, pz*1e-3);
        if isnan(Bx) || isnan(Bz); Bx=0; Bz=0; end

        % â = unit toward P2 tip (in APDL xz)
        dxt = tip_xz(1) - px;  dzt = tip_xz(2) - pz;
        dist_t = sqrt(dxt^2 + dzt^2);
        ahx = dxt/dist_t;  ahz = dzt/dist_t;

        E = [nhat_xz(1) ahx; nhat_xz(2) ahz];
        lambda = E \ [Bx; Bz];

        B_xz(j,:)   = [Bx, Bz];
        a_all(j)    = lambda(1);
        b_all(j)    = lambda(2);
        nh_all(j,:) = nhat_xz';
        ah_all(j,:) = [ahx, ahz];
        pos_all(j,:) = [px, pz];
    end

    %% --- Disc-averaged B·n̂ at sensor end (Ø0.3 mm = R 0.15 mm) ---
    % Disc plane perpendicular to n̂. In 3D: n̂_3d = (nhat_x, 0, nhat_z), and
    % the in-plane axes are v1=(0,1,0) and v2 = n̂ × v1 (in xz plane).
    nhat_3d = [nhat_xz(1); 0; nhat_xz(2)];
    v1_3d   = [0; 1; 0];                                % y-axis
    v2_3d   = cross(nhat_3d, v1_3d);
    v2_3d   = v2_3d / norm(v2_3d);

    N_grid = 41;
    [Ug, Vg] = meshgrid(linspace(-R_disc_mm, R_disc_mm, N_grid), ...
                        linspace(-R_disc_mm, R_disc_mm, N_grid));
    mask = sqrt(Ug.^2 + Vg.^2) <= R_disc_mm;
    Px = epx*1e-3 + Ug*1e-3 * v1_3d(1) + Vg*1e-3 * v2_3d(1);
    Py =  0       + Ug*1e-3 * v1_3d(2) + Vg*1e-3 * v2_3d(2);
    Pz = epz*1e-3 + Ug*1e-3 * v1_3d(3) + Vg*1e-3 * v2_3d(3);

    Bx_disc = F_bx(Px(:), Py(:), Pz(:));
    By_disc = F_by(Px(:), Py(:), Pz(:));
    Bz_disc = F_bz(Px(:), Py(:), Pz(:));

    valid = ~isnan(Bx_disc) & mask(:);
    Bn_disc = Bx_disc(valid)*nhat_3d(1) + By_disc(valid)*nhat_3d(2) + Bz_disc(valid)*nhat_3d(3);
    Bn_avg  = mean(Bn_disc);
    fprintf('Sensor disc avg B·n̂ = %.4f T (%.2f mT), N pts = %d\n', ...
            Bn_avg, Bn_avg*1e3, nnz(valid));
    fprintf('Centerline point B·n̂ = %.4f T\n', a_all(end));

    %% --- Visual scales ---
    Bmag = sqrt(B_xz(:,1).^2 + B_xz(:,2).^2);
    scale_B = 1.5 / max(Bmag) * 0.45;
    scale_c = 1.5 / max(Bmag);

    %% --- Arrows + â guides ---
    h_B = []; h_n = []; h_a = []; h_d = [];
    for j = 1:N
        x0 = pos_all(j,1);  z0 = pos_all(j,2);
        h_dash = plot([x0, tip_xz(1)], [z0, tip_xz(2)], ...
                      '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.6);
        vbx = b_all(j) * ah_all(j,1) * scale_c;
        vbz = b_all(j) * ah_all(j,2) * scale_c;
        h_beta = quiver(x0, z0, vbx, vbz, 0, 'Color', [0.9 0.2 0.2], ...
                        'LineWidth', 1.6, 'MaxHeadSize', 0.6, 'AutoScale', 'off');
        vax = a_all(j) * nh_all(j,1) * scale_c;
        vaz = a_all(j) * nh_all(j,2) * scale_c;
        h_alpha = quiver(x0, z0, vax, vaz, 0, 'Color', [0 0.45 0.95], ...
                         'LineWidth', 1.4, 'MaxHeadSize', 1.5, 'AutoScale', 'off');
        vBx = B_xz(j,1) * scale_B;  vBz = B_xz(j,2) * scale_B;
        h_full = quiver(x0, z0, vBx, vBz, 0, 'Color', 'k', 'LineWidth', 1.8, ...
                        'MaxHeadSize', 0.9, 'AutoScale', 'off');
        if j == 1; h_B=h_full; h_n=h_alpha; h_a=h_beta; h_d=h_dash; end
    end

    %% --- Sensor-end disc-averaged value box (upper-right, inside zoom window) ---
    text(-2.95, -10.25, ...
         sprintf('$\\langle\\mathbf{B}\\!\\cdot\\!\\hat{\\mathbf{n}}\\rangle_{\\mathrm{disc}} = %.4f$ T', Bn_avg), ...
         'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold', ...
         'Color', [0 0.35 0.85], 'BackgroundColor', [1 1 1 0.85], ...
         'EdgeColor', [0 0.35 0.85], 'Margin', 3, ...
         'HorizontalAlignment', 'left');

    %% --- Legend + axes ---
    legend([h_B, h_n, h_a, h_d], ...
        {sprintf('$\\mathbf{B}$ (FEM, scale: 1\\,mm = %.3f\\,T)', 1/scale_B), ...
         sprintf('$\\alpha\\,\\hat{\\mathbf{n}}$ (scale: 1\\,mm = %.3f\\,T)', 1/scale_c), ...
         sprintf('$\\beta\\,\\hat{\\mathbf{a}}$ (scale: 1\\,mm = %.3f\\,T)', 1/scale_c), ...
         '$\hat{\mathbf{a}}$ direction (toward P2 tip)'}, ...
        'Interpreter', 'latex', 'FontSize', 11, 'Location', 'northwest');

    axis equal; grid on; box on;
    xlim([-5, -2]); ylim([-10.5, -8]);
    set(gca, 'FontSize', 11);
    xlabel('x [mm]', 'FontSize', 12);
    ylabel('z [mm]', 'FontSize', 12);
    title('P2 baseline (Long2016 verbatim) -- APDL frame, upper pole tilted $+36.59^\circ$', ...
          'FontSize', 13, 'Interpreter', 'latex');

    %% --- Save ---
    out_path = fullfile(out_dir, 'P2_geom_vectors_nhat_ahat.png');
    exportgraphics(fig, out_path, 'Resolution', 300);
    fprintf('Saved: %s\n', out_path);
end
