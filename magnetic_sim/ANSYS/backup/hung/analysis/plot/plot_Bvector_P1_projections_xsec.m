%% plot_Bvector_P1_projections_xsec.m
%  Cross-section (axis-perpendicular) view of B decomposition inside P1.
%  Cutting plane is perpendicular to P1 axis at station s; viewer looks along -a_hat
%  (from tip toward pole body), seeing pole as a circle of radius R_surf.
%
%    r_hat : outward surface normal
%      cylinder -> pure in-plane radial
%      cone     -> cos(θ_c)·radial − sin(θ_c)·a_hat (has small out-of-plane part)
%    n_hat : toward P1 tip = -a_hat  (fully out-of-plane in this view)
%
%  Visualization:
%    α component -> in-plane arrow (radial), length + colour ∝ |α|
%    β component -> filled circle marker (⊙ = out of page toward tip), size + colour ∝ |β|
%
%  Outputs:
%    magnetic_sim/hung/figures/coil1/Bvector_P1_projections_RoundFillet_l500_xsec_cyl.png   (s=25)
%    magnetic_sim/hung/figures/coil1/Bvector_P1_projections_RoundFillet_l500_xsec_cone.png  (s=10)

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));

%% 1. Load full 3D FEM data + build 3D interpolants
fprintf('Loading RoundFillet l=500 full-model data...\n');
d = import_ansys_data( ...
    fullfile('..','..','results','coil1','round_filleted_conv'), ...
    'all', 'coil1');
N_total = length(d.bx);
fprintf('  total nodes: %d\n', N_total);

fprintf('Building 3D scatteredInterpolants (this takes ~30 s)...\n');
tic;
F_Bx_3d = scatteredInterpolant(d.x, d.y, d.z, d.bx, 'linear', 'nearest');
F_By_3d = scatteredInterpolant(d.x, d.y, d.z, d.by, 'linear', 'nearest');
F_Bz_3d = scatteredInterpolant(d.x, d.y, d.z, d.bz, 'linear', 'nearest');
fprintf('  built in %.1f s\n', toc);

%% 2. Pole geometry (3D)
tip_3d = [0.408e-3; 0; -0.289e-3];              % m
tilt  = -5.71;
a_hat_3d  = [cosd(tilt); 0; sind(tilt)];        % axial tip->end
n_perp_3d = [-sind(tilt); 0; cosd(tilt)];       % perp in xz plane
u_hat_3d  = [0; 1; 0];                          % ANSYS +y
v_hat_3d  = n_perp_3d;                          % in-plane basis (v = n_perp)

POLE_R        = 3.175e-3;                       % m
POLE_CONE_LEN = 15.875e-3;
CONE_HALF_ANG = 11.31;                          % deg

%% 3. Plot configurations
plots = struct();
plots(1).name      = 'cylinder';
plots(1).section   = 'P1 cylinder';
plots(1).s_anchor  = 25e-3;                     % m
plots(1).r_surf    = POLE_R;
plots(1).is_cone   = false;
plots(1).xylim     = 6.0;                       % mm half-range
plots(1).filename  = 'Bvector_P1_projections_RoundFillet_l500_xsec_cyl.png';

plots(2).name      = 'cone';
plots(2).section   = 'P1 cone';
plots(2).s_anchor  = 10e-3;
plots(2).r_surf    = 10e-3 * tand(CONE_HALF_ANG);   % ≈ 2.0 mm
plots(2).is_cone   = true;
plots(2).xylim     = 6.0;
plots(2).filename  = 'Bvector_P1_projections_RoundFillet_l500_xsec_cone.png';

%% 4. Loop over sections
for k = 1:length(plots)
    p = plots(k);
    fprintf('\n=== %s (axis anchor s=%.1f mm) ===\n', p.section, p.s_anchor*1e3);

    % Cross-section origin (3D, meters)
    center_3d = tip_3d + p.s_anchor * a_hat_3d;

    % --- Sample points: 12 angles (every 30°) × 2 radii (mid + surface) ---
    rs = p.r_surf;
    angles_deg = 0:30:330;                      % 12 angles
    radii_frac = [0.5, 1.0];                    % mid, surface
    samples_uv = zeros(0, 2);
    for th = angles_deg
        for rf = radii_frac
            r = rf * rs;
            samples_uv(end+1, :) = [r*cosd(th), r*sind(th)]; %#ok<AGROW>
        end
    end
    N_samp = size(samples_uv, 1);
    fprintf('  %d angular x %d radial = %d samples\n', ...
        length(angles_deg), length(radii_frac), N_samp);

    % --- Convert to 3D positions then query 3D interpolants ---
    pos_3d = zeros(N_samp, 3);
    for i = 1:N_samp
        pos_3d(i, :) = (center_3d + samples_uv(i,1)*u_hat_3d + samples_uv(i,2)*v_hat_3d)';
    end
    Bx = F_Bx_3d(pos_3d(:,1), pos_3d(:,2), pos_3d(:,3));
    By = F_By_3d(pos_3d(:,1), pos_3d(:,2), pos_3d(:,3));
    Bz = F_Bz_3d(pos_3d(:,1), pos_3d(:,2), pos_3d(:,3));

    % --- Compute α via 3D oblique decomposition (same math as side view) ---
    %   B_3d  ≈  α·r̂_3d + β·n̂_3d   (non-orthogonal basis, n̂ per-sample to tip)
    %   [1 c; c 1] [α;β] = [B·r̂; B·n̂]    where c = r̂·n̂
    r_local = sqrt(samples_uv(:,1).^2 + samples_uv(:,2).^2);
    alpha_T = zeros(N_samp, 1);
    rhat_inplane_uv = zeros(N_samp, 2);
    for i = 1:N_samp
        B_3d = [Bx(i); By(i); Bz(i)];
        sample_3d = center_3d + samples_uv(i,1)*u_hat_3d + samples_uv(i,2)*v_hat_3d;

        rhat_inplane_uv(i,:) = samples_uv(i,:) / r_local(i);
        rhat_3d_inplane = rhat_inplane_uv(i,1)*u_hat_3d + rhat_inplane_uv(i,2)*v_hat_3d;
        if p.is_cone
            r_hat_3d = cosd(CONE_HALF_ANG)*rhat_3d_inplane - sind(CONE_HALF_ANG)*a_hat_3d;
        else
            r_hat_3d = rhat_3d_inplane;
        end

        d_to_tip = tip_3d - sample_3d;
        n_hat_3d = d_to_tip / norm(d_to_tip);

        c       = dot(r_hat_3d, n_hat_3d);
        B_dot_r = dot(B_3d, r_hat_3d);
        B_dot_n = dot(B_3d, n_hat_3d);
        alpha_T(i) = (B_dot_r - c*B_dot_n) / (1 - c^2);
    end
    alpha_mT = alpha_T * 1e3;

    % --- Plot ---
    fig = figure('Position', [50 50 1000 1100], 'Color', 'w');
    hold on;

    uv_mm = samples_uv * 1e3;
    r_surf_mm = p.r_surf * 1e3;

    max_mag   = max(abs(alpha_mT));
    Ncolor    = 256;
    cmap      = jet(Ncolor);
    idx_a = min(Ncolor, max(1, round(abs(alpha_mT)/max_mag * (Ncolor-1)) + 1));

    % Arrow scale — max α → 1.2 mm arrow
    L_max_mm_arrow = 1.2;
    scale_arrow    = L_max_mm_arrow / max_mag;

    % --- α arrows (drawn along in-plane radial direction) ---
    for i = 1:N_samp
        sgn = sign(alpha_mT(i)); if sgn == 0, sgn = 1; end
        arrow_u = sgn * abs(alpha_mT(i)) * rhat_inplane_uv(i,1) * scale_arrow;
        arrow_v = sgn * abs(alpha_mT(i)) * rhat_inplane_uv(i,2) * scale_arrow;
        quiver(uv_mm(i,1), uv_mm(i,2), arrow_u, arrow_v, 0, ...
            'Color', cmap(idx_a(i),:), 'LineWidth', 2.0, 'MaxHeadSize', 0.6);
    end

    % Sample dots
    plot(uv_mm(:,1), uv_mm(:,2), 'k.', 'MarkerSize', 6);

    % Pole outline circle + axis center marker
    th = linspace(0, 2*pi, 200);
    plot(r_surf_mm*cos(th), r_surf_mm*sin(th), 'k-', 'LineWidth', 1.8);
    plot(0, 0, '+', 'Color', [0.3 0.3 0.3], 'MarkerSize', 10, 'LineWidth', 1.2);

    % Colorbar
    colormap(gca, jet(Ncolor));
    clim([0, max_mag]);
    cb = colorbar;
    cb.Label.String   = '|\alpha|  [mT]';
    cb.Label.FontSize = 12;
    cb.FontSize       = 11;

    % Text annotation
    tx =  0.95 * p.xylim;
    ty = -0.95 * p.xylim;
    text(tx, ty, '\alpha r̂  :  ⟂ pole surface  (from B = \alpha r̂ + \beta n̂)', ...
        'FontSize', 12, 'Interpreter', 'tex', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'EdgeColor', 'k', 'BackgroundColor', 'w', 'Margin', 6);

    xlabel('u [mm]', 'FontSize', 13);
    ylabel('v [mm]', 'FontSize', 13);
    title(sprintf('\\alpha r̂ component  at s=%g mm (%s, cross-section, Coil1)', ...
        p.s_anchor*1e3, p.section), ...
        'FontSize', 14, 'FontWeight', 'bold', 'Interpreter','tex');

    axis equal;
    xlim([-p.xylim, p.xylim]); ylim([-p.xylim, p.xylim]);
    set(gca, 'FontSize', 12, 'GridAlpha', 0.3);
    grid on; box on;

    % --- Save ---
    out_path = fullfile('..','..','figures','coil1', p.filename);
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, out_path, '-dpng', '-r150');
    fprintf('  Saved: %s\n', out_path);

    % --- Report (surface ring only, 12 angles) ---
    fprintf('  α at surface ring (12 angles):\n');
    for i = 2:2:N_samp      % every other sample = surface points
        fprintf('    θ=%3.0f°  (u=%+5.2f, v=%+5.2f) mm  α=%+7.2f mT\n', ...
            angles_deg(i/2), uv_mm(i,1), uv_mm(i,2), alpha_mT(i));
    end
    fprintf('  α range: [%+7.2f, %+7.2f] mT, |max|=%.2f\n', ...
        min(alpha_mT), max(alpha_mT), max(abs(alpha_mT)));
end
