%% plot_Bvector_P1_projections.m
%  Oblique (non-orthogonal) basis decomposition of B at selected stations inside P1.
%    r_hat : outward normal to pole surface (pure radial in cylinder; tilted in cone)
%    n_hat : from SAMPLE toward P1 tip (per-sample, varies across cross-section)
%  Samples taken along r_hat from the axis anchor to the surface — so for cone the
%  sample line is SLANTED (perpendicular to the tilted cone surface).
%  Arrow length is proportional to |coefficient|; colour encodes magnitude (redundant).
%
%  Two figures produced (one per section):
%    magnetic_sim/hung/figures/coil1/Bvector_P1_projections_RoundFillet_l500.png       (cylinder, s=25)
%    magnetic_sim/hung/figures/coil1/Bvector_P1_projections_RoundFillet_l500_cone.png  (cone,     s=10)

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));

%% 1. Load full-model FEM data + build 2D interpolants
fprintf('Loading RoundFillet l=500 full-model data...\n');
d = import_ansys_data( ...
    fullfile('..','..','results','coil1','round_filleted_conv'), ...
    'all', 'coil1');
N_total = length(d.bx);
fprintf('  total nodes: %d\n', N_total);

slab_half = 100e-6;
mask = abs(d.y) < slab_half;
N_slab = sum(mask);
fprintf('  slab |y|<%g um: %d nodes\n', slab_half*1e6, N_slab);

F_Bx = scatteredInterpolant(d.x(mask), d.z(mask), d.bx(mask), 'linear', 'nearest');
F_Bz = scatteredInterpolant(d.x(mask), d.z(mask), d.bz(mask), 'linear', 'nearest');

%% 2. Pole geometry + local unit vectors (xz plane, mm)
tip_x = 0.408;  tip_z = -0.289;
tilt  = -5.71;
a_hat  = [cosd(tilt); sind(tilt)];
n_perp = [-sind(tilt); cosd(tilt)];

POLE_R        = 3.175;
POLE_CONE_LEN = 15.875;
POLE_TOT_LEN  = 43.0;
CONE_HALF_ANG = 11.31;

r_hat_cyl_up  = n_perp;
r_hat_cyl_dn  = -n_perp;
r_hat_cone_up = cosd(CONE_HALF_ANG)*n_perp - sind(CONE_HALF_ANG)*a_hat;
r_hat_cone_dn = -cosd(CONE_HALF_ANG)*n_perp - sind(CONE_HALF_ANG)*a_hat;

%% 3. Pole silhouette (shared overlay)
s_c = linspace(0, POLE_CONE_LEN, 80);
s_y = linspace(POLE_CONE_LEN, POLE_TOT_LEN, 20);
n_c = s_c * tand(CONE_HALF_ANG);
n_y = POLE_R * ones(size(s_y));
s_all = [s_c s_y];
n_all = [n_c n_y];
cx = tip_x + s_all * a_hat(1);
cz = tip_z + s_all * a_hat(2);
ux = cx + n_all * n_perp(1);  uz = cz + n_all * n_perp(2);
lx = cx - n_all * n_perp(1);  lz = cz - n_all * n_perp(2);
sil_x = [ux fliplr(lx) ux(1)];
sil_z = [uz fliplr(lz) uz(1)];
end_x = tip_x + POLE_TOT_LEN * a_hat(1);
end_z = tip_z + POLE_TOT_LEN * a_hat(2);

%% 4. Plot configurations — one per section
N_rad = 2;                                      % axis + mid + surface

plots = struct();
% --- Cylinder: t_max = POLE_R (r_hat ⟂ axis) ---
plots(1).name      = 'cylinder';
plots(1).section   = 'P1 cylinder';
plots(1).s_anchor  = 25;
plots(1).rhat_up   = r_hat_cyl_up;
plots(1).rhat_dn   = r_hat_cyl_dn;
plots(1).t_max     = POLE_R;
plots(1).xlim_     = [22 29];
plots(1).ylim_     = [-8 2.5];
plots(1).filename  = 'Bvector_P1_projections_RoundFillet_l500.png';

% --- Cone: t_max = s_anchor * sin(θ_c)   (slanted sample line ⟂ cone surface) ---
plots(2).name      = 'cone';
plots(2).section   = 'P1 cone';
plots(2).s_anchor  = 10;
plots(2).rhat_up   = r_hat_cone_up;
plots(2).rhat_dn   = r_hat_cone_dn;
plots(2).t_max     = 10 * sind(CONE_HALF_ANG);   % ≈ 1.96 mm
plots(2).xlim_     = [7 14];
plots(2).ylim_     = [-5 3];
plots(2).filename  = 'Bvector_P1_projections_RoundFillet_l500_cone.png';

%% 5. Loop over plots
for k = 1:length(plots)
    p = plots(k);
    fprintf('\n=== %s (axis anchor s=%g mm, t_max=%.3f mm) ===\n', ...
        p.section, p.s_anchor, p.t_max);

    % --- Build samples along r_hat from axis anchor to surface ---
    anchor_x = tip_x + p.s_anchor*a_hat(1);
    anchor_z = tip_z + p.s_anchor*a_hat(2);

    samples = zeros(0, 4);
    for r_local = linspace(0, p.t_max, N_rad + 1)
        x_mm = anchor_x + r_local*p.rhat_up(1);
        z_mm = anchor_z + r_local*p.rhat_up(2);
        samples(end+1, :) = [x_mm, z_mm, p.rhat_up(1), p.rhat_up(2)]; %#ok<AGROW>
    end
    for r_local = linspace(p.t_max/N_rad, p.t_max, N_rad)
        x_mm = anchor_x + r_local*p.rhat_dn(1);
        z_mm = anchor_z + r_local*p.rhat_dn(2);
        samples(end+1, :) = [x_mm, z_mm, p.rhat_dn(1), p.rhat_dn(2)]; %#ok<AGROW>
    end

    N_samp = size(samples, 1);
    pos_x = samples(:, 1);
    pos_z = samples(:, 2);
    rhx   = samples(:, 3);
    rhz   = samples(:, 4);

    % --- Per-sample n_hat = unit vector to tip ---
    dx_to_tip = tip_x - pos_x;
    dz_to_tip = tip_z - pos_z;
    dist_tip  = sqrt(dx_to_tip.^2 + dz_to_tip.^2);
    nhx = dx_to_tip ./ dist_tip;
    nhz = dz_to_tip ./ dist_tip;

    % --- Query B + solve oblique decomposition ---
    Bx = F_Bx(pos_x*1e-3, pos_z*1e-3);
    Bz = F_Bz(pos_x*1e-3, pos_z*1e-3);
    D        = rhx.*nhz - rhz.*nhx;
    alpha    = (Bx.*nhz - Bz.*nhx) ./ D;
    beta     = (Bz.*rhx - Bx.*rhz) ./ D;
    alpha_mT = alpha * 1e3;
    beta_mT  = beta  * 1e3;
    ang_rn   = acosd(rhx.*nhx + rhz.*nhz);

    % --- Plot ---
    fig = figure('Position', [50 50 1000 1100], 'Color', 'w');

    L_max_mm = 1.5;
    max_mag  = max([max(abs(alpha_mT)), max(abs(beta_mT))]);
    scale    = L_max_mm / max_mag;               % shared mm/mT -> true magnitude ratio

    r_dx = alpha_mT .* rhx * scale;
    r_dz = alpha_mT .* rhz * scale;
    n_dx = beta_mT  .* nhx * scale;
    n_dz = beta_mT  .* nhz * scale;

    Ncolor = 256;
    cmap   = jet(Ncolor);
    idx_a = min(Ncolor, max(1, round(abs(alpha_mT)/max_mag * (Ncolor-1)) + 1));
    idx_b = min(Ncolor, max(1, round(abs(beta_mT)/max_mag  * (Ncolor-1)) + 1));

    hold on;

    % Grey reference rays sample -> tip
    for i = 1:N_samp
        plot([pos_x(i), tip_x], [pos_z(i), tip_z], ':', ...
            'Color', [0.80 0.80 0.80], 'LineWidth', 0.4);
    end

    % β arrows (thin) under
    for i = 1:N_samp
        quiver(pos_x(i), pos_z(i), n_dx(i), n_dz(i), 0, ...
            'Color', cmap(idx_b(i),:), 'LineWidth', 1.3, 'MaxHeadSize', 0.45);
    end
    % α arrows (thick) on top
    for i = 1:N_samp
        quiver(pos_x(i), pos_z(i), r_dx(i), r_dz(i), 0, ...
            'Color', cmap(idx_a(i),:), 'LineWidth', 2.8, 'MaxHeadSize', 0.8);
    end

    plot(pos_x, pos_z, 'k.', 'MarkerSize', 10);
    plot(sil_x, sil_z, 'k-', 'LineWidth', 1.6);
    plot([tip_x end_x], [tip_z end_z], '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0);

    colormap(gca, jet(Ncolor));
    clim([0, max_mag]);
    cb = colorbar;
    cb.Label.String   = 'magnitude [mT]';
    cb.Label.FontSize = 12;
    cb.FontSize       = 11;

    tx = p.xlim_(2) - 0.03*(p.xlim_(2) - p.xlim_(1));
    ty = p.ylim_(1) + 0.03*(p.ylim_(2) - p.ylim_(1));
    text(tx, ty, {'\alpha r̂  :  ⟂ pole surface', '\beta n̂  :  toward P1 tip'}, ...
        'FontSize', 12, 'Interpreter', 'tex', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'EdgeColor', 'k', 'BackgroundColor', 'w', 'Margin', 6);

    xlabel('x [mm]', 'FontSize', 13);
    ylabel('z [mm]', 'FontSize', 13);
    title(sprintf('B = \\alpha r̂ + \\beta n̂   at s=%g mm (%s, Coil1)', p.s_anchor, p.section), ...
        'FontSize', 14, 'FontWeight', 'bold', 'Interpreter','tex');

    axis equal;
    xlim(p.xlim_); ylim(p.ylim_);
    set(gca, 'FontSize', 12, 'GridAlpha', 0.3);
    grid on; box on;

    % --- Save ---
    out_path = fullfile('..','..','figures','coil1', p.filename);
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, out_path, '-dpng', '-r150');
    fprintf('  Saved: %s\n', out_path);

    fprintf('  N_samp=%d,  scale=%.4f mm/mT (L_max=%.2f mm)\n', N_samp, scale, L_max_mm);
    fprintf('  alpha (r_hat coef): min=%+7.2f mT, max=%+7.2f mT, |max|=%.2f\n', ...
        min(alpha_mT), max(alpha_mT), max(abs(alpha_mT)));
    fprintf('  beta  (n_hat coef): min=%+7.2f mT, max=%+7.2f mT, |max|=%.2f\n', ...
        min(beta_mT),  max(beta_mT),  max(abs(beta_mT)));
    fprintf('  r_hat angle to n_hat: [%.2f, %.2f] deg\n', min(ang_rn), max(ang_rn));
end
