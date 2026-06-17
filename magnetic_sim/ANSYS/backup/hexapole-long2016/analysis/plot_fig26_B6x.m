%% plot_fig26_B6x.m — Reproduce dissertation Fig 2.6 for Method [B-6x]
%
%  Fig 2.6: Validation of the hexapole magnetic field model
%    (a) comparison of magnetic induction vectors (FEM vs model)
%    (b) normalized norms of error vectors (percent fitting error)
%    (c) definition of fitting error (text annotation)

%% 1. Load [B-6x] fit results
fit = load(fullfile('..', 'data', 'all6_bias_fit.mat'));
cnst = mt_constants();

fprintf('Loaded [B-6x] fit: ell=%.1f um, R_a=%.3e, C=%.4e\n', ...
    fit.ell_fixed*1e6, fit.R_a, fit.C);

%% 2. Load FEM data (same processing as fit script)
K_I = eye(6) - ones(6)/6;
apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];
cube_half = 50e-6;

% Build I_diss
I_diss = zeros(6,1);
for k = 1:6
    paper_idx = apdl_to_paper_idx(k);
    coil_sign = 2*cnst.pole_is_lower(paper_idx) - 1;
    I_diss(paper_idx) = coil_sign;
end
w = K_I * I_diss;

% Superposition of 6 single-coil data
d1 = import_ansys_data(fullfile('..', 'results', 'coil1'), 'wp', 'coil1');
[air_mask, ~] = filter_iron_nodes(d1.x, d1.y, d1.z, cnst, struct('visualize', false));
z_wp = d1.z - cnst.SPH_OFST;
mask = air_mask & abs(d1.x) < cube_half & abs(d1.y) < cube_half & abs(z_wp) < cube_half;

p_wp = [d1.x(mask), d1.y(mask), z_wp(mask)];
N = sum(mask);

bx_sum = zeros(N, 1);
by_sum = zeros(N, 1);
bz_sum = zeros(N, 1);
for k = 1:6
    coil_name = sprintf('coil%d', k);
    dk = import_ansys_data(fullfile('..', 'results', coil_name), 'wp', coil_name);
    bx_sum = bx_sum + dk.bx(mask);
    by_sum = by_sum + dk.by(mask);
    bz_sum = bz_sum + dk.bz(mask);
end
b_fem = [bx_sum; by_sum; bz_sum];

fprintf('Loaded %d FEM nodes in 100 um cube\n', N);

%% 3. Compute model B-field
pos = fit.pos;
k_m = cnst.k_m;
bx_m = zeros(N,1); by_m = zeros(N,1); bz_m = zeros(N,1);
for i = 1:6
    dx = p_wp(:,1) - pos(1,i);
    dy = p_wp(:,2) - pos(2,i);
    dz = p_wp(:,3) - pos(3,i);
    r3 = (dx.^2 + dy.^2 + dz.^2).^(3/2);
    bx_m = bx_m + (-w(i)) * dx ./ r3;
    by_m = by_m + (-w(i)) * dy ./ r3;
    bz_m = bz_m + (-w(i)) * dz ./ r3;
end
bx_model = fit.C * k_m * bx_m;
by_model = fit.C * k_m * by_m;
bz_model = fit.C * k_m * bz_m;

%% 4. Compute per-node error
bx_f = bx_sum; by_f = by_sum; bz_f = bz_sum;
err_x = bx_model - bx_f;
err_y = by_model - by_f;
err_z = bz_model - bz_f;
err_mag = sqrt(err_x.^2 + err_y.^2 + err_z.^2);
fem_mag = sqrt(bx_f.^2 + by_f.^2 + bz_f.^2);
valid = fem_mag > 1e-6 * max(fem_mag);
pct_err = 100 * err_mag ./ fem_mag;
pct_err(~valid) = 0;

fprintf('Error stats: mean=%.2f%%, median=%.2f%%, 95th=%.2f%%, max=%.2f%%\n', ...
    mean(pct_err(valid)), median(pct_err(valid)), ...
    prctile(pct_err(valid), 95), max(pct_err(valid)));

%% 5. Create figure (matching dissertation Fig 2.6 style)
fig = figure('Position', [100 100 1100 500], 'Color', 'w');

% --- Panel (a): 3D vector comparison ---
ax1 = subplot(1, 2, 1);

% Convert to um for display
x_um = p_wp(:,1) * 1e6;
y_um = p_wp(:,2) * 1e6;
z_um = p_wp(:,3) * 1e6;

% Subsample for clarity (every 3rd point)
skip = 3;
idx = 1:skip:N;

% Scale factor for arrows
scale_fem = 1e6;  % will be normalized by quiver auto-scaling

% FEM vectors (blue)
q1 = quiver3(x_um(idx), y_um(idx), z_um(idx), ...
    bx_f(idx)*scale_fem, by_f(idx)*scale_fem, bz_f(idx)*scale_fem, ...
    1.5, 'Color', [0 0.2 0.8], 'LineWidth', 0.6);
hold on;

% Model vectors (red)
q2 = quiver3(x_um(idx), y_um(idx), z_um(idx), ...
    bx_model(idx)*scale_fem, by_model(idx)*scale_fem, bz_model(idx)*scale_fem, ...
    1.5, 'Color', [0.85 0 0], 'LineWidth', 0.6);

xlabel('x (\mum)', 'FontSize', 11);
ylabel('y (\mum)', 'FontSize', 11);
zlabel('z (\mum)', 'FontSize', 11);
legend([q1, q2], {'Flux Density (FEM)', 'Flux Density (fit)'}, ...
    'Location', 'northeast', 'FontSize', 9);
title('(a)', 'FontSize', 12, 'FontWeight', 'bold');
axis equal;
xlim([-50 50]); ylim([-50 50]); zlim([-50 50]);
view(-37.5, 30);
grid on;
set(ax1, 'FontSize', 10);

% --- Panel (b): Percent fitting error ---
ax2 = subplot(1, 2, 2);

% Sort by index (just use sequential)
plot(1:sum(valid), pct_err(valid), '.', 'Color', [0 0.2 0.8], 'MarkerSize', 4);
hold on;

% Add mean line
yline(mean(pct_err(valid)), '--r', sprintf('Mean = %.2f%%', mean(pct_err(valid))), ...
    'LineWidth', 1.2, 'FontSize', 9, 'LabelHorizontalAlignment', 'left');

xlabel('Index of Points', 'FontSize', 11);
ylabel('Percent Fitting Error (%)', 'FontSize', 11);
title('(b)', 'FontSize', 12, 'FontWeight', 'bold');
set(ax2, 'FontSize', 10);
grid on;
xlim([0 sum(valid)]);

% Add annotation for error definition and parameters
dim = [0.55 0.02 0.4 0.15];
str = {sprintf('(c) Error %% = |\\DeltaB| / |B_{FEM}| x 100'), ...
    sprintf('[B-6x]: ell = %.1f um,  R_a = %.2e A/Wb', fit.ell_fixed*1e6, fit.R_a), ...
    sprintf('Mean err = %.2f%%,  Max = %.2f%%', mean(pct_err(valid)), max(pct_err(valid)))};
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', ...
    'BackgroundColor', [1 1 0.9], 'EdgeColor', [0.5 0.5 0.5], ...
    'FontSize', 9, 'Interpreter', 'tex');

% Main title
sgtitle('Validation of the Hexapole Magnetic Field Model — Method [B-6x]', ...
    'FontSize', 13, 'FontWeight', 'bold');

%% 6. Save figure
saveas(fig, fullfile('..', 'figures', 'fig26_B6x_validation.png'));
fprintf('\nFigure saved to figures/fig26_B6x_validation.png\n');
