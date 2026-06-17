%% plot_Bvector_P1_axis.m
%  B-field vector plot on ANSYS y=0 plane (contains P1 axis)
%  RoundFillet l=500 model, coil1 excited, 4-pass NREFINE (839k nodes)
%
%  View: ANSYS coordinates (x horizontal, z vertical)
%  Range: x = 0..45 mm (covers WP to P1 pole end), z = -10..+10 mm
%  P1 axis: tilted 5.71° below horizontal, from tip (+0.408, 0, -0.289) mm
%           extends outward along (0.995, 0, -0.0995) to end (+43.2, 0, -4.57) mm
%
%  Data source: magnetic_sim/hung/results/coil1/round_filleted_conv/
%  Output:      magnetic_sim/hung/figures/coil1/Bvector_xz_P1_RoundFillet_l500.png

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));

%% 1. Load full-model FEM nodes + B-field
fprintf('Loading RoundFillet l=500 full-model data...\n');
d = import_ansys_data( ...
    fullfile('..','..','results','coil1','round_filleted_conv'), ...
    'all', 'coil1');
N_total = length(d.bx);
fprintf('  total nodes: %d\n', N_total);

%% 2. Slab filter: keep nodes with |y| < slab_half for y=0 plane approx
slab_half = 100e-6;                       % 100 um half-width
mask = abs(d.y) < slab_half;
N_slab = sum(mask);
fprintf('  slab |y|<%g um: %d nodes\n', slab_half*1e6, N_slab);

x_slab = d.x(mask);
z_slab = d.z(mask);
Bx_slab = d.bx(mask);
Bz_slab = d.bz(mask);
By_slab = d.by(mask);
Bmag_slab = sqrt(Bx_slab.^2 + By_slab.^2 + Bz_slab.^2);

%% 3. Build 2D interpolants (scattered) for Bx, Bz, |B| on y=0 plane
fprintf('Building 2D scatteredInterpolants...\n');
tic;
F_Bx   = scatteredInterpolant(x_slab, z_slab, Bx_slab, 'linear', 'nearest');
F_Bz   = scatteredInterpolant(x_slab, z_slab, Bz_slab, 'linear', 'nearest');
F_Bmag = scatteredInterpolant(x_slab, z_slab, Bmag_slab, 'linear', 'nearest');
fprintf('  built in %.1f s\n', toc);

%% 4. Sparse grid for quiver arrows (only one grid needed now)
Nx_q = 35; Nz_q = 15;                       % reduced density
x_q = linspace(0, 45e-3, Nx_q);
z_q = linspace(-10e-3, 10e-3, Nz_q);
[Xq, Zq] = meshgrid(x_q, z_q);

fprintf('Querying Bx, Bz on %d x %d quiver grid...\n', Nx_q, Nz_q);
tic;
Bxq = F_Bx(Xq, Zq);
Bzq = F_Bz(Xq, Zq);
fprintf('  done in %.1f s\n', toc);

% Per-arrow magnitude (|B_xz|, 2D component on cut plane)
Bmag_q = sqrt(Bxq.^2 + Bzq.^2) * 1e3;   % mT
Bnorm  = sqrt(Bxq.^2 + Bzq.^2);
Bnorm(Bnorm<1e-12) = 1e-12;

%% 5a. Compute pole silhouette polygon (used for masking arrows inside iron)
tip_x = 0.408; tip_z = -0.289;
tilt  = -5.71;
a_hat = [cosd(tilt); sind(tilt)];
n_hat = [-sind(tilt); cosd(tilt)];
POLE_R        = 3.175;
POLE_CONE_LEN = 15.875;
POLE_TOT_LEN  = 43.0;
CONE_HALF_ANG = 11.31;

s_cone = linspace(0, POLE_CONE_LEN, 80);
s_cyl  = linspace(POLE_CONE_LEN, POLE_TOT_LEN, 20);
n_cone = s_cone * tand(CONE_HALF_ANG);
n_cyl  = POLE_R * ones(size(s_cyl));
s_all  = [s_cone s_cyl];
n_all  = [n_cone n_cyl];

cx = tip_x + s_all * a_hat(1);
cz = tip_z + s_all * a_hat(2);
ux = cx + n_all * n_hat(1);  uz = cz + n_all * n_hat(2);
lx = cx - n_all * n_hat(1);  lz = cz - n_all * n_hat(2);
sil_x = [ux fliplr(lx) ux(1)];
sil_z = [uz fliplr(lz) uz(1)];

%% 5b. Plot — pure quiver, arrows coloured by magnitude, UNIFORM length
%      (arrows inside pole are kept — user wants to see B flow in iron)
fig = figure('Position',[100 100 1800 900], 'Color','w');

% Unit vectors (all arrows same physical length)
Ux = Bxq ./ Bnorm;
Uz = Bzq ./ Bnorm;

% Arrow length in axis units (mm)
arrow_len = 0.7;
Ux_scaled = Ux * arrow_len;
Uz_scaled = Uz * arrow_len;

% Colour binning
Ncolor = 64;
cmap   = jet(Ncolor);
clim_mT = [0, 50];
idx = round((Bmag_q - clim_mT(1)) / (clim_mT(2) - clim_mT(1)) * (Ncolor-1)) + 1;
idx = max(1, min(Ncolor, idx));
idx(isnan(Bmag_q)) = NaN;

hold on;
% Loop over colour bins — scale=0 means "use U,V as literal displacements"
for i = 1:Ncolor
    mask = (idx == i);
    if any(mask(:))
        quiver(Xq(mask)*1e3, Zq(mask)*1e3, ...
               Ux_scaled(mask), Uz_scaled(mask), ...
               0, 'Color', cmap(i,:), 'LineWidth', 1.1, ...
               'MaxHeadSize', 0.9);
    end
end

colormap(jet(256));
clim(clim_mT);
cb = colorbar;
cb.Label.String   = '|B_{xz}| [mT]  (arrow magnitude)';
cb.Label.FontSize = 14;
cb.FontSize       = 12;

% --- P1 pole silhouette on y=0 plane (uses sil_x, sil_z computed above) ---
% Outline only — arrows inside pole remain visible
plot(sil_x, sil_z, 'k-', 'LineWidth', 1.6);
% P1 axis dashed centerline
end_x = tip_x + POLE_TOT_LEN * a_hat(1);
end_z = tip_z + POLE_TOT_LEN * a_hat(2);
plot([tip_x end_x], [tip_z end_z], '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0);

% Markers
plot(tip_x, tip_z, 'k+', 'MarkerSize', 12, 'LineWidth', 2);
text(tip_x+0.3, tip_z+0.7, 'P1 tip', 'Color','k', 'FontSize',11, 'FontWeight','bold');
plot(0, 0, 'k+', 'MarkerSize', 14, 'LineWidth', 2);
text(0.5, 0.7, 'WP', 'Color','k', 'FontSize',11, 'FontWeight','bold');

xlabel('x [mm]  (ANSYS)', 'FontSize', 14);
ylabel('z [mm]  (ANSYS)', 'FontSize', 14);
title('B-field vectors on y=0 plane (along P1 axis) — RoundFillet l=500\mum, Coil1', ...
    'FontSize', 14, 'FontWeight', 'bold', 'Interpreter','tex');
axis equal;
xlim([0 45]); ylim([-10 10]);
set(gca, 'FontSize', 12, 'GridAlpha', 0.3);
grid on; box on;

% Node count annotation
text(1, -9.5, sprintf('N_{slab}=%s  (|y|<%.0f\\mum, from %s total)', ...
    addcommas(N_slab), slab_half*1e6, addcommas(N_total)), ...
    'FontSize', 10, 'Color', [0.2 0.2 0.2], ...
    'BackgroundColor', [1 1 1 0.8], 'Margin', 2);

%% 6. Save
out_path = fullfile('..','..','figures','coil1','Bvector_xz_P1_RoundFillet_l500.png');
set(fig, 'PaperPositionMode', 'auto');
print(fig, out_path, '-dpng', '-r150');
fprintf('\nSaved: %s\n', out_path);

% Report WP |B|
B_WP = F_Bmag(0, 0) * 1e3;
fprintf('WP |B| = %.3f mT\n', B_WP);

%% --- helper ---
function s = addcommas(n)
    str = num2str(n);
    L = length(str);
    s = '';
    for i = 1:L
        s = [s str(i)];
        if mod(L - i, 3) == 0 && i < L
            s = [s ','];
        end
    end
end
