%% plot_Bcontour_xaza_baseline.m
%  |B| contour on the x_a-z_a plane (y_a = 0) for Coil1 BASELINE run
%  Baseline: D-shape pole (no fillet), no NREFINE (Apr 2 run).
%  Same style/view as plot_Bcontour_xaza.m for direct comparison.
%
%  Actuator frame: x_a -> P1, y_a -> P3, z_a -> P5
%
%  Data source: magnetic_sim/hung/results/coil1/baseline/  (368,686 nodes)
%  Output:      magnetic_sim/hung/figures/coil1/Bcontour_xaza_Dshape_baseline.png

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));

%% 1. Actuator frame basis (ANSYS coords)
theta = 54.74 * pi/180;
s = sin(theta); c = cos(theta);
xa_hat = [ s*cosd(0);    s*sind(0);    -c ];
ya_hat = [ s*cosd(120);  s*sind(120);  -c ];
za_hat = [ s*cosd(60);   s*sind(60);   +c ];
R = [xa_hat'; ya_hat'; za_hat'];

%% 2. Load baseline (no NREFINE) nodes + B-field, rotate to actuator frame
d = import_ansys_data(fullfile('..','..','results','coil1','baseline'), 'all', 'coil1');

P_act = R * [d.x, d.y, d.z]';
x_a = P_act(1,:)';
y_a = P_act(2,:)';
z_a = P_act(3,:)';

bsum = sqrt(d.bx.^2 + d.by.^2 + d.bz.^2);
N_total = length(bsum);

%% 3. Node-count report
r_ansys = sqrt(d.x.^2 + d.y.^2 + d.z.^2);
N_r300  = sum(r_ansys < 300e-6);
N_cube50 = sum(abs(d.x)<50e-6 & abs(d.y)<50e-6 & abs(d.z)<50e-6);

fprintf('=== Baseline (no NREFINE) node count ===\n');
fprintf('  Total nodes (whole model)        : %d\n', N_total);
fprintf('  R < 300 um (figure visual core)  : %d\n', N_r300);
fprintf('  +/-50 um cube                    : %d\n', N_cube50);

%% 4. 3D linear interpolant, query exact y_a = 0 plane
%  Restrict to nodes inside the air sphere (R < 0.6 mm in ANSYS coords) to
%  exclude iron-body nodes whose very-high |B| otherwise leaks into the plot
%  window via Delaunay tets that span air+iron.
r_ansys_all = sqrt(d.x.^2 + d.y.^2 + d.z.^2);
air_mask = r_ansys_all < 0.6e-3;
fprintf('\nNodes inside air sphere (R<0.6mm) used for interp: %d\n', sum(air_mask));

fprintf('Building 3D scatteredInterpolant...\n');
tic;
F3D = scatteredInterpolant(x_a(air_mask), y_a(air_mask), z_a(air_mask), bsum(air_mask), ...
    'linear', 'nearest');
fprintf('  built in %.1f s\n', toc);

Ngrid = 401;
xi = linspace(-300e-6, 300e-6, Ngrid);
zi = linspace(-300e-6, 300e-6, Ngrid);
[Xi, Zi] = meshgrid(xi, zi);

fprintf('Querying %d x %d = %d points at y_a=0 ...\n', Ngrid, Ngrid, Ngrid^2);
tic;
Bi_raw = F3D(Xi, zeros(size(Xi)), Zi) * 1e3;
fprintf('  queried in %.1f s\n', toc);

B_WP = F3D(0, 0, 0) * 1e3;
fprintf('WP |B| = %.2f mT\n', B_WP);

%% 5. Light Gaussian smoothing (same as smooth version)
Bi = imgaussfilt(Bi_raw, 1.0);

%% 6. Plot (IDENTICAL style/view as plot_Bcontour_xaza.m)
fig = figure('Position',[100 100 950 820], 'Color','w');
contourf(Xi*1e6, Zi*1e6, Bi, 24, 'LineStyle','none');
colormap(jet(256));
caxis([4 20]);                          % match smooth version colorbar range
cb = colorbar;
cb.Label.String   = '|B| [mT]';
cb.Label.FontSize = 14;
cb.FontSize       = 12;

hold on;
plot(0, 0, 'k+', 'MarkerSize', 16, 'LineWidth', 2.5);
text(30, -30, sprintf('%.2f mT', B_WP), ...
    'FontSize', 13, 'FontWeight', 'bold', ...
    'BackgroundColor', [0.95 1.0 0.95], 'EdgeColor', 'none', 'Margin', 3);

% Node count annotation (same placement as smooth version)
text(-290, -280, sprintf('N_{nodes} = %s', addcommas(N_total)), ...
    'FontSize', 10, 'Color', [0.2 0.2 0.2], ...
    'BackgroundColor', [1 1 1 0.7], 'EdgeColor', 'none', 'Margin', 2);

xlabel('x_a [\mum]', 'FontSize', 14);
ylabel('z_a [\mum]', 'FontSize', 14);
title('|B| on x_a-z_a plane (y_a=0) — D-shape (no fillet, no NREFINE)', ...
    'FontSize', 15, 'FontWeight', 'bold', 'Interpreter','tex');
axis equal;
xlim([-300 300]); ylim([-300 300]);
set(gca, 'FontSize', 12, 'GridAlpha', 0.3);
grid on; box on;

%% 7. Save
out_path = fullfile('..','..','figures','coil1','Bcontour_xaza_Dshape_baseline.png');
set(fig, 'PaperPositionMode', 'auto');
print(fig, out_path, '-dpng', '-r150');
fprintf('\nSaved: %s\n', out_path);

%% --- Local helper ---
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
