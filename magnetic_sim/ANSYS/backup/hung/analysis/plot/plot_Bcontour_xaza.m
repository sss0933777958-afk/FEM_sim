%% plot_Bcontour_xaza.m
%  |B| contour on the x_a-z_a plane (y_a = 0) for Coil1 (D-shape + R20um fillet)
%  Actuator frame (x_a, y_a, z_a) aligned with pole directions:
%    x_a -> P1 (lower, azimuth 0 deg)
%    y_a -> P3 (lower, azimuth 120 deg)
%    z_a -> P5 (upper, azimuth 60 deg)
%  Magic angle 54.74 deg guarantees orthogonality; right-handed (x_a x y_a = z_a).
%
%  Smooth version: 3D scatteredInterpolant at y_a=0 (no slab), N=401 grid,
%                  + light Gaussian smoothing (sigma=1 px).
%
%  Data source: magnetic_sim/ANSYS/backup/hung/results/coil1/filleted/  (3-pass NREFINE)
%  Output:      magnetic_sim/ANSYS/backup/hung/figures/coil1/Bcontour_xaza_Dfillet_smooth.png

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));

%% 1. Actuator frame basis (ANSYS coords)
theta = 54.74 * pi/180;   % magic angle (acosd(1/sqrt(3)))
s = sin(theta); c = cos(theta);

xa_hat = [ s*cosd(0);    s*sind(0);    -c ];   % P1  (lower)
ya_hat = [ s*cosd(120);  s*sind(120);  -c ];   % P3  (lower)
za_hat = [ s*cosd(60);   s*sind(60);   +c ];   % P5  (upper)

% Rotation matrix ANSYS -> actuator (rows = basis vectors)
R = [xa_hat'; ya_hat'; za_hat'];

%% 2. Load Coil1 nodes + B-field, rotate to actuator frame
d = import_ansys_data(fullfile('..','..','results','coil1','filleted'), 'all', 'coil1');

P_act = R * [d.x, d.y, d.z]';      % 3 x N
x_a = P_act(1,:)';
y_a = P_act(2,:)';
z_a = P_act(3,:)';

bsum = sqrt(d.bx.^2 + d.by.^2 + d.bz.^2);

N_total = length(bsum);

%% 3. Node-count report
r_ansys = sqrt(d.x.^2 + d.y.^2 + d.z.^2);
N_r300  = sum(r_ansys < 300e-6);
N_cube50 = sum(abs(d.x)<50e-6 & abs(d.y)<50e-6 & abs(d.z)<50e-6);

fprintf('=== Node count summary ===\n');
fprintf('  Total nodes (whole model)        : %d\n', N_total);
fprintf('  R < 300 um (figure visual core)  : %d\n', N_r300);
fprintf('  +/-50 um cube (fitting region)   : %d\n', N_cube50);

%% 4. 3D scatteredInterpolant, query exact y_a = 0 plane
fprintf('\nBuilding 3D scatteredInterpolant (Delaunay on %d nodes)...\n', N_total);
tic;
F3D = scatteredInterpolant(x_a, y_a, z_a, bsum, 'linear', 'nearest');
fprintf('  built in %.1f s\n', toc);

Ngrid = 401;
xi = linspace(-300e-6, 300e-6, Ngrid);
zi = linspace(-300e-6, 300e-6, Ngrid);
[Xi, Zi] = meshgrid(xi, zi);

fprintf('Querying %d x %d = %d points at y_a=0 ...\n', Ngrid, Ngrid, Ngrid^2);
tic;
Bi_raw = F3D(Xi, zeros(size(Xi)), Zi) * 1e3;   % mT
fprintf('  queried in %.1f s\n', toc);

B_WP = F3D(0, 0, 0) * 1e3;
fprintf('WP |B| = %.2f mT\n', B_WP);

%% 5. Light Gaussian smoothing (sigma = 1 pixel ~ 1.5 um real space)
Bi = imgaussfilt(Bi_raw, 1.0);

%% 6. Plot (match reference style)
fig = figure('Position',[100 100 950 820], 'Color','w');
contourf(Xi*1e6, Zi*1e6, Bi, 24, 'LineStyle','none');
colormap(jet(256));
cb = colorbar;
cb.Label.String   = '|B| [mT]';
cb.Label.FontSize = 14;
cb.FontSize       = 12;

hold on;
plot(0, 0, 'k+', 'MarkerSize', 16, 'LineWidth', 2.5);
text(30, -30, sprintf('%.2f mT', B_WP), ...
    'FontSize', 13, 'FontWeight', 'bold', ...
    'BackgroundColor', [0.95 1.0 0.95], 'EdgeColor', 'none', 'Margin', 3);

% Node count annotation (bottom-left, small)
text(-290, -280, sprintf('N_{nodes} = %s', addcommas(N_total)), ...
    'FontSize', 10, 'Color', [0.2 0.2 0.2], ...
    'BackgroundColor', [1 1 1 0.7], 'EdgeColor', 'none', 'Margin', 2);

xlabel('x_a [\mum]', 'FontSize', 14);
ylabel('z_a [\mum]', 'FontSize', 14);
title('|B| on x_a-z_a plane (y_a=0) — D-shape + R20\mum fillet', ...
    'FontSize', 15, 'FontWeight', 'bold', 'Interpreter','tex');
axis equal;
xlim([-300 300]); ylim([-300 300]);
set(gca, 'FontSize', 12, 'GridAlpha', 0.3);
grid on; box on;

%% 7. Save (new filename, do NOT overwrite original)
out_path = fullfile('..','..','figures','coil1','Bcontour_xaza_Dfillet_smooth.png');
set(fig, 'PaperPositionMode', 'auto');
print(fig, out_path, '-dpng', '-r150');
fprintf('\nSaved: %s\n', out_path);

%% --- Local helper ---
function s = addcommas(n)
    % format integer with thousands separators (526645 -> '526,645')
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
