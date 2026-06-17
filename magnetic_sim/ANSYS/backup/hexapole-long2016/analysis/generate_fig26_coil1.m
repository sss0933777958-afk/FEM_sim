%% generate_fig26_coil1.m — Reproduce fig26_B6x_clean style for coil1 [J]
%  Panel (a): 3D arrows - FEM (blue) vs Model (red)
%  Panel (b): Scatter - per-node fitting error with mean dashed line

%% 1. Load data
cnst = mt_constants();
K_I = eye(6) - ones(6)/6;

d = import_ansys_data(fullfile('..', 'results', 'coil1'), 'wp', 'coil1');
[air_mask, ~] = filter_iron_nodes(d.x, d.y, d.z, cnst, struct('visualize', false));
z_wp = d.z - cnst.SPH_OFST;

cube_half = 50e-6;
mask_cube = air_mask & abs(d.x) < cube_half & abs(d.y) < cube_half & abs(z_wp) < cube_half;
p_cube = [d.x(mask_cube), d.y(mask_cube), z_wp(mask_cube)];
bx_cube = d.bx(mask_cube); by_cube = d.by(mask_cube); bz_cube = d.bz(mask_cube);
N_cube = sum(mask_cube);

I_vec = zeros(6,1); I_vec(1) = 1;
w = K_I * I_vec;

J = load(fullfile('..', 'data', 'joint_6coil_19param_fit.mat'));
pos_J = J.best.pos;
b_unit = eval_charge_field(pos_J, w, cnst.k_m, p_cube);
C_k1 = J.best.C_k(1);
bx_mod = C_k1 * b_unit(1:N_cube);
by_mod = C_k1 * b_unit(N_cube+1:2*N_cube);
bz_mod = C_k1 * b_unit(2*N_cube+1:3*N_cube);

fprintf('Cube nodes: %d\n', N_cube);

%% 2. Grid subsample
spacing = 18e-6;
xg = -cube_half:spacing:cube_half;
yg = -cube_half:spacing:cube_half;
zg = -cube_half:spacing:cube_half;

cell_ids = zeros(N_cube, 3);
for ii = 1:N_cube
    [~, cell_ids(ii,1)] = min(abs(xg - p_cube(ii,1)));
    [~, cell_ids(ii,2)] = min(abs(yg - p_cube(ii,2)));
    [~, cell_ids(ii,3)] = min(abs(zg - p_cube(ii,3)));
end
[~, sel] = unique(cell_ids, 'rows', 'first');
fprintf('Grid-subsampled: %d nodes (spacing=%d um)\n', length(sel), spacing*1e6);

px = p_cube(sel,1)*1e6; py = p_cube(sel,2)*1e6; pz = p_cube(sel,3)*1e6;
ufx = bx_cube(sel); ufy = by_cube(sel); ufz = bz_cube(sel);
umx = bx_mod(sel);  umy = by_mod(sel);  umz = bz_mod(sel);

% Normalize to uniform arrow length
arrow_len = 10;  % um
for ii = 1:length(sel)
    mag_f = norm([ufx(ii), ufy(ii), ufz(ii)]);
    if mag_f > 0
        ufx(ii) = arrow_len*ufx(ii)/mag_f;
        ufy(ii) = arrow_len*ufy(ii)/mag_f;
        ufz(ii) = arrow_len*ufz(ii)/mag_f;
    end
    mag_m = norm([umx(ii), umy(ii), umz(ii)]);
    if mag_m > 0
        umx(ii) = arrow_len*umx(ii)/mag_m;
        umy(ii) = arrow_len*umy(ii)/mag_m;
        umz(ii) = arrow_len*umz(ii)/mag_m;
    end
end

%% 3. Panel (b) data
b_fem_mag = sqrt(bx_cube.^2 + by_cube.^2 + bz_cube.^2);
err_mag = sqrt((bx_mod-bx_cube).^2 + (by_mod-by_cube).^2 + (bz_mod-bz_cube).^2);
rel_err = err_mag ./ b_fem_mag * 100;
valid = b_fem_mag > 1e-6 * max(b_fem_mag);
rel_err(~valid) = NaN;
mean_err = mean(rel_err(valid), 'omitnan');
fprintf('Mean error: %.2f%%, Max error: %.2f%%\n', mean_err, max(rel_err(valid)));

%% 4. Create figure — two near-square subplots
fig = figure('Position', [50, 50, 1400, 680], 'Color', 'w');

col_fem   = [0 0.447 0.741];   % blue
col_model = [1 0 0];           % red

%% Panel (a) — near-square 3D quiver
ax1 = axes('Position', [0.04, 0.06, 0.46, 0.84]);
hold on;

% Arrow parameters — match reference: thin shaft, moderate head
head_len = 3.0;   % um
head_wid = 1.3;   % um half-width
lw = 1.0;         % shaft line width

% Draw interleaved
for ii = 1:length(sel)
    draw_arrow(ax1, px(ii),py(ii),pz(ii), ufx(ii),ufy(ii),ufz(ii), col_fem, head_len,head_wid,lw);
    draw_arrow(ax1, px(ii),py(ii),pz(ii), umx(ii),umy(ii),umz(ii), col_model, head_len,head_wid,lw);
end

% Legend
h1 = plot3(NaN,NaN,NaN, '-', 'Color', col_fem, 'LineWidth', 2.5);
h2 = plot3(NaN,NaN,NaN, '-', 'Color', col_model, 'LineWidth', 2.5);
legend([h1 h2], {'FEM','Model'}, 'Location','northeast', ...
    'FontSize', 13, 'FontWeight', 'bold', 'LineWidth', 1.5);

hold off;
grid on;
set(ax1, 'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.7);
xlim([-55 55]); ylim([-55 55]); zlim([-55 55]);
set(ax1, 'XTick', -50:50:50, 'YTick', -50:50:50, 'ZTick', -50:25:50);
xlabel('x (\mum)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('y (\mum)', 'FontSize', 14, 'FontWeight', 'bold');
zlabel('z (\mum)', 'FontSize', 14, 'FontWeight', 'bold');
set(ax1, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 2.5, 'Box', 'on');
view(135, 20);

% (a) label
annotation(fig, 'textbox', [0.01, 0.91, 0.05, 0.08], ...
    'String', '(a)', 'FontSize', 20, 'FontWeight', 'bold', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left');

%% Panel (b) — near-square scatter
ax2 = axes('Position', [0.58, 0.12, 0.38, 0.78]);
hold on;

scatter(1:sum(valid), rel_err(valid), 10, [0.1 0.35 0.75], 'filled', 'MarkerFaceAlpha', 0.6);
yline(mean_err, '--', 'Color', [0.9 0.2 0.1], 'LineWidth', 2.5);

hold off;
grid on;
set(ax2, 'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.7);
xlabel('Index of Points', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Fitting Error (%)', 'FontSize', 14, 'FontWeight', 'bold');
set(ax2, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 2.5, 'Box', 'on');
xlim([0 N_cube+10]);

% (b) label
annotation(fig, 'textbox', [0.55, 0.91, 0.05, 0.08], ...
    'String', '(b)', 'FontSize', 20, 'FontWeight', 'bold', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left');

%% Save
out_path = fullfile('..', 'figures', 'fig26_coil1_J_final.png');
exportgraphics(fig, out_path, 'Resolution', 200);
fprintf('Saved to %s\n', out_path);
close(fig);

%% ===== Local functions =====

function b = eval_charge_field(pos, w, k_m, p_wp)
    N = size(p_wp, 1);
    bx = zeros(N,1); by = zeros(N,1); bz = zeros(N,1);
    for i = 1:6
        dx = p_wp(:,1)-pos(1,i); dy = p_wp(:,2)-pos(2,i); dz = p_wp(:,3)-pos(3,i);
        r3 = (dx.^2+dy.^2+dz.^2).^(3/2);
        bx = bx+(-w(i))*dx./r3; by = by+(-w(i))*dy./r3; bz = bz+(-w(i))*dz./r3;
    end
    b = k_m * [bx; by; bz];
end

function draw_arrow(ax, ox, oy, oz, ux, uy, uz, col, hl, hw, lw)
    tx = ox+ux; ty = oy+uy; tz = oz+uz;
    mag = sqrt(ux^2+uy^2+uz^2);
    if mag < 1e-10, return; end
    dx = ux/mag; dy = uy/mag; dz = uz/mag;
    % Shaft
    sx = tx-hl*dx; sy = ty-hl*dy; sz = tz-hl*dz;
    line([ox sx],[oy sy],[oz sz], 'Color',col, 'LineWidth',lw, 'Parent',ax);
    % Triangular arrowhead
    if abs(dz) < 0.9
        perp = cross([dx dy dz],[0 0 1]);
    else
        perp = cross([dx dy dz],[1 0 0]);
    end
    perp = perp/norm(perp);
    p1 = [sx+hw*perp(1), sy+hw*perp(2), sz+hw*perp(3)];
    p2 = [sx-hw*perp(1), sy-hw*perp(2), sz-hw*perp(3)];
    patch([p1(1) p2(1) tx],[p1(2) p2(2) ty],[p1(3) p2(3) tz], ...
        col, 'EdgeColor',col, 'FaceAlpha',1, 'LineWidth',0.3, 'Parent',ax);
end
