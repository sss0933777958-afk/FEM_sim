function plot_P1_topview_mT()
% PLOT_P1_TOPVIEW_MT
%   Top view (xy plane) B-field vector distribution — coil1 (P1 active)
%   from smrt 4 baseline halfcut, replicating hexapole-long2016 fig 2.3(a)
%   style.
%
%   Key approach (per hexapole-long2016 generate_figures_2_3.m):
%     - NOT a fixed-z slice.  Instead, MAX-PROJECTION through z:
%       split xy domain into 80x80 cells, pick node with HIGHEST |B| per cell.
%     - This automatically catches iron-channeled flux regardless of z.
%     - Power-law scaling B^0.25 for arrow length (compress 1000:1 range).
%     - 28 magnitude bins for color/linewidth (turbo colormap).
%     - Colorbar 0-1 Tesla (matches hexapole-long2016 range).

    %% --- Style ---
    FONT_NAME   = 'Helvetica';
    FONT_LABEL  = 12;
    FONT_TITLE  = 13;
    FONT_CB     = 11;
    FONT_ANNOT  = 10;
    LINE_MAIN   = 1.5;
    LINE_THIN   = 0.8;
    DPI         = 300;

    %% --- Paths ---
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');

    res_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data\coil1\standard';
    out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    %% --- Load coil1 (P1 active) RAW FEM (full dataset, no slicing) ---
    fprintf('Loading coil1 smrt 4 (P1 active) full dataset...\n');
    d = import_ansys_data(res_dir, 'all', 'coil1');
    d.bsum = sqrt(d.bx.^2 + d.by.^2 + d.bz.^2);
    fprintf('  Nodes: %d\n', length(d.x));
    fprintf('  |B| range: [%.3g, %.3g] T\n', min(d.bsum), max(d.bsum));

    %% --- Geometry (mm) ---
    YOKE_IN_R  = 42;
    YOKE_OUT_R = 53;
    YOKE_MID_R = (YOKE_IN_R + YOKE_OUT_R)/2;
    PROT_R     = 5;

    %% --- Spatial binning: 80x80 grid, pick max-|B| node per cell ---
    grid_cells = 80;
    x_mm = d.x * 1e3;
    y_mm = d.y * 1e3;
    x_edges = linspace(-80, 80, grid_cells+1);
    y_edges = linspace(-80, 80, grid_cells+1);

    ix_bin = discretize(x_mm, x_edges);
    iy_bin = discretize(y_mm, y_edges);
    valid_bin = ~isnan(ix_bin) & ~isnan(iy_bin);
    cell_id = zeros(length(x_mm), 1);
    cell_id(valid_bin) = (ix_bin(valid_bin)-1)*grid_cells + iy_bin(valid_bin);

    idx = zeros(grid_cells*grid_cells, 1);
    n_found = 0;
    unique_cells = unique(cell_id(cell_id > 0));
    for k = 1:length(unique_cells)
        cell_nodes = find(cell_id == unique_cells(k));
        [~, best] = max(d.bsum(cell_nodes));
        n_found = n_found + 1;
        idx(n_found) = cell_nodes(best);
    end
    idx = idx(1:n_found)';
    fprintf('  Max-projection arrows (1 per cell): %d\n', n_found);

    % Threshold: keep nodes with |B| > 0.1 mT
    valid = d.bsum(idx) > 1e-4;
    idx = idx(valid);
    fprintf('  After |B| > 0.1 mT threshold: %d\n', length(idx));

    x_q = d.x(idx)*1e3;
    y_q = d.y(idx)*1e3;
    bmag_q = d.bsum(idx);

    % Power-law arrow length B^0.25 (compresses 1000:1 -> 6:1 dynamic range)
    arrow_max = 4.0;   % mm (at |B|=1T)
    bxy = sqrt(d.bx(idx).^2 + d.by(idx).^2);
    bxy(bxy == 0) = 1e-10;
    scale_factor = arrow_max * bmag_q.^0.25 ./ bxy;
    bx_q = d.bx(idx) .* scale_factor;
    by_q = d.by(idx) .* scale_factor;

    %% --- Plot ---
    fig = figure('Position', [50 50 850 750], 'Color', 'w');
    set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LABEL);
    hold on;

    % 28 magnitude bins with varying linewidth (matches reference style)
    n_bins = 28;
    bmax = max(bmag_q);
    edges = linspace(0, bmax, n_bins + 1);
    cmap_a = turbo(n_bins);
    lw_range = [0.40, 1.75];

    for k = 1:n_bins
        in_bin = bmag_q >= edges(k) & bmag_q < edges(k+1);
        if k == n_bins
            in_bin = in_bin | (bmag_q >= edges(end));
        end
        if any(in_bin)
            lw = lw_range(1) + (k-1)/(n_bins-1) * (lw_range(2) - lw_range(1));
            quiver(x_q(in_bin), y_q(in_bin), bx_q(in_bin), by_q(in_bin), 0, ...
                'Color', cmap_a(k,:), 'LineWidth', lw, 'MaxHeadSize', 0.3);
        end
    end

    % Yoke inner/outer rings (red dashed)
    theta_c = linspace(0, 2*pi, 200);
    plot(YOKE_IN_R *cos(theta_c), YOKE_IN_R *sin(theta_c), 'r--', 'LineWidth', LINE_THIN);
    plot(YOKE_OUT_R*cos(theta_c), YOKE_OUT_R*sin(theta_c), 'r--', 'LineWidth', LINE_THIN);

    % 6 post circles labeled P1-P6
    pole_data = {'P1', 0; 'P3', 120; 'P6', 240; 'P5', 60; 'P2', 180; 'P4', 300};
    for k = 1:size(pole_data,1)
        ang = pole_data{k,2} * pi/180;
        cx = YOKE_MID_R * cos(ang);
        cy = YOKE_MID_R * sin(ang);
        plot(cx + PROT_R*cos(theta_c), cy + PROT_R*sin(theta_c), ...
            'k-', 'LineWidth', LINE_MAIN);
        text(cx*1.15, cy*1.15, pole_data{k,1}, ...
            'HorizontalAlignment', 'center', 'FontSize', FONT_ANNOT, ...
            'FontWeight', 'bold', 'FontName', FONT_NAME);
    end

    % WP center crosshair
    plot(0, 0, 'k+', 'MarkerSize', 14, 'LineWidth', 2);
    hold off;

    % Colorbar 0-1 Tesla (matches hexapole-long2016 fig 2.3a)
    colormap(turbo);
    clim([0, 1.0]);
    cb = colorbar;
    ylabel(cb, 'Tesla', 'FontSize', FONT_CB);

    axis equal; grid on;
    set(gca, 'Color', 'w');
    xlim([-80 80]); ylim([-80 80]);
    xlabel('x [mm]', 'FontSize', FONT_LABEL);
    ylabel('y [mm]', 'FontSize', FONT_LABEL);
    title('Top view: B-field vector distribution (Coil1 = P1, unit: Tesla)', ...
        'FontSize', FONT_TITLE);

    %% --- Save ---
    out_path = fullfile(out_dir, 'P1_Bvector_topview_mT.png');
    exportgraphics(fig, out_path, 'Resolution', DPI);
    fprintf('\nSaved: %s\n', out_path);
end
