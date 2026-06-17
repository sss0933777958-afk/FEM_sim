function plot_P1_topview_Bcontour()
% PLOT_P1_TOPVIEW_BCONTOUR  |B| contour on the P1 pole-local xa-ya plane
% for long2016_hexapole_halfcut, P1-only excitation.
%
% Reads FEM result `magnetic_sim/ANSYS/main/ANSYS_data/long2016_hexapole_halfcut/coil1/` re-run at
% I=1.0 A, TURNS=70, SmartSize 2 (5/17 04:44 finish). No scaling required.
%
% Frame:
%   alpha = atan(sqrt(2)) = 54.7356°  (Long Fei fixed)
%   P1 lower pole at azimuth 0°, tip in xz plane at +x −z
%   xa = sin(α)·X − cos(α)·Z   (along P1 axis, positive toward pole tip)
%   ya = Y                      (perpendicular, horizontal)
%   za = cos(α)·X + sin(α)·Z    (3rd axis, used only for slab thickness)
%
% Style mirrors magnetic_sim/ANSYS/main/matlab/kuo_quadrupole/fit/plot_P1_topview_Bcontour.m.

    %% ---- Style ----
    FONT_NAME = 'Helvetica';
    FONT_LBL  = 13;
    FONT_TTL  = 14;
    FONT_CB   = 12;
    DPI       = 300;

    %% ---- Paths ----
    kuo_root  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main';
    hung_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\hung';
    addpath(fullfile(hung_root, 'analysis', 'core'));

    res_dir = fullfile(kuo_root, 'results', 'long2016_hexapole_halfcut', 'coil1');
    out_dir = fullfile(kuo_root, 'figures', 'long2016_hexapole_halfcut');
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    %% ---- Constants ----
    I_TGT   = 1.0;              % FEM was solved at 1.0 A per coil — no scaling
    PROT_H  = 7e-3;             % protrusion height (matches sim script line 49)
    R_norm_z = 500e-6 / sqrt(3);% WP z-offset from yoke (matches sim script line 42)
    SPH_OFST = -PROT_H - 6e-3 + R_norm_z;  % WP global z (matches Long Fei convention)

    %% ---- Load FEM data (readlines + sscanf for both files) ----
    fprintf('Loading coil1 FEM data (CURR=%.2f A, SmartSize 2)...\n', I_TGT);
    coord_file  = fullfile(res_dir, 'coil1_coord_all.dat');
    bfield_file = fullfile(res_dir, 'coil1_bfield_all.dat');

    % --- coord: 4 cols per row (node x y z). Filter to data lines.
    fprintf('  reading %s ...\n', coord_file);
    L = readlines(coord_file);
    coord_data = zeros(numel(L), 4);
    nc = 0;
    for k = 1:numel(L)
        s = char(L(k));
        nums = sscanf(s, '%f');
        if numel(nums) == 4 && nums(1) == round(nums(1)) && nums(1) > 0
            nc = nc + 1;
            coord_data(nc, :) = nums.';
        end
    end
    coord = coord_data(1:nc, :);
    fprintf('  coord: %d rows\n', size(coord,1));

    % --- bfield: 5 cols per row (node BX BY BZ |B|). Handle concatenated negatives.
    fprintf('  reading %s ...\n', bfield_file);
    L = readlines(bfield_file);
    bfield_data = zeros(numel(L), 5);
    nb = 0;
    for k = 1:numel(L)
        s = char(L(k));
        s = regexprep(s, '(\d)([-+])', '$1 $2');
        nums = sscanf(s, '%f');
        if numel(nums) == 5 && nums(1) == round(nums(1)) && nums(1) > 0
            nb = nb + 1;
            bfield_data(nb, :) = nums.';
        end
    end
    bfield = bfield_data(1:nb, :);
    fprintf('  bfield: %d rows\n', size(bfield,1));

    % --- Merge on node ID
    [~, ic, ib] = intersect(coord(:,1), bfield(:,1));
    fprintf('  Matched %d nodes\n', numel(ic));
    d.x    = coord(ic, 2);
    d.y    = coord(ic, 3);
    d.z    = coord(ic, 4);
    d.bx   = bfield(ib, 2);
    d.by   = bfield(ib, 3);
    d.bz   = bfield(ib, 4);
    d.bsum = bfield(ib, 5);

    % Merge on node ID
    [node_ids, ic, ib] = intersect(coord(:,1), bfield(:,1));
    fprintf('  Matched %d nodes\n', length(node_ids));
    d.x    = coord(ic, 2);
    d.y    = coord(ic, 3);
    d.z    = coord(ic, 4);
    d.bx   = bfield(ib, 2);
    d.by   = bfield(ib, 3);
    d.bz   = bfield(ib, 4);
    d.bsum = bfield(ib, 5);

    % [FIXED] Align with Long Fei generate_figures_2_4.m: shift z so WP is at z=0,
    % then use GLOBAL xy directly (NOT pole-local rotation). The previous rotation
    % bug placed "WP center" at global yoke center (0,0,0), which is ~12.7 mm AWAY
    % from the real WP at (0, 0, SPH_OFST=-12.7 mm), giving 60× too-weak |B|.
    x_a =  d.x;                  % global x (P1 lower pole at azimuth 0° → +x)
    y_a =  d.y;                  % global y
    z_a =  d.z - SPH_OFST;      % WP-shifted z; z_a=0 is WP plane

    % |B| is rotation-invariant
    bsum = d.bsum;

    %% ---- Diagnostic: distance distribution of nodes from WP ----
    r_wp = sqrt(d.x.^2 + d.y.^2 + d.z.^2);
    fprintf('  Total nodes: %d\n', numel(r_wp));
    fprintf('  Nodes within 100 µm of WP:  %d\n', sum(r_wp < 100e-6));
    fprintf('  Nodes within 500 µm of WP:  %d\n', sum(r_wp < 500e-6));
    fprintf('  Nodes within  1 mm of WP:   %d\n', sum(r_wp < 1e-3));
    fprintf('  Nodes within  7 mm of WP:   %d\n', sum(r_wp < 7e-3));

    %% ---- Slice z ≈ WP (50 µm slab matching Long Fei) ----
    z_tol = 50e-6;                    % ±50 µm slab around WP plane (matches Long Fei)
    mask  = abs(z_a) < z_tol;
    fprintf('  Slice nodes (|z_a - WP| < %.0f µm): %d\n', z_tol*1e6, sum(mask));
    if sum(mask) < 50
        z_tol = 100e-6;
        mask = abs(z_a) < z_tol;
        fprintf('  Expanded to %.0f µm slab: %d nodes\n', z_tol*1e6, sum(mask));
    end

    xs = x_a(mask)*1e6;     % µm
    ys = y_a(mask)*1e6;
    bs = bsum(mask)*1e3;    % mT

    %% ---- Plot range: ±300 µm (match Long Fei Fig 2.4) ----
    xlim_v = [-300, 300];
    ylim_v = [-300, 300];

    in_range = xs >= xlim_v(1) & xs <= xlim_v(2) & ...
               ys >= ylim_v(1) & ys <= ylim_v(2);
    fprintf('  Nodes in plot range ±%d µm: %d\n', xlim_v(2), sum(in_range));
    xs = xs(in_range); ys = ys(in_range); bs = bs(in_range);

    %% ---- Interpolate to regular grid for contourf ----
    grid_n = 120;
    [Xg, Yg] = meshgrid(linspace(xlim_v(1), xlim_v(2), grid_n), ...
                         linspace(ylim_v(1), ylim_v(2), grid_n));
    F  = scatteredInterpolant(xs, ys, bs, 'linear', 'linear');
    Bg = F(Xg, Yg);

    %% ---- Plot ----
    fig = figure('Position', [60 60 900 760], 'Color', 'w');
    set(fig, 'DefaultAxesFontName', FONT_NAME);

    bmin = max(0.0, floor(min(Bg(:))*10)/10);
    bmax = ceil(max(Bg(:))*10)/10;
    levels = linspace(bmin, bmax, 25);

    contourf(Xg, Yg, Bg, levels, 'LineStyle', 'none');
    colormap(jet);
    hold on;

    % Center mark + |B| value annotation
    [~, ic] = min(abs(Xg(1,:)));
    [~, ir] = min(abs(Yg(:,1)));
    B_center = Bg(ir, ic);
    plot(0, 0, 'k+', 'MarkerSize', 14, 'LineWidth', 2);
    text(75, -37, sprintf('%.2f mT', B_center), ...
         'FontSize', 14, 'FontWeight', 'bold', ...
         'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 3);

    % N_nodes annotation
    text(xlim_v(1) + 45, ylim_v(1) + 45, ...
         sprintf('N_{nodes} = %s', regexprep(num2str(length(xs)), ...
                 '(\d)(?=(\d{3})+$)', '$1,')), ...
         'FontSize', 11, ...
         'BackgroundColor', [1 1 1 0.85], 'EdgeColor', [0.5 0.5 0.5]);

    cb = colorbar;
    ylabel(cb, '|B| [mT]', 'FontSize', FONT_CB);
    set(cb, 'FontSize', FONT_CB);

    axis equal;
    grid on;
    xlim(xlim_v); ylim(ylim_v);
    xlabel(['x [' char(956) 'm]   (P1 axis at azimuth 0\circ \rightarrow)'], 'FontSize', FONT_LBL);
    ylabel(['y [' char(956) 'm]'], 'FontSize', FONT_LBL);
    title(sprintf('|B|  on  x-y  plane  (z=WP),  P1 only, I = %.1f A', I_TGT), ...
          'FontSize', FONT_TTL);

    out_path = fullfile(out_dir, 'Bcontour_xaya_P1_1p0A.png');
    exportgraphics(fig, out_path, 'Resolution', DPI);
    fprintf('Saved: %s\n', out_path);
end
