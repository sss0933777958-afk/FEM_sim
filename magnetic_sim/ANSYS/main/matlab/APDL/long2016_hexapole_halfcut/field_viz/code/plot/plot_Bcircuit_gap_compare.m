function plot_Bcircuit_gap_compare()
% PLOT_BCIRCUIT_GAP_COMPARE
%   Magnetic-circuit arrow (quiver) side view on the y=0 slice, coil1 (P1
%   excitation), for the two base-gap variants gap_300um and gap400 (=400um).
%   Purpose: let the user SEE how the flux (arrows) routes through the iron /
%   across the base gap differently for the two gap sizes, and where the
%   self-pole (P1) and opposing-pole (P2) Hall sensors sit.
%
%   Real FEM nodes only (NO interpolation): per plot-real-nodes rule, bin the
%   xz view into cells and pick the highest-|B| REAL node per cell (same
%   sanctioned approach as plot_P1_topview_mT.m, adapted to an xz side view
%   with a thin |y|<Y_SLAB slab).
%
%   Shared color scale across the two variants (figure-style comparison rule):
%   Pass 1 loads both and computes a common CAP; Pass 2 renders both with the
%   same clim so the comparison is fair.
%
%   Style = bold-framed (option 1): large bold fonts, 4-side box, grid off,
%   halved ticks, units in (), turbo colormap, |B| (T) colorbar.

    %% --- Paths ---
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath(['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\' ...
             'Calibration_using_FEM_modeling\voltage_base\code\main_function']);

    DATA_ROOT = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\' ...
                 'long2016_hexapole_halfcut\data\coil1'];
    OUT_DIR   = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
                 'long2016_hexapole_halfcut\field_viz\figures'];
    if ~exist(OUT_DIR, 'dir'); mkdir(OUT_DIR); end

    cnst = mt_constants();
    SPH_OFST = cnst.SPH_OFST * 1e3;   % -12.711 mm (APDL frame, WP z)

    %% --- View + sampling params ---
    XMIN = -52.5; XMAX = 52.5;        % mm (APDL x)
    ZMIN = -16.5; ZMAX =  2.0;        % mm (APDL z)
    Y_SLAB    = 2.0;                  % mm : |y| < Y_SLAB slab about y=0
    GRID_X    = 92;                   % cells across x (~1.14 mm)
    GRID_Z    = 16;                   % cells across z (~1.16 mm)
    ARROW_MAX = 1.45;                 % mm : arrow length at full magnitude
    THRESH_T  = 1e-4;                 % drop |B| < 0.1 mT nodes

    variants = { 'gap_300um', 'gap400',    'Bcircuit_gap300um.png', 'gap 300 \mum';
                 'gap400',    'gap400',    'Bcircuit_gap400um.png', 'gap 400 \mum' };
    % (col1 = data folder, col3 = out png, col4 = title tag) — variant list

    %% --- Sensor centers (P1 self, P2 opposing) in APDL frame [mm] ---
    [sp, sn] = build_sensor_geometry(cnst);      % 3x6 WP frame [m]
    sp_apdl = sp * 1e3;  sp_apdl(3,:) = sp_apdl(3,:) + SPH_OFST;   % -> APDL mm
    % sn is a unit normal (dimensionless) — same in both frames

    %% ================= PASS 1: load both, compute shared CAP =================
    D = struct('name',{},'bsum',{});
    for v = 1:size(variants,1)
        res_dir = fullfile(DATA_ROOT, variants{v,1});
        fprintf('\n[Pass1] Loading %s ...\n', variants{v,1});
        d = import_ansys_data(res_dir, 'all', 'coil1');
        d.bsum = sqrt(d.bx.^2 + d.by.^2 + d.bz.^2);
        fprintf('  FINGERPRINT %-10s : matched %d nodes, |B|max = %.4f T\n', ...
                variants{v,1}, numel(d.x), max(d.bsum));
        D(v).name = variants{v,1};
        D(v).d    = d;
    end
    CAP = max( prctile(D(1).d.bsum, 92), prctile(D(2).d.bsum, 92) );
    fprintf('\nShared color CAP (max of 92nd pct) = %.4f T\n', CAP);

    %% ================= PASS 2: render each variant ==========================
    for v = 1:size(variants,1)
        d = D(v).d;
        render_one(d, cnst, sp_apdl, sn, SPH_OFST, CAP, ...
                   XMIN, XMAX, ZMIN, ZMAX, Y_SLAB, GRID_X, GRID_Z, ...
                   ARROW_MAX, THRESH_T, variants{v,4}, ...
                   fullfile(OUT_DIR, variants{v,3}));
    end
end

% =========================================================================
function render_one(d, cnst, sp_apdl, sn, SPH_OFST, CAP, ...
                    XMIN, XMAX, ZMIN, ZMAX, Y_SLAB, GRID_X, GRID_Z, ...
                    ARROW_MAX, THRESH_T, title_tag, out_path)

    x_mm = d.x*1e3;  y_mm = d.y*1e3;  z_mm = d.z*1e3;

    % --- y=0 slab + view-box restriction ---
    sel = abs(y_mm) < Y_SLAB & ...
          x_mm >= XMIN & x_mm <= XMAX & z_mm >= ZMIN & z_mm <= ZMAX & ...
          d.bsum > THRESH_T;
    xs = x_mm(sel); zs = z_mm(sel);
    bxs = d.bx(sel); bzs = d.bz(sel); bss = d.bsum(sel);
    fprintf('[render] %-12s : %d nodes in slab/view\n', title_tag, numel(xs));

    % --- per-cell max-|B| REAL node (no interpolation) ---
    x_edges = linspace(XMIN, XMAX, GRID_X+1);
    z_edges = linspace(ZMIN, ZMAX, GRID_Z+1);
    ix = discretize(xs, x_edges);  iz = discretize(zs, z_edges);
    ok = ~isnan(ix) & ~isnan(iz);
    xs=xs(ok); zs=zs(ok); bxs=bxs(ok); bzs=bzs(ok); bss=bss(ok);
    ix=ix(ok); iz=iz(ok);
    cell_id = (ix-1)*GRID_Z + iz;
    [bs_sorted, order] = sort(bss, 'descend');       %#ok<ASGLU>
    cell_sorted = cell_id(order);
    [~, ia] = unique(cell_sorted, 'stable');         % first = max-|B| per cell
    idx = order(ia);
    xq = xs(idx); zq = zs(idx);
    bxq = bxs(idx); bzq = bzs(idx); bmq = bss(idx);
    fprintf('           arrows (1 per cell) = %d\n', numel(xq));

    % --- arrow vectors: direction-emphasis, power-law length compression ---
    bxz = sqrt(bxq.^2 + bzq.^2);  bxz(bxz==0) = 1e-12;
    Lnorm = min((bmq./CAP).^0.25, 1.2);              % compress 1000:1 range
    uq = (bxq./bxz) .* Lnorm * ARROW_MAX;
    wq = (bzq./bxz) .* Lnorm * ARROW_MAX;

    %% --- Display units (Unit Reference Sheet): length um, field mT ---
    UL = 1e3;   % mm -> um  (length axes in um per unit rule)
    UF = 1e3;   % T  -> mT  (field / colorbar in mT per unit rule)

    %% --- Figure ---
    fig = figure('Position',[40 60 1650 560], 'Color','w');
    ax = axes(fig); hold(ax,'on');

    CAP_mT  = CAP * UF;
    nb = 24;  edges_b = linspace(0, CAP_mT, nb+1);  cmap = turbo(nb);
    bmq_mT = bmq * UF;
    lw_rng = [0.5, 2.2];
    for k = 1:nb
        if k < nb
            in = bmq_mT >= edges_b(k) & bmq_mT < edges_b(k+1);
        else
            in = bmq_mT >= edges_b(k);               % top bin catches >= CAP
        end
        if any(in)
            lw = lw_rng(1) + (k-1)/(nb-1)*(lw_rng(2)-lw_rng(1));
            quiver(ax, xq(in)*UL, zq(in)*UL, uq(in)*UL, wq(in)*UL, 0, ...
                   'Color', cmap(k,:), 'LineWidth', lw, 'MaxHeadSize', 0.5);
        end
    end

    % --- WP crosshair ---
    plot(ax, 0, SPH_OFST*UL, 'k+', 'MarkerSize', 16, 'LineWidth', 2.4);
    text(ax, 1.2*UL, (SPH_OFST+1.0)*UL, 'WP', 'FontSize', 15, 'FontWeight','bold');

    % --- Sensors: P1 (self, col1) + P2 (opposing, col2) ---
    draw_sensor(ax, sp_apdl(:,1)*UL, sn(:,1), 'P1', UL);
    draw_sensor(ax, sp_apdl(:,2)*UL, sn(:,2), 'P2', UL);

    hold(ax,'off');

    % --- colorbar + bold-framed (option 1) style ---
    colormap(ax, turbo);  clim(ax, [0 CAP_mT]);
    cb = colorbar(ax);
    cb.Label.String = '|B| (mT)';  cb.Label.FontWeight = 'bold';
    cb.FontSize = 16;  cb.FontWeight = 'bold';
    cb.Ticks = cb.Ticks(1:2:end);

    axis(ax,'equal');
    xlim(ax,[XMIN XMAX]*UL);  ylim(ax,[ZMIN ZMAX]*UL);
    set(ax, 'FontSize',16, 'FontWeight','bold', 'LineWidth',2, 'TickLength',[.012 .012]);
    box(ax,'on');  grid(ax,'off');
    xt = get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));
    yt = get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    xlabel(ax, 'x (\mum)', 'FontWeight','bold', 'Interpreter','tex');
    ylabel(ax, 'z (\mum)', 'FontWeight','bold', 'Interpreter','tex');
    title(ax, sprintf('Coil1 (P1) magnetic circuit, y=0 side view  |  %s', title_tag), ...
          'FontSize', 15, 'FontWeight','bold', 'Interpreter','tex');

    exportgraphics(fig, out_path, 'Resolution', 150);
    fprintf('           Saved: %s\n', out_path);
    close(fig);
end

% =========================================================================
function draw_sensor(ax, p, n, name, U)
% Draw a small sensor panel (square) + n+ normal arrow + label at (x,z).
% p already in display units; U = length display scale (mm->um) for sizing.
    px = p(1); pz = p(3);
    % panel oriented along the in-plane component of n+ (perpendicular drawn)
    L = 1.6*U;  W = 0.5*U;                   % panel size in display units
    % tangent (perp to n+ in xz): rotate n by 90 deg
    nx = n(1); nz = n(3);  nn = hypot(nx,nz);  if nn==0, nn=1; end
    nx=nx/nn; nz=nz/nn;
    tx = -nz; tz = nx;                       % tangent
    c = [px pz];
    P = [ c + ( L/2)*[tx tz] + ( W/2)*[nx nz];
          c + (-L/2)*[tx tz] + ( W/2)*[nx nz];
          c + (-L/2)*[tx tz] + (-W/2)*[nx nz];
          c + ( L/2)*[tx tz] + (-W/2)*[nx nz] ];
    fill(ax, P(:,1), P(:,2), [0.92 0.15 0.15], 'EdgeColor','k', 'LineWidth',1.4);
    % n+ arrow
    aL = 1.8*U;
    plot(ax, [px px+aL*nx], [pz pz+aL*nz], '-', 'Color',[0.92 0.15 0.15], 'LineWidth',2.4);
    plot(ax, px+aL*nx, pz+aL*nz, '^', 'MarkerSize',8, ...
         'MarkerFaceColor',[0.92 0.15 0.15], 'MarkerEdgeColor',[0.92 0.15 0.15]);
    text(ax, px+2.4*U*nx, pz+2.4*U*nz, name, 'FontSize',15, 'FontWeight','bold', ...
         'Color',[0.75 0.05 0.05], 'HorizontalAlignment','center');
end
