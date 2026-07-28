function plot_charge_circuit_sideview()
% PLOT_CHARGE_CIRCUIT_SIDEVIEW
%   Magnetic-circuit arrow (quiver) side view on the y=0 slice for coil1 (P1
%   self-excitation), no_gap FEM, zoomed near the working point, with the fitted
%   P1 charge_location overlaid as a PINK dot. One figure per package:
%     fix    -> single-parameter charge  -> current_base/figures/charge_circuit_sideview.png
%     no_fix -> 18-param bias charge      -> current_base/figures/charge_circuit_sideview.png
%
%   User convention (2026-07-14):
%     * PLOT CANVAS = measure/APDL frame (window x:-2..2 mm, z:-14..-12 mm, WP at
%       z=SPH_OFST~=-12.711 mm -> NOT the origin, NO WP marker; upright poles).
%     * The quantity being shown (charge_location) is COMPUTED in the actuator frame
%       (the fits), then mapped back to measure coords for display:
%         P_wp = R_act.' * charge_act ;  x_meas = P_wp(1) ; z_meas = P_wp(3) + SPH_OFST
%
%   Real FEM nodes only (NO interpolation): bin the xz view into cells, keep the
%   highest-|B| REAL node per cell (same sanctioned approach as
%   field_viz/plot_Bcircuit_gap_compare.m).  Style = bold-framed (option 1).

    %% --- Paths (import_ansys_data / mt_constants live in backup analysis) ---
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');

    CAL = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
           'long2016_hexapole_halfcut\Calibration_using_FEM_modeling'];
    DATA_ROOT = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\' ...
                 'long2016_hexapole_halfcut\data\no_gap\coil1'];

    %% --- View + sampling params (measure/APDL frame, mm) ---
    XMIN = -2.0; XMAX =  2.0;          % mm
    ZMIN = -14.0; ZMAX = -12.0;        % mm
    Y_SLAB    = 0.3;                   % mm : |y| < Y_SLAB slab about y=0
    GRID_X    = 34;                    % cells across x (~0.12 mm)
    GRID_Z    = 17;                    % cells across z (~0.12 mm)
    THRESH_T  = 1e-4;                  % drop |B| < 0.1 mT nodes

    %% --- Constants + actuator rotation (measure/WP <-> actuator) ---
    cnst = mt_constants();
    SPH_OFST = cnst.SPH_OFST;                           % m (WP z in APDL frame, ~ -12.711e-3)
    tip   = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];   % 3x6 (WP frame)
    dhat  = tip ./ vecnorm(tip);
    R_act = [dhat(:,1), dhat(:,3), dhat(:,5)].';        % measure(WP) -> actuator rotation

    %% --- Load no_gap coil1 field ONCE (shared by both packages) ---
    fprintf('Loading no_gap coil1 field ...\n');
    d = import_ansys_data(DATA_ROOT, 'all', 'coil1');
    d.bsum = sqrt(d.bx.^2 + d.by.^2 + d.bz.^2);
    fprintf('  FINGERPRINT no_gap coil1 : matched %d nodes, |B|max = %.4f T\n', numel(d.x), max(d.bsum));

    %% --- Shared color CAP so both figures share the same scale ---
    CAP = prctile(d.bsum, 92);

    %% --- Per package: load charge_act_P1, map to measure, render ---
    PKGS = { 'fix',    fullfile(CAL,'current_base','data','field_err_hist_no_gap.mat'),    fullfile(CAL,'current_base','figures','eighteen_param','charge_circuit_sideview.png'),    'single-parameter';
             'no_fix', fullfile(CAL,'current_base','data','field_err_hist_no_gap.mat'), fullfile(CAL,'current_base','figures','eighteen_param','charge_circuit_sideview.png'), '18-param bias' };

    for p = 1:size(PKGS,1)
        S = load(PKGS{p,2});                            % has charge_act_P1 (3x1 actuator, m)
        Pw = R_act.' * S.charge_act_P1(:);              % actuator -> WP frame (m)
        cx_mm = Pw(1)*1e3;                              % measure x [mm]
        cz_mm = Pw(3)*1e3 + SPH_OFST*1e3;               % measure z [mm] (add back SPH_OFST)
        fprintf('[%s] P1 charge (measure): x=%.3f mm, z=%.3f mm\n', PKGS{p,1}, cx_mm, cz_mm);

        render_one(d, CAP, XMIN, XMAX, ZMIN, ZMAX, Y_SLAB, GRID_X, GRID_Z, THRESH_T, ...
                   cx_mm, cz_mm, PKGS{p,4}, PKGS{p,3}, cnst, SPH_OFST);
    end
end

% =========================================================================
function render_one(d, CAP, XMIN, XMAX, ZMIN, ZMAX, Y_SLAB, GRID_X, GRID_Z, THRESH_T, ...
                    cx_mm, cz_mm, title_tag, out_path, cnst, SPH_OFST)

    x_mm = d.x*1e3;  y_mm = d.y*1e3;  z_mm = d.z*1e3;

    % --- y=0 slab + view-box restriction ---
    sel = abs(y_mm) < Y_SLAB & ...
          x_mm >= XMIN & x_mm <= XMAX & z_mm >= ZMIN & z_mm <= ZMAX & ...
          d.bsum > THRESH_T;
    xs = x_mm(sel); zs = z_mm(sel);
    bxs = -d.bx(sel); bzs = -d.bz(sel); bss = d.bsum(sel);   % source convention: coil1=P1 (lower/sink) -> B=-B_FEM (arrows emanate from P1 tip)
    fprintf('[render] %-16s : %d nodes in slab/view\n', title_tag, numel(xs));

    % --- per-cell max-|B| REAL node (no interpolation) ---
    x_edges = linspace(XMIN, XMAX, GRID_X+1);
    z_edges = linspace(ZMIN, ZMAX, GRID_Z+1);
    ix = discretize(xs, x_edges);  iz = discretize(zs, z_edges);
    ok = ~isnan(ix) & ~isnan(iz);
    xs=xs(ok); zs=zs(ok); bxs=bxs(ok); bzs=bzs(ok); bss=bss(ok);
    ix=ix(ok); iz=iz(ok);
    cell_id = (ix-1)*GRID_Z + iz;
    [~, order] = sort(bss, 'descend');
    cell_sorted = cell_id(order);
    [~, ia] = unique(cell_sorted, 'stable');           % first = max-|B| per cell
    idx = order(ia);
    xq = xs(idx); zq = zs(idx);
    bxq = bxs(idx); bzq = bzs(idx); bmq = bss(idx);
    fprintf('           arrows (1 per cell) = %d\n', numel(xq));

    % --- arrow vectors: direction-emphasis, power-law length compression ---
    ARROW_MAX = 0.9*(XMAX-XMIN)/GRID_X;                 % mm : full-magnitude arrow ~= cell size
    bxz = sqrt(bxq.^2 + bzq.^2);  bxz(bxz==0) = 1e-12;
    Lnorm = min((bmq./CAP).^0.25, 1.2);                 % compress large dynamic range
    uq = (bxq./bxz) .* Lnorm * ARROW_MAX;
    wq = (bzq./bxz) .* Lnorm * ARROW_MAX;

    %% --- Display units (Unit Reference Sheet): length um, field mT ---
    UL = 1e3;   % mm -> um
    UF = 1e3;   % T  -> mT

    %% --- Figure ---
    fig = figure('Position',[60 80 900 620], 'Color','w');
    ax = axes(fig); hold(ax,'on');

    CAP_mT = CAP * UF;
    nb = 24;  edges_b = linspace(0, CAP_mT, nb+1);  cmap = turbo(nb);
    bmq_mT = bmq * UF;
    lw_rng = [0.5, 2.2];
    for k = 1:nb
        if k < nb
            in = bmq_mT >= edges_b(k) & bmq_mT < edges_b(k+1);
        else
            in = bmq_mT >= edges_b(k);
        end
        if any(in)
            lw = lw_rng(1) + (k-1)/(nb-1)*(lw_rng(2)-lw_rng(1));
            quiver(ax, xq(in)*UL, zq(in)*UL, uq(in)*UL, wq(in)*UL, 0, ...
                   'Color', cmap(k,:), 'LineWidth', lw, 'MaxHeadSize', 0.5);
        end
    end

    % --- P1 & P2 pole outlines (black; measure frame, um) ---
    draw_pole_outlines(ax, cnst, SPH_OFST, UL);

    % --- P1 charge_location : PINK dot (computed in actuator, shown in measure) ---
    plot(ax, cx_mm*UL, cz_mm*UL, 'o', 'MarkerSize',10, ...
         'MarkerFaceColor',[1 0.45 0.75], 'MarkerEdgeColor',[0.80 0.20 0.50], 'LineWidth',1.4);

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
    ax.XAxis.Exponent = 0;  ax.YAxis.Exponent = 0;     % no x10^4 offset (collides with title); show full um values
    xlabel(ax, 'x (\mum)', 'FontWeight','bold', 'Interpreter','tex');
    ylabel(ax, 'z (\mum)', 'FontWeight','bold', 'Interpreter','tex');
    % (title removed per user request)

    exportgraphics(fig, out_path, 'Resolution', 300);
    fprintf('           Saved: %s\n', out_path);
    close(fig);
end

% =========================================================================
function draw_pole_outlines(ax, cnst, SPH_OFST, UL)
% Draw P1 (lower, schematic wedge) + P2 (upper, full cone via pole_axis) outlines
% in the measure frame, display units = um (multiply mm by UL). Black, LineWidth 2.
% Geometry reused verbatim from plot_P1_circuit_charge.m / plot_P2_circuit_charge.m.
    SPHmm    = SPH_OFST*1e3;                                  % mm
    R_norm_xy= cnst.R_norm_xy*1e3;  R_norm_z = cnst.R_norm_z*1e3;
    POLE_R   = cnst.POLE_R*1e3;     POLE_LEN = cnst.POLE_CONE_LEN*1e3;

    % --- P1 (lower pole, angle 0deg): schematic outline (plot_P1_circuit_charge) ---
    Z_CONE = SPHmm - R_norm_z;
    x_tip  = R_norm_xy;  x_base = x_tip + POLE_LEN;  x_end = (42+53)/2 - 10;   % =37.5 mm
    plot(ax,[x_tip x_base]*UL,[Z_CONE Z_CONE-POLE_R]*UL,'k-','LineWidth',2);   % cone slant
    plot(ax,[x_base x_end]*UL,[Z_CONE-POLE_R Z_CONE-POLE_R]*UL,'k-','LineWidth',2); % bottom edge
    plot(ax,[x_end x_end]*UL,[Z_CONE-POLE_R Z_CONE]*UL,'k-','LineWidth',2);    % right edge
    plot(ax,[x_end x_tip]*UL,[Z_CONE Z_CONE]*UL,'k-','LineWidth',2);           % flat top

    % --- P2 (upper pole, angle 180deg): full cone via pole_axis (plot_P2_circuit_charge) ---
    Tx = cnst.pole_tip_x(2)*1e3;  Tz = cnst.pole_tip_z_wp(2)*1e3 + SPHmm;      % P2 tip (measure mm)
    ax2 = cnst.pole_axis(:,2);  ux = ax2(1); uz = ax2(3);  un = hypot(ux,uz);  ux=ux/un; uz=uz/un;
    ha = atan(POLE_R/POLE_LEN);  ca = cos(ha);  sa = sin(ha);
    g1 = [ux*ca-uz*sa,  ux*sa+uz*ca];                        % axis rotated +ha
    g2 = [ux*ca+uz*sa, -ux*sa+uz*ca];                        % axis rotated -ha
    slant = POLE_LEN/ca;
    B1 = [Tx Tz] + slant*g1;  B2 = [Tx Tz] + slant*g2;
    plot(ax,[Tx B1(1)]*UL,[Tz B1(2)]*UL,'k-','LineWidth',2);
    plot(ax,[Tx B2(1)]*UL,[Tz B2(2)]*UL,'k-','LineWidth',2);
    plot(ax,[B1(1) B2(1)]*UL,[B1(2) B2(2)]*UL,'k-','LineWidth',2);
end
