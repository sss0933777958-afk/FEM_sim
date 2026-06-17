function plot_P1_circuit_arrows()
% PLOT_P1_CIRCUIT_ARROWS
%   Schematic side-view of P1 lower pole magnetic circuit (APDL xz-plane).
%   Shows steel structure (cone + horizontal arm + vertical post + yoke top piece),
%   SOURC36 coil ring (cross-section), halfcut plane, and flux arrows for:
%     RED  : main flux path (coil -> post -> arm -> cone -> tip -> air -> return)
%     BLUE : leakage flux (coil -> air -> back to coil, no use of tip)
%
%   Pure schematic — no FEM data, drawn from MT_Sim_P1.txt geometry.
%
%   Geometry (APDL frame, mm, lower pole P1 at azimuth 0°):
%     YOKE      : x in [42, 53],     z in [ 0,  2]   (radial disc)
%     POST      : x in [42.5, 52.5], z in [-7,  0]   (cylinder cross-section)
%     COIL ring : x in [-]39.5/52.5 +- 3,  z in [-7,  0]  (cross-section both sides)
%     ARM block : x in [0, 37.5],    z in [-13, -7]
%     CONE (halfcut, lower half): x in [0.408, 15.4], z in [-13-R(x), -13]
%     TIP       : x = 0.408,         z = -13 (cone axis)
%     WP        : x = 0,             z = -12.71

    %% --- Style ---
    FONT_NAME = 'Helvetica';
    FONT_LBL  = 13;
    FONT_TTL  = 15;
    FONT_ANNOT = 11;
    LW_STEEL  = 1.6;
    LW_HALFCUT = 1.0;
    LW_ARROW  = 2.2;
    DPI       = 300;

    COL_STEEL  = [0.78 0.78 0.82];   % gray
    COL_COIL   = [0.88 0.55 0.10];   % orange
    COL_AIR    = [1.00 1.00 1.00];
    COL_MAIN   = [0.85 0.10 0.10];   % red (main flux loop)
    COL_LEAK   = [0.10 0.30 0.85];   % blue (leakage flux)
    COL_YOKE_RET = [0.10 0.60 0.15]; % green (return through air-gap to other poles)
    COL_HALFCUT = [0.45 0.45 0.45];  % dashed cut plane

    %% --- Geometry (APDL frame, all in mm) ---
    YOKE_IN_R  = 42;
    YOKE_OUT_R = 53;
    YOKE_MID_R = (YOKE_IN_R + YOKE_OUT_R)/2;   % 47.5
    YOKE_H     = 2;
    PROT_R     = 5;
    PROT_H     = 7;
    POLE_R     = 3;
    POLE_TIP_R = 0.04;
    CONE_LEN   = 15;
    R_norm     = 0.5;
    R_norm_xy  = R_norm * sqrt(2/3);      % 0.408
    R_norm_z   = R_norm / sqrt(3);        % 0.289
    SPH_OFST   = -PROT_H - 6 + R_norm_z;  % -12.711
    Z_CONE     = -PROT_H - 6;             % -13 (cone axis z)

    COIL_IN_R   = 5;
    COIL_OUT_R  = 8;
    COIL_H      = PROT_H;                 % 7
    COIL_Z_CTR  = -PROT_H + COIL_H/2;     % -3.5

    % Pole base center (X0_LOW)
    X0_LOW = YOKE_MID_R;                  % 47.5
    POST_X1 = X0_LOW - PROT_R;            % 42.5
    POST_X2 = X0_LOW + PROT_R;            % 52.5

    % Horizontal arm
    ARM_X1 = 0;
    ARM_X2 = X0_LOW - 10;                 % 37.5
    ARM_Z1 = Z_CONE;                      % -13
    ARM_Z2 = -PROT_H;                     % -7

    % Cone (halfcut, lower half)
    x_tip    = R_norm_xy;                 % 0.408
    x_base   = x_tip + CONE_LEN;          % 15.408

    %% --- Set up figure ---
    fig = figure('Position', [50 50 1500 750], 'Color', 'w');
    set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LBL);
    ax = axes('Position', [0.04 0.10 0.94 0.83]);
    hold(ax, 'on');

    %% --- Draw STEEL pieces (light gray patches with black outline) ---
    % YOKE (top piece, schematic — show only portion near P1)
    yoke_x = [POST_X2-1, YOKE_OUT_R+1, YOKE_OUT_R+1, POST_X2-1];
    yoke_z = [0, 0, YOKE_H, YOKE_H];
    patch(ax, yoke_x, yoke_z, COL_STEEL, 'EdgeColor', 'k', 'LineWidth', LW_STEEL);
    text(ax, (POST_X2+YOKE_OUT_R)/2-2, YOKE_H/2, 'Yoke', ...
        'FontSize', FONT_ANNOT, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

    % Yoke (showing it continues to the LEFT toward other poles abstractly)
    yoke_left_x = [-3, POST_X1+1, POST_X1+1, -3];
    yoke_left_z = [0, 0, YOKE_H, YOKE_H];
    patch(ax, yoke_left_x, yoke_left_z, COL_STEEL, 'EdgeColor', 'k', ...
        'LineWidth', LW_STEEL, 'LineStyle', '--');
    text(ax, 10, YOKE_H/2+0.4, '(yoke continues -> adjacent poles)', ...
        'FontSize', FONT_ANNOT-2, 'HorizontalAlignment', 'center', ...
        'Color', [0.4 0.4 0.4], 'FontAngle', 'italic');

    % POST (vertical cylinder cross-section at y=0)
    post_x = [POST_X1, POST_X2, POST_X2, POST_X1];
    post_z = [ARM_Z2, ARM_Z2, 0, 0];
    patch(ax, post_x, post_z, COL_STEEL, 'EdgeColor', 'k', 'LineWidth', LW_STEEL);
    text(ax, X0_LOW, -3.5, 'Post', 'FontSize', FONT_ANNOT, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');

    % HORIZONTAL ARM (block at x=[0, 37.5], z=[-13, -7])
    arm_x = [ARM_X1, ARM_X2, ARM_X2, ARM_X1];
    arm_z = [ARM_Z1, ARM_Z1, ARM_Z2, ARM_Z2];
    patch(ax, arm_x, arm_z, COL_STEEL, 'EdgeColor', 'k', 'LineWidth', LW_STEEL);
    text(ax, 20, -10, 'Horizontal arm (block)', ...
        'FontSize', FONT_ANNOT, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

    % CONE (halfcut, lower half — triangle below z=-13)
    % Side view: top edge from (x_tip, Z_CONE) to (x_base, Z_CONE),
    % bottom edge from (x_tip, Z_CONE-POLE_TIP_R) to (x_base, Z_CONE-POLE_R)
    cone_x = [x_tip, x_base, x_base, x_tip];
    cone_z = [Z_CONE, Z_CONE, Z_CONE-POLE_R, Z_CONE-POLE_TIP_R];
    patch(ax, cone_x, cone_z, COL_STEEL, 'EdgeColor', 'k', 'LineWidth', LW_STEEL);
    text(ax, 8, -14.6, 'Cone (halfcut)', ...
        'FontSize', FONT_ANNOT, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

    %% --- Draw COIL ring cross-section (2 rectangles, left + right of post) ---
    coil_z = [-PROT_H, 0, 0, -PROT_H];
    % Left coil cross-section (between post and -x side)
    coil_L_x = [X0_LOW-COIL_OUT_R, X0_LOW-COIL_IN_R, X0_LOW-COIL_IN_R, X0_LOW-COIL_OUT_R];
    patch(ax, coil_L_x, coil_z, COL_COIL, 'EdgeColor', 'k', 'LineWidth', 1.2);
    text(ax, X0_LOW-(COIL_IN_R+COIL_OUT_R)/2, -3.5, '⊙', ...
        'FontSize', 22, 'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold');
    % Right coil cross-section
    coil_R_x = [X0_LOW+COIL_IN_R, X0_LOW+COIL_OUT_R, X0_LOW+COIL_OUT_R, X0_LOW+COIL_IN_R];
    patch(ax, coil_R_x, coil_z, COL_COIL, 'EdgeColor', 'k', 'LineWidth', 1.2);
    text(ax, X0_LOW+(COIL_IN_R+COIL_OUT_R)/2, -3.5, '⊗', ...
        'FontSize', 22, 'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold');
    text(ax, X0_LOW, -8.5, 'SOURC36 coil ring', ...
        'FontSize', FONT_ANNOT, 'HorizontalAlignment', 'center', ...
        'Color', COL_COIL*0.7, 'FontWeight', 'bold');

    %% --- Draw HALFCUT plane (dashed horizontal line at z = -13) ---
    plot(ax, [-3 60], [Z_CONE Z_CONE], '--', 'Color', COL_HALFCUT, ...
        'LineWidth', LW_HALFCUT);
    text(ax, 58, Z_CONE+0.3, 'halfcut plane (z=-13)', ...
        'FontSize', FONT_ANNOT-1, 'Color', COL_HALFCUT, ...
        'HorizontalAlignment', 'right', 'FontAngle', 'italic');

    %% --- WP marker ---
    plot(ax, 0, SPH_OFST, 'k+', 'MarkerSize', 14, 'LineWidth', 2);
    text(ax, -0.5, SPH_OFST+0.7, 'WP', 'FontSize', FONT_ANNOT+1, ...
        'HorizontalAlignment', 'right', 'FontWeight', 'bold');

    % Tip marker
    plot(ax, x_tip, Z_CONE, 'r.', 'MarkerSize', 14);
    text(ax, x_tip+0.5, Z_CONE-0.8, 'P1 tip', 'FontSize', FONT_ANNOT, ...
        'Color', [0.7 0 0], 'FontWeight', 'bold');

    %% --- Draw FLUX ARROWS ---

    % === MAIN PATH (RED): coil -> post -> arm -> cone -> tip -> air -> [abstract return]
    % Inside steel (post going down, arm going left, cone going to tip)
    draw_arrow(ax, [X0_LOW, -1], [X0_LOW, -6.3], COL_MAIN, LW_ARROW);   % post down
    draw_arrow(ax, [X0_LOW-1, ARM_Z1/2+ARM_Z2/2], [POST_X1-1, ARM_Z1/2+ARM_Z2/2-0.5], COL_MAIN, LW_ARROW);  % post->arm corner
    draw_arrow(ax, [30, -10], [16, -10], COL_MAIN, LW_ARROW);    % arm leftward
    draw_arrow(ax, [12, -13.5], [3, -12.95], COL_MAIN, LW_ARROW);    % cone toward tip
    draw_arrow(ax, [x_tip+0.3, Z_CONE+0.3], [-1.5, Z_CONE+1.6], COL_MAIN, LW_ARROW);  % tip -> air toward WP

    % Air-gap to "other poles" (abstract loop in air going up-and-around-back)
    % Curved path from WP area going up to yoke region
    t = linspace(0, 1, 50);
    arc_x = -2 + 12*t.^2;
    arc_z = SPH_OFST + (Z_CONE - SPH_OFST) + 14*t;   % rises from -12.71 toward yoke
    arc_x_clip = arc_x(arc_z <= -0.5);
    arc_z_clip = arc_z(arc_z <= -0.5);
    plot(ax, arc_x_clip, arc_z_clip, '-', 'Color', COL_YOKE_RET, 'LineWidth', LW_ARROW);
    draw_arrow(ax, [arc_x_clip(end-3), arc_z_clip(end-3)], ...
        [arc_x_clip(end), arc_z_clip(end)], COL_YOKE_RET, LW_ARROW);
    text(ax, 4, -7, sprintf('Air gap return\n(via adjacent poles)'), ...
        'FontSize', FONT_ANNOT-1, 'Color', COL_YOKE_RET, 'FontWeight', 'bold');

    % Continue green path through yoke back to coil (going right along yoke)
    draw_arrow(ax, [POST_X1-3, YOKE_H/2], [POST_X1-0.5, YOKE_H/2], COL_YOKE_RET, LW_ARROW);

    % === LEAKAGE PATH (BLUE): coil -> radiates outward into air -> short-circuits back
    % Sideways leakage out the arm into surrounding air
    draw_arrow(ax, [25, -7.2], [22, -4], COL_LEAK, LW_ARROW);   % flux leak from top of arm
    draw_arrow(ax, [25, -12.8], [22, -16], COL_LEAK, LW_ARROW); % flux leak from bottom of arm (halfcut side)
    draw_arrow(ax, [30, -7.2], [27, -4], COL_LEAK, LW_ARROW);
    draw_arrow(ax, [35, -7.2], [32, -4], COL_LEAK, LW_ARROW);

    % Curved leakage loop: small loop near coil (illustrative, blue arc going from
    % near coil cross-section out into air and back)
    cx = X0_LOW - COIL_OUT_R - 2; cz = -3.5;
    th = linspace(-pi/2, 3*pi/2, 40);
    r_loop = 1.8;
    plot(ax, cx + r_loop*cos(th), cz + r_loop*sin(th), '-', ...
        'Color', COL_LEAK, 'LineWidth', LW_ARROW);
    % Small arrowhead on the loop
    draw_arrow(ax, [cx-r_loop*0.95, cz+r_loop*0.3], ...
        [cx-r_loop*0.85, cz+r_loop*0.7], COL_LEAK, LW_ARROW);

    text(ax, 26, -2.0, sprintf('Leakage flux\n(short-circuit\nin air)'), ...
        'FontSize', FONT_ANNOT-1, 'Color', COL_LEAK, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');

    %% --- Axis settings ---
    axis(ax, 'equal');
    xlim(ax, [-4 60]);
    ylim(ax, [-17 4]);
    xlabel(ax, 'x  [mm]  (APDL frame, lower pole P1 along +x from WP)', ...
        'FontSize', FONT_LBL, 'Interpreter', 'tex');
    ylabel(ax, 'z  [mm]', 'FontSize', FONT_LBL);
    title(ax, sprintf(['P1 lower pole magnetic circuit — side view (y=0, APDL frame)\n' ...
        '\\color[rgb]{0.85,0.10,0.10}━ main flux (steel)  ' ...
        '\\color[rgb]{0.10,0.60,0.15}━ air-gap return  ' ...
        '\\color[rgb]{0.10,0.30,0.85}━ leakage flux']), ...
        'FontSize', FONT_TTL, 'Interpreter', 'tex');
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'GridAlpha', 0.18, 'Layer', 'top');

    %% --- Save ---
    out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    out_path = fullfile(out_dir, 'P1_magnetic_circuit_schematic.png');
    exportgraphics(fig, out_path, 'Resolution', DPI);
    fprintf('Saved: %s\n', out_path);
end


function draw_arrow(ax, p0, p1, color, lw)
% Simple data-space arrow using quiver with no scaling.
    quiver(ax, p0(1), p0(2), p1(1)-p0(1), p1(2)-p0(2), 0, ...
        'Color', color, 'LineWidth', lw, 'MaxHeadSize', 1.5, ...
        'AutoScale', 'off');
end
