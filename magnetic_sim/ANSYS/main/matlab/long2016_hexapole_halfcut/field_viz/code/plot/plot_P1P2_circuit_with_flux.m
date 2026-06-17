function plot_P1P2_circuit_with_flux()
% PLOT_P1P2_CIRCUIT_WITH_FLUX
%   兩張獨立圖,各 2-panel:
%     Top:  pole geometry side view (in pole-local frame, s = pole-axis arc)
%     Bot:  Φ(s) curve (from existing _smrt4.mat)
%   x axes(s)兩 panel 對齊。
%
%   P1:half disc (D-shape, lower halfcut)
%   P2:full disc (upper full cone)

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\hexapole-long2016\analysis');
    cnst = mt_constants();

    %% 幾何常數
    R_norm_xy    = cnst.R_norm_xy * 1e3;        % 0.408 mm
    POLE_R       = cnst.POLE_R * 1e3;           % 3 mm
    POLE_TIP_R   = cnst.POLE_TIP_R * 1e3;       % 0.040 mm
    cone_len     = cnst.POLE_CONE_LEN * 1e3;    % 15 mm
    s_tip        = R_norm_xy;                    % 0.408
    s_base       = s_tip + cone_len;             % 15.408
    s_block_end  = s_base + 8;                   % cone end + ~8 mm block
    s_xlim       = [-2, 25];

    data_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\flux_profile';
    fig_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';

    %% ===============================================================
    %% Figure 1: P1 (lower, half-disc D-shape)
    %% ===============================================================
    p1 = load(fullfile(data_dir, 'P1_flux_profile_smrt4.mat'));
    plot_one_pole('P1', p1.x_query_mm, p1.Phi_uWb, ...
        s_tip, s_base, s_block_end, POLE_R, POLE_TIP_R, s_xlim, ...
        true, ...                              % halfcut = true (D-shape)
        fullfile(fig_dir, 'P1_circuit_and_flux.png'));

    %% ===============================================================
    %% Figure 2: P2 (upper, full disc, in pole-local frame)
    %% ===============================================================
    p2 = load(fullfile(data_dir, 'P2_flux_profile_smrt4.mat'));
    plot_one_pole('P2', p2.s_query_mm, p2.Phi_uWb, ...
        s_tip, s_base, s_block_end, POLE_R, POLE_TIP_R, s_xlim, ...
        false, ...                             % full cone (not halfcut)
        fullfile(fig_dir, 'P2_circuit_and_flux.png'));
end


function plot_one_pole(label, s_arr, Phi_arr, ...
                       s_tip, s_base, s_block_end, POLE_R, POLE_TIP_R, ...
                       s_xlim, is_halfcut, out_path)
    fig = figure('Position', [50 50 1300 750], 'Color', 'w');
    set(fig, 'DefaultAxesFontName', 'Helvetica', 'DefaultAxesFontSize', 11);

    %% ----- Panel A:磁路側視(pole-local frame) -----
    ax1 = subplot(2, 1, 1);
    hold on; box on;

    % Cone outline (in pole-local: s = axial, r = radial)
    % Tip at (s_tip, 0), cone walls expand to ±POLE_R at s_base
    % For halfcut (P1): show only lower half (r ≤ 0)
    % For full   (P2): show both halves
    if is_halfcut
        % Lower half cone outline
        s_outline = [s_tip,      s_tip,       s_base,      s_block_end, s_block_end, s_tip,        s_tip];
        r_outline = [0,         -POLE_TIP_R, -POLE_R,      -POLE_R,      0,           0,            0];
        fill_color = [0.75 0.75 0.78];
        fill(s_outline, r_outline, fill_color, 'EdgeColor', [0.25 0.25 0.30], 'LineWidth', 1.8);
        % 標示 milled flat
        plot([s_tip s_block_end], [0 0], '-', 'Color', [0.85 0.15 0.15], 'LineWidth', 1.5);
        text((s_tip+s_base)/2, 0.55, 'milled flat (z = Z_{cone})', ...
             'FontSize', 10, 'Color', [0.85 0.15 0.15], 'Interpreter', 'tex', ...
             'HorizontalAlignment', 'center');
    else
        % Full cone (both halves)
        s_outline = [s_tip,       s_tip,        s_base,       s_block_end,  s_block_end, s_base,       s_tip];
        r_outline = [POLE_TIP_R, -POLE_TIP_R,  -POLE_R,       -POLE_R,       POLE_R,      POLE_R,       POLE_TIP_R];
        fill_color = [0.75 0.78 0.85];
        fill(s_outline, r_outline, fill_color, 'EdgeColor', [0.10 0.40 0.85], 'LineWidth', 1.8);
    end

    % Block boundary (鋼到 yoke)
    plot([s_base s_base], [-POLE_R, POLE_R*0.5], 'k--', 'LineWidth', 1.0);

    % ===== Flux arrows inside the pole (cone end → tip → WP, 朝 -x) =====
    % Several arrows along the cone center line + one outside tip toward WP
    flux_arrow_color = [0.90 0.40 0.10];   % orange
    if is_halfcut
        arrow_y = -POLE_R/2;          % half-disc 中線
    else
        arrow_y = 0;                  % full cone 中心
    end
    arrow_xs = [22, 18, 13, 8, 3];
    for k = 1:length(arrow_xs)
        x_arr = arrow_xs(k);
        % R(x) at this position (for arrow length)
        if x_arr >= s_base
            len_arrow = 2.5;
        else
            len_arrow = 2.5;
        end
        quiver(x_arr, arrow_y, -len_arrow, 0, 0, ...
               'Color', flux_arrow_color, 'LineWidth', 2.8, ...
               'MaxHeadSize', 0.7, 'AutoScale', 'off');
    end
    % Arrow exiting tip toward WP (in air gap)
    quiver(s_tip, 0, -1.5, 0, 0, 'Color', flux_arrow_color, ...
           'LineWidth', 3.0, 'MaxHeadSize', 1.0, 'AutoScale', 'off');
    % Annotate flux direction
    text(11, POLE_R + 0.5, '\bf flux  \Phi  (yoke \rightarrow tip \rightarrow WP)', ...
         'Color', flux_arrow_color, 'FontSize', 12, ...
         'HorizontalAlignment', 'center', 'Interpreter', 'tex');

    % WP + tip 標籤
    plot(0, 0, 'k+', 'MarkerSize', 14, 'LineWidth', 2);
    text(0, 0.7, 'WP', 'FontSize', 12, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center');
    plot(s_tip, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    text(s_tip, -0.6, sprintf('%s tip\n(x=%.2f)', label, s_tip), ...
         'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center');

    % cone end annotation
    text(s_base, POLE_R + 0.4, sprintf('cone end\n(x=%.1f)', s_base), ...
         'FontSize', 10, 'Color', [0.10 0.40 0.85], 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center');

    % R(x) taper label
    text((s_tip + s_base)/2, -POLE_R - 0.7, ...
         sprintf('R(x) taper: %.2f mm \\rightarrow %.1f mm', POLE_TIP_R, POLE_R), ...
         'FontSize', 10, 'HorizontalAlignment', 'center', 'Interpreter', 'tex');

    % Title
    if is_halfcut
        ttl = sprintf('A. %s 磁路側視(下極 halfcut,D-shape)', label);
    else
        ttl = sprintf('A. %s 磁路側視(上極 full cone)', label);
    end
    title(ttl, 'FontSize', 13);

    xlim(s_xlim);
    ylim([-POLE_R - 1.2, POLE_R + 1.2]);
    ylabel('r [mm]', 'FontSize', 12);
    grid on;
    set(ax1, 'XTickLabel', []);   % 隱藏 x 標籤(共享軸)

    %% ----- Panel B:Φ(s) curve -----
    ax2 = subplot(2, 1, 2);
    hold on; box on;

    plot(s_arr, Phi_arr, 'b-', 'LineWidth', 2.2);

    % 標 WP / tip / cone end / block 對應位置(虛線)
    ymax = max(Phi_arr) * 1.10;
    ymin = -ymax * 0.05;
    ylim([ymin, ymax]);

    plot([0 0], [ymin ymax], 'k--', 'LineWidth', 1.0);
    plot([s_tip s_tip], [ymin ymax], 'r--', 'LineWidth', 1.0);
    plot([s_base s_base], [ymin ymax], '--', 'Color', [0.10 0.40 0.85], ...
         'LineWidth', 1.0);

    text(0,      ymax*0.95, 'WP',       'FontSize', 10, 'HorizontalAlignment', 'left');
    text(s_tip,  ymax*0.85, sprintf('%s tip', label), 'FontSize', 10, 'Color', 'r');
    text(s_base, ymax*0.85, 'cone end', 'FontSize', 10, 'Color', [0.10 0.40 0.85]);

    % Mark Φ_max
    [phi_max, k_max] = max(Phi_arr);
    plot(s_arr(k_max), phi_max, 'bo', 'MarkerSize', 9, 'MarkerFaceColor', 'b');
    text(s_arr(k_max) + 0.4, phi_max + 0.05, ...
         sprintf('\\Phi_{max} = %.2f \\muWb @ x = %.1f mm', phi_max, s_arr(k_max)), ...
         'FontSize', 11, 'FontWeight', 'bold', 'Color', 'b', 'Interpreter', 'tex');

    xlim(s_xlim);
    xlabel('x [mm]   (pole axis arc-length from WP)', 'FontSize', 12);
    ylabel('\Phi(x) [\muWb]', 'FontSize', 12, 'Interpreter', 'tex');
    title(sprintf('B. %s axial flux Φ(x) = ∫ B_{axial} dA  (smrt 4,toward-WP positive)', label), ...
          'FontSize', 13);
    grid on;

    %% Link x-axes
    linkaxes([ax1, ax2], 'x');

    %% Save
    exportgraphics(fig, out_path, 'Resolution', 200);
    fprintf('Saved: %s\n', out_path);
end
