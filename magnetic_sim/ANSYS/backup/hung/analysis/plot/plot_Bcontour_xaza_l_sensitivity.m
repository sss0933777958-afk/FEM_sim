%% plot_Bcontour_xaza_l_sensitivity.m
%  |B| contour on the x_a-z_a plane (y_a = 0) for Coil1.
%  Produces TWO figures — baseline (l=500 µm) and variant (l=250 µm) —
%  in the exact style of the reference figure, for l-sensitivity comparison.
%
%  Data sources:
%    l=500: magnetic_sim/hung/results/coil1/filleted_conv/       (4-pass NREFINE, 732k nodes)
%    l=250: magnetic_sim/hung/results/coil1/filleted_l250_conv/  (4-pass NREFINE, 718k nodes)
%
%  Actuator frame (x_a, y_a, z_a) aligned with pole directions:
%    x_a -> P1 (lower, azimuth 0 deg)
%    y_a -> P3 (lower, azimuth 120 deg)
%    z_a -> P5 (upper, azimuth 60 deg)
%
%  Outputs:
%    magnetic_sim/hung/figures/coil1/Bcontour_xaza_Dfillet_l500.png
%    magnetic_sim/hung/figures/coil1/Bcontour_xaza_Dfillet_l250.png

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));

%% Config for the two runs
%  window_um = 0.6 * l keeps the 6 pole tips just outside the figure window
%  so each figure shows only the smooth air-domain field (style matches ref).
runs = struct( ...
    'l_um',      {500, 250}, ...
    'window_um', {300, 150}, ...
    'data_dir',  {fullfile('..','..','results','coil1','filleted_conv'), ...
                  fullfile('..','..','results','coil1','filleted_l250_conv')}, ...
    'out_name',  {'Bcontour_xaza_Dfillet_l500.png', ...
                  'Bcontour_xaza_Dfillet_l250.png'} ...
);

%% Actuator frame basis (ANSYS coords)
theta = 54.74 * pi/180;   % magic angle
s = sin(theta); c = cos(theta);

xa_hat = [ s*cosd(0);    s*sind(0);    -c ];   % P1 (lower)
ya_hat = [ s*cosd(120);  s*sind(120);  -c ];   % P3 (lower)
za_hat = [ s*cosd(60);   s*sind(60);   +c ];   % P5 (upper)

R = [xa_hat'; ya_hat'; za_hat'];

%% Grid resolution (built per-run since window differs)
Ngrid = 401;

%% Loop over two l values
summary = cell(length(runs), 4);

for k = 1:length(runs)
    l_um      = runs(k).l_um;
    win_um    = runs(k).window_um;
    data_dir  = runs(k).data_dir;
    out_name  = runs(k).out_name;

    fprintf('\n=== Processing l = %d um (window +/-%d um) ===\n', l_um, win_um);

    % Per-run interpolation grid
    win_m = win_um * 1e-6;
    xi = linspace(-win_m, win_m, Ngrid);
    zi = linspace(-win_m, win_m, Ngrid);
    [Xi, Zi] = meshgrid(xi, zi);

    % Load ANSYS nodes + B-field, rotate to actuator frame
    d = import_ansys_data(data_dir, 'all', 'coil1');
    P_act = R * [d.x, d.y, d.z]';
    x_a = P_act(1,:)';
    y_a = P_act(2,:)';
    z_a = P_act(3,:)';
    bsum = sqrt(d.bx.^2 + d.by.^2 + d.bz.^2);
    N_total = length(bsum);

    % Node-count report
    r_ansys  = sqrt(d.x.^2 + d.y.^2 + d.z.^2);
    N_win    = sum(r_ansys < win_m);
    N_cube50 = sum(abs(d.x)<50e-6 & abs(d.y)<50e-6 & abs(d.z)<50e-6);
    fprintf('  Total nodes (whole model)        : %d\n', N_total);
    fprintf('  R < %d um (figure visual core)   : %d\n', win_um, N_win);
    fprintf('  +/-50 um cube (fitting region)   : %d\n', N_cube50);

    % 3D scatteredInterpolant at y_a=0
    fprintf('  Building scatteredInterpolant...\n');
    tic;
    F3D = scatteredInterpolant(x_a, y_a, z_a, bsum, 'linear', 'nearest');
    fprintf('    built in %.1f s\n', toc);

    tic;
    Bi_raw = F3D(Xi, zeros(size(Xi)), Zi) * 1e3;   % mT
    fprintf('    queried in %.1f s\n', toc);

    B_WP = F3D(0, 0, 0) * 1e3;
    fprintf('  WP |B| = %.3f mT\n', B_WP);

    % Light Gaussian smoothing (sigma = 1 pixel ~ 1.5 um)
    % (base-MATLAB replacement for imgaussfilt — no Image Processing Toolbox needed)
    sigma  = 1.0;
    khalf  = ceil(3*sigma);
    xk     = -khalf:khalf;
    g1     = exp(-xk.^2 / (2*sigma^2));
    g1     = g1 / sum(g1);
    Bi     = conv2(g1, g1, Bi_raw, 'same');

    % Plot (match reference style)
    % Auto-scale per figure (matches reference style: smooth gradient end-to-end)
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

    text(-win_um*0.97, -win_um*0.93, sprintf('N_{nodes} = %s', addcommas(N_total)), ...
        'FontSize', 10, 'Color', [0.2 0.2 0.2], ...
        'BackgroundColor', [1 1 1 0.7], 'EdgeColor', 'none', 'Margin', 2);

    xlabel('x_a [\mum]', 'FontSize', 14);
    ylabel('z_a [\mum]', 'FontSize', 14);
    title(sprintf('|B| on x_a-z_a plane, l=%d\\mum', l_um), ...
        'FontSize', 15, 'FontWeight', 'bold', 'Interpreter','tex');
    axis equal;
    xlim([-win_um win_um]); ylim([-win_um win_um]);
    set(gca, 'FontSize', 12, 'GridAlpha', 0.3);
    grid on; box on;

    % Save
    out_path = fullfile('..','..','figures','coil1', out_name);
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, out_path, '-dpng', '-r150');
    fprintf('  Saved: %s\n', out_path);

    summary{k,1} = l_um;
    summary{k,2} = N_total;
    summary{k,3} = B_WP;
    summary{k,4} = out_name;

    close(fig);
end

%% Summary
fprintf('\n=== l-sensitivity summary (Coil1) ===\n');
fprintf('  l=%d um : WP |B| = %.3f mT , N_nodes = %d -> %s\n', ...
    summary{1,1}, summary{1,3}, summary{1,2}, summary{1,4});
fprintf('  l=%d um : WP |B| = %.3f mT , N_nodes = %d -> %s\n', ...
    summary{2,1}, summary{2,3}, summary{2,2}, summary{2,4});
ratio = summary{2,3} / summary{1,3};
fprintf('  Ratio l=250/l=500 |B| : %.3fx\n', ratio);

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
