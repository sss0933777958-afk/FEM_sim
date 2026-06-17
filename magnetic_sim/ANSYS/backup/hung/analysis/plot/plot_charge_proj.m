%% plot_charge_proj.m — Lower/Upper charge projection diagram
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));
cnst = mt_constants();
alpha = cnst.alpha;
TILT_DN = 5.71;
TILT_UP = 35.0;
dhat_ang = 90 - alpha*180/pi;  % 35.26 deg from horizontal

fig = figure('Position', [50 50 1600 700]);

%% === Left: Lower (P1) ===
subplot(1,2,1);
hold on; axis equal; grid on;

% Parameters
R = 500;  % tip distance [um]
s_lower = 485;
theta_lower = 29.6;
l_lower = 922;

% Directions (in XZ plane, angles from +X horizontal)
dhat_deg = -dhat_ang;  % negative because lower goes to -Z
axis_deg = -TILT_DN;

dhat_vec = [cosd(dhat_deg), sind(dhat_deg)];
axis_vec = [cosd(axis_deg), sind(axis_deg)];

% Tip position
tip = R * dhat_vec;

% Charge on d_hat
charge = l_lower * dhat_vec;

% Physical charge along pole axis
phys_charge = tip + s_lower * axis_vec;

% Projection point (physical charge projected onto d_hat)
proj_dist = dot(phys_charge, dhat_vec);
proj_point = proj_dist * dhat_vec;

% Draw pole body (thick gray line from tip along axis)
pole_end = tip + 600 * axis_vec;
plot([tip(1) pole_end(1)], [tip(2) pole_end(2)], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 8);

% d_hat line (blue, from origin past charge)
dhat_end = 1100 * dhat_vec;
plot([0 dhat_end(1)], [0 dhat_end(2)], 'b-', 'LineWidth', 1.5);

% Origin
plot(0, 0, 'k+', 'MarkerSize', 18, 'LineWidth', 3);
text(15, 25, 'WP', 'FontSize', 13, 'FontWeight', 'bold');

% Tip
plot(tip(1), tip(2), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
text(tip(1)-60, tip(2)+30, 'tip', 'FontSize', 12, 'Color', 'r', 'FontWeight', 'bold');

% Physical charge (along pole axis)
plot(phys_charge(1), phys_charge(2), 'gs', 'MarkerSize', 14, 'MarkerFaceColor', 'g', 'LineWidth', 2);
text(phys_charge(1)+20, phys_charge(2)+25, sprintf('s = %d um', s_lower), 'FontSize', 11, 'Color', [0 0.6 0], 'FontWeight', 'bold');

% Charge on d_hat (model)
plot(charge(1), charge(2), 'b*', 'MarkerSize', 16, 'LineWidth', 2);
text(charge(1)-80, charge(2)-40, sprintf('l = %d um', l_lower), 'FontSize', 11, 'Color', 'b', 'FontWeight', 'bold');

% Projection dashed line (from physical charge to d_hat)
plot([phys_charge(1) proj_point(1)], [phys_charge(2) proj_point(2)], 'k--', 'LineWidth', 1.2);

% s arrow along pole axis
mid_s = tip + s_lower/2 * axis_vec;
text(mid_s(1)+15, mid_s(2)+20, sprintf('s = %d', s_lower), 'FontSize', 11, 'Color', [0 0.6 0]);

% Depth annotation: 500 along d_hat
mid_R = R/2 * dhat_vec;
text(mid_R(1)-20, mid_R(2)+30, 'R = 500', 'FontSize', 10, 'Color', [0.4 0.4 0.4]);

% Angle arc theta
arc_r = 200;
th1 = dhat_deg; th2 = axis_deg;
th_arc = linspace(min(th1,th2), max(th1,th2), 30);
plot(tip(1)+arc_r*cosd(th_arc), tip(2)+arc_r*sind(th_arc), 'k-', 'LineWidth', 1.5);
mid_th = mean([th1 th2]);
text(tip(1)+arc_r*1.1*cosd(mid_th), tip(2)+arc_r*1.1*sind(mid_th)-10, ...
    sprintf('\\theta = %.1f\\circ', theta_lower), 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');

% cos(theta) annotation
text(650, -450, sprintf('l = 500 + s \\times cos(%.1f\\circ)', theta_lower), 'FontSize', 12, 'Interpreter', 'tex');
text(650, -490, sprintf('  = 500 + %d \\times %.2f', s_lower, cosd(theta_lower)), 'FontSize', 12, 'Interpreter', 'tex');
text(650, -530, sprintf('  = %d \\mum', l_lower), 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');

% Labels
xlabel('x [\mum]', 'FontSize', 13);
ylabel('z [\mum]', 'FontSize', 13);
title('Lower (P1): \theta = 29.6\circ', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'tex');
xlim([-50 1100]); ylim([-600 100]);
set(gca, 'FontSize', 12);

% Legend
h1 = plot(nan,nan,'b-','LineWidth',2);
h2 = plot(nan,nan,'-','Color',[0.7 0.7 0.7],'LineWidth',5);
h3 = plot(nan,nan,'gs','MarkerSize',10,'MarkerFaceColor','g');
h4 = plot(nan,nan,'b*','MarkerSize',12,'LineWidth',2);
legend([h1 h2 h3 h4], 'd\_hat', 'pole body', 'physical depth s', 'model charge l', ...
    'Location', 'northwest', 'FontSize', 10);

%% === Right: Upper (P2) ===
subplot(1,2,2);
hold on; axis equal; grid on;

s_upper = 418;
theta_upper = 0.3;
l_upper = 918;

% Upper goes to -X, +Z
dhat_deg_u = 180 - dhat_ang;  % ~144.74 deg
axis_deg_u = 180 - TILT_UP;  % 145 deg

dhat_vec_u = [cosd(dhat_deg_u), sind(dhat_deg_u)];
axis_vec_u = [cosd(axis_deg_u), sind(axis_deg_u)];

tip_u = R * dhat_vec_u;
charge_u = l_upper * dhat_vec_u;
phys_charge_u = tip_u + s_upper * axis_vec_u;
proj_point_u = dot(phys_charge_u, dhat_vec_u) * dhat_vec_u;

% Pole body
pole_end_u = tip_u + 600 * axis_vec_u;
plot([tip_u(1) pole_end_u(1)], [tip_u(2) pole_end_u(2)], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 8);

% d_hat
dhat_end_u = 1100 * dhat_vec_u;
plot([0 dhat_end_u(1)], [0 dhat_end_u(2)], 'b-', 'LineWidth', 1.5);

% Origin
plot(0, 0, 'k+', 'MarkerSize', 18, 'LineWidth', 3);
text(15, -25, 'WP', 'FontSize', 13, 'FontWeight', 'bold');

% Tip
plot(tip_u(1), tip_u(2), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
text(tip_u(1)+30, tip_u(2)-30, 'tip', 'FontSize', 12, 'Color', 'r', 'FontWeight', 'bold');

% Physical charge
plot(phys_charge_u(1), phys_charge_u(2), 'gs', 'MarkerSize', 14, 'MarkerFaceColor', 'g', 'LineWidth', 2);

% Model charge
plot(charge_u(1), charge_u(2), 'b*', 'MarkerSize', 16, 'LineWidth', 2);
text(charge_u(1)+20, charge_u(2)-30, sprintf('l = %d um', l_upper), 'FontSize', 11, 'Color', 'b', 'FontWeight', 'bold');

% s label
mid_s_u = tip_u + s_upper/2 * axis_vec_u;
text(mid_s_u(1)+20, mid_s_u(2)-25, sprintf('s = %d', s_upper), 'FontSize', 11, 'Color', [0 0.6 0], 'FontWeight', 'bold');

% Angle arc (very small, just show label)
text(tip_u(1)+100, tip_u(2)+50, sprintf('\\theta = %.1f\\circ', theta_upper), ...
    'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');
text(tip_u(1)+100, tip_u(2)+20, '(nearly aligned)', 'FontSize', 10, 'Color', [0.5 0.5 0.5]);

% Formula
text(-1050, 550, sprintf('l = 500 + s \\times cos(%.1f\\circ)', theta_upper), 'FontSize', 12, 'Interpreter', 'tex');
text(-1050, 510, sprintf('  = 500 + %d \\times %.2f', s_upper, cosd(theta_upper)), 'FontSize', 12, 'Interpreter', 'tex');
text(-1050, 470, sprintf('  = %d \\mum', l_upper), 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');

xlabel('x [\mum]', 'FontSize', 13);
ylabel('z [\mum]', 'FontSize', 13);
title('Upper (P2): \theta = 0.3\circ', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'tex');
xlim([-1100 50]); ylim([-100 600]);
set(gca, 'FontSize', 12);

h1 = plot(nan,nan,'b-','LineWidth',2);
h2 = plot(nan,nan,'-','Color',[0.7 0.7 0.7],'LineWidth',5);
h3 = plot(nan,nan,'gs','MarkerSize',10,'MarkerFaceColor','g');
h4 = plot(nan,nan,'b*','MarkerSize',12,'LineWidth',2);
legend([h1 h2 h3 h4], 'd\_hat', 'pole body', 'physical depth s', 'model charge l', ...
    'Location', 'northwest', 'FontSize', 10);

set(fig, 'PaperPositionMode', 'auto');
print(fig, fullfile('..','..','figures','analytic','charge_projection_lower_upper.png'), '-dpng', '-r200');
fprintf('Saved.\n');
