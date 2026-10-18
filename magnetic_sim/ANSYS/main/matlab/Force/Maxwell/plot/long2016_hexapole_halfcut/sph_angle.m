% sph_angle.m -- schematic of the force direction r_hat in the actuator frame (x_a, y_a, z_a drawn upright):
% azimuth theta measured in the x_a-y_a plane from x_a toward y_a, elevation phi from that plane toward z_a,
% r_hat = [cos(phi)cos(theta), cos(phi)sin(theta), sin(phi)].  No data is loaded (drawing only).
% Output: figures/long2016_hexapole_halfcut/hexapole/sph_angle.png
clearvars;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(here));   % .../Force/Maxwell
OUT  = fullfile(CALROOT, 'figures', 'long2016_hexapole_halfcut', 'hexapole', 'sph_angle.png');

% ---- style (figure-style.md) ---------------------------------------------------------
FS = 60;  CANV = 14.5;  LWA = 5.0;  LWV = 7.0;  LWD = 3.0;  LWC = 4.0;
CB = [0.05 0.10 0.95];  CR = [0.85 0.10 0.10];  GR = [0.45 0.45 0.45];
TH = 50;  PH = 35;                              % example direction [deg]
LA = 1.0;  LV = 0.9;   RT = 0.35;  RP = 0.55;   % axis length, vector length, theta / phi arc radii
AZ = 60;   EL = 20;                             % camera ~80 deg off r_hat in azimuth -> r_hat seen side-on

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes(fig, 'Position', [0.02 0.02 0.96 0.96]);  hold(ax, 'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
axis(ax, 'off');  daspect(ax, [1 1 1]);  view(ax, AZ, EL);
lim = [-0.25 1.15];  xlim(ax, lim);  ylim(ax, lim);  zlim(ax, lim);
cv = campos(ax) - camtarget(ax);  cv = cv / norm(cv);              % viewing direction, for arrow heads

% axes
E = eye(3);  ALAB = {'x_a', 'y_a', 'z_a'};
for k = 1:3
    arrow3_(ax, [0 0 0], LA*E(k,:), cv, 'k', LWA, 0.07);
    p = (LA + 0.09) * E(k,:);
    text(ax, p(1), p(2), p(3), ['$\mathbf{' ALAB{k} '}$'], 'Interpreter','latex', 'FontSize',FS, ...
         'HorizontalAlignment','center', 'VerticalAlignment','middle');
end

% direction, its projection on the x_a-y_a plane and the drop line
r  = [cosd(PH)*cosd(TH), cosd(PH)*sind(TH), sind(PH)];
rp = [cosd(TH), sind(TH), 0];
plot3(ax, [0 LV*rp(1)]*cosd(PH), [0 LV*rp(2)]*cosd(PH), [0 0], '--', 'Color',GR, 'LineWidth',LWD);
plot3(ax, LV*r(1)*[1 1], LV*r(2)*[1 1], [0 LV*r(3)], '--', 'Color',GR, 'LineWidth',LWD);
arrow3_(ax, [0 0 0], LV*r, cv, CR, LWV, 0.08);
p = (LV + 0.08) * r;
text(ax, p(1), p(2), p(3), '$\mathbf{\hat{r}}$', 'Interpreter','latex', 'FontSize',FS, 'Color',CR, ...
     'HorizontalAlignment','left', 'VerticalAlignment','bottom');

% theta arc: in the x_a-y_a plane, from x_a to the projection
t  = linspace(0, TH, 60);
A  = RT * [cosd(t); sind(t); zeros(size(t))].';
arc3_(ax, A, cv, CB, LWC);
p  = 1.35*RT * [cosd(TH/2), sind(TH/2), 0];
text(ax, p(1), p(2), p(3), '$\mathbf{\theta}$', 'Interpreter','latex', 'FontSize',FS, 'Color',CB, ...
     'HorizontalAlignment','center', 'VerticalAlignment','middle');

% phi arc: in the vertical plane through r_hat, from the projection up to r_hat
u  = linspace(0, PH, 60);
A  = RP * [cosd(u)*cosd(TH); cosd(u)*sind(TH); sind(u)].';
arc3_(ax, A, cv, CB, LWC);
p  = 1.2*RP * [cosd(0.75*PH)*cosd(TH), cosd(0.75*PH)*sind(TH), sind(0.75*PH)];   % outside the arc, clear of y_a
text(ax, p(1), p(2), p(3), '$\mathbf{\phi}$', 'Interpreter','latex', 'FontSize',FS, 'Color',CB, ...
     'HorizontalAlignment','left', 'VerticalAlignment','middle');

set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
print(fig, OUT, '-dpng', '-r200');
close(fig);
fprintf('wrote %s\n', OUT);

% ---------------------------------------------------------------------------------------
function arrow3_(ax, p0, p1, cv, col, lw, hl)
% Line p0 -> p1 with a flat triangular head of length hl, facing the camera.
    d = (p1 - p0) / norm(p1 - p0);
    w = cross(d, cv);  if norm(w) < 1e-6, w = cross(d, [0 0 1]); end
    w = w / norm(w);
    b = p1 - hl*d;
    plot3(ax, [p0(1) b(1)], [p0(2) b(2)], [p0(3) b(3)], '-', 'Color',col, 'LineWidth',lw);
    T = [p1; b + 0.45*hl*w; b - 0.45*hl*w];
    patch(ax, T(:,1), T(:,2), T(:,3), col, 'EdgeColor',col, 'LineWidth',1);
end

function arc3_(ax, A, cv, col, lw)
% Arc through the points A (n x 3) with an arrow head at its end.
    hl = 0.06;  d = A(end,:) - A(end-3,:);  d = d / norm(d);
    n  = size(A,1);  k = find(vecnorm(A - A(end,:), 2, 2) > hl, 1, 'last');  k = min(max(k,2), n);
    plot3(ax, A(1:k,1), A(1:k,2), A(1:k,3), '-', 'Color',col, 'LineWidth',lw);
    w = cross(d, cv);  w = w / norm(w);
    b = A(end,:) - hl*d;
    T = [A(end,:); b + 0.45*hl*w; b - 0.45*hl*w];
    patch(ax, T(:,1), T(:,2), T(:,3), col, 'EdgeColor',col, 'LineWidth',1);
end
