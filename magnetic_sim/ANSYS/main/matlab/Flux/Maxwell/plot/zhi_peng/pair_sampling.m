%% pair_sampling.m -- Zhi-peng pair pole: calibration sampling, side view x-z at y = 0.
%  SAMP = 'xa'  (default): the ORIGINAL (2026-10-01) calibration sampling.
%    Reconstructed from the 10-01 hand-over (fitting_pair.m was deleted on 10-03):
%    origin = centre between the tips, Maxwell global (0, 0, 0.2887) mm (config/zhi_peng/R500 SPH_OFST)
%    samples on the actuator axis xa = direction from the centre to the P1 tip corner (magic angle,
%    35.26 deg below horizontal in the x-z plane), 2*Nr+1 = 11 points (Nr = 5), +-150 um, 30 um apart
%    R = 150 um calibration ball
%    Output: figures/zhi_peng/pair_pole/sampling.png
%  SAMP = 'tip': the NTU-style axis sampling used by utils/scripts/zhi_peng/zp_axis_fit.m (SAMP = 'tip')
%    origin = P1 tip apex at mid-tongue height (config/zhi_peng/Pair_pole WP), pole axis +x (into P1)
%    17 points, x = -20 .. -340 um every 20 um, z = 0; arrow = +l_hat direction (into P1)
%    Output: figures/zhi_peng/pair_pole/sampling_tip.png
%  Geometry (config/zhi_peng/Pair_pole, y = 0 section near the tips = tongues only):
%    P1 tongue x >= 408.0 um, z 0 .. 0.178 mm; P2 tongue x <= -408.3 um, z 0.402 .. 0.580 mm (AEDT)
if ~exist('SAMP','var'), SAMP = 'xa'; end
clearvars -except SAMP;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(here));   % .../Flux/Maxwell
FIG  = fullfile(CALROOT, 'figures', 'zhi_peng', 'pair_pole');

FS = 60; FSLAB = 46; LWBOX = 5; CANV = 14.5; LWD = 5;
CB = [0.05 0.10 0.95];  CR = [0.85 0.10 0.10];  CG = [0.45 0.45 0.45];

switch SAMP
    case 'xa'
        OUT = fullfile(FIG, 'sampling.png');
        ZC = 288.675;                                         % centre height [um] (global z)
        P1 = struct('x0',  408.0, 'z', [0, 178] - ZC);        % P1 tongue, extends to +x
        P2 = struct('x0', -408.3, 'z', [402, 580] - ZC);      % P2 tongue, extends to -x
        NR = 5;  R = 150;  MS = 14;
        u  = [sqrt(2/3), -1/sqrt(3)];                         % xa in the x-z plane (towards the P1 tip corner)
        s  = (-NR:NR) * (R/NR);                               % 11 points, 30 um apart
        xk = s * u(1);  zk = s * u(2);
        tk = [-300 0 300];  sx = tk(2) - tk(1);  xr = [tk(1)-sx, tk(end)+sx];
        tz = tk;  yr = xr;
    case 'tip'
        OUT = fullfile(FIG, 'sampling_tip.png');
        ZC = 89.0;                                            % origin height = P1 mid-tongue [um] (global z)
        P1 = struct('x0',    0.0, 'z', [0, 178] - ZC);        % P1 tongue, apex at the origin
        P2 = struct('x0', -816.3, 'z', [402, 580] - ZC);      % P2 tongue, 816.3 um to -x, 402 um higher
        DX = 20;  NP = 17;  MS = 8;
        xk = -(1:NP) * DX;  zk = zeros(size(xk));
        tk = [-800 -400 0];  sx = tk(2) - tk(1);  xr = [tk(1)-sx, tk(end)+sx];
        tz = [-400 0 400];  yr = [tz(1)-sx, tz(end)+sx];
    otherwise
        error('unknown SAMP ''%s''', SAMP);
end

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax = axes(fig, 'Position', [0.19 0.15 0.74 0.78]); hold(ax, 'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
fill(ax, [P1.x0 xr(2) xr(2) P1.x0], [P1.z(1) P1.z(1) P1.z(2) P1.z(2)], CG, 'EdgeColor','k', 'LineWidth', 2);
fill(ax, [xr(1) P2.x0 P2.x0 xr(1)], [P2.z(1) P2.z(1) P2.z(2) P2.z(2)], CG, 'EdgeColor','k', 'LineWidth', 2);
switch SAMP
    case 'xa'
        t = linspace(0, 2*pi, 400);
        plot(ax, R*cos(t), R*sin(t), '--', 'Color', CR, 'LineWidth', LWD);
        L = 1.5*xr(2);
        plot(ax, [-L L]*u(1), [-L L]*u(2), ':', 'Color', [0 0 0], 'LineWidth', 3);
        plot(ax, xk, zk, 'o', 'Color', CB, 'MarkerFaceColor', CB, 'MarkerSize', MS);
        plot(ax, 0, 0, 'o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'MarkerSize', 10);   % origin
    case 'tip'
        plot(ax, xr, [0 0], ':', 'Color', [0 0 0], 'LineWidth', 3);
        plot(ax, xk, zk, 'o', 'Color', CB, 'MarkerFaceColor', CB, 'MarkerSize', MS);
        % +l_hat definition: arrow from the tip (origin) toward +x, drawn above the P1 tongue
        ZA = 200;  XA = 300;  HL = 60;  HW = 26;  LWA = 5;
        plot(ax, [0 0], [P1.z(2) ZA+30], '-', 'Color', [0 0 0], 'LineWidth', 2);
        plot(ax, [0 XA-HL], [ZA ZA], '-', 'Color', [0 0 0], 'LineWidth', LWA);
        fill(ax, [XA-HL XA XA-HL], [ZA-HW ZA ZA+HW], [0 0 0], 'EdgeColor', 'none');
        text(ax, XA/2, ZA+12, '$\hat{\ell}$', 'Interpreter','latex', 'FontSize', FSLAB+14, ...
             'HorizontalAlignment','center', 'VerticalAlignment','bottom');
end
axis(ax, 'equal');  xlim(ax, xr);  ylim(ax, yr);
set(ax, 'XTick', tk, 'YTick', tz);
box(ax, 'on');
set(ax, 'FontSize', FS, 'FontWeight','bold', 'LineWidth', LWBOX, 'TickDir','in', 'Layer','top');
yoff = yr(1) - 0.035*diff(yr);
for xv = xr
    text(ax, xv, yoff, sprintf('%g', xv), 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
         'FontSize', FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{x\;(\mu m)}$', 'Interpreter','latex', 'FontSize', FSLAB);
ylabel(ax, '$\mathbf{z\;(\mu m)}$', 'Interpreter','latex', 'FontSize', FSLAB);
set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
print(fig, OUT, '-dpng', '-r200');
fprintf('saved %s\n', OUT);
