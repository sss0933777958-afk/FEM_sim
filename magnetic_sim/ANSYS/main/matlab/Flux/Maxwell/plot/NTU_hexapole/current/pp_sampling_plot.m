% pp_sampling_plot.m -- side view (x-z, y = 0) of the NTU pair_pole calibration samples.
% Frame: origin = right (P1) tip apex, pole axis +x, z from the plate mid-plane; left (P2) tip at x = -1000 um.
% Arrow: +l_hat direction (charge position l_hat on +x, inside the plate).
% Samples: SAMP = 'tip'  : 17 points on the axis, 20 um apart, starting 20 um in front of the tip (x = -20 .. -340 um)
%                  -> pair_pole/sampling.png
%          SAMP = 'c500' : 17 points on +-150 um around a centre 500 um in front of the tip (x = -650 .. -350 um),
%                  red dashed circle = the +-150 um range -> pair_pole/sampling_c500.png
if ~exist('SAMP','var'), SAMP = 'tip'; end
clearvars -except SAMP;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
DAT = fullfile(CALROOT, 'utils', 'data', 'NTU_hexapole');  FIG = fullfile(CALROOT, 'figures', 'NTU_hexapole');
NP = 17;  TH = 250;  XL = -1000;                    % sample count, plate thickness [um], XL = left (P2) tip [um]
switch SAMP
    case 'tip'
        OUT = fullfile(FIG, 'pair_pole', 'sampling.png');
        DX = 20;  xk = -(1:NP) * DX;                    % sample step [um]
    case 'c500'
        OUT = fullfile(FIG, 'pair_pole', 'sampling_c500.png');
        XC = -500;  R = 150;  xk = XC + linspace(-R, R, NP);   % centre, half range [um]
    otherwise
        error('unknown SAMP ''%s''', SAMP);
end
CR = [0.85 0.10 0.10];

FS = 60; FSLAB = 46; LWBOX = 5; CANV = 14.5;
tk = [-800 -400 0];  s = tk(2) - tk(1);  xr = [tk(1)-s, tk(end)+s];
tz = [-400 0 400];   yr = [tz(1)-s, tz(end)+s];
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax = axes(fig, 'Position', [0.19 0.15 0.74 0.78]); hold(ax, 'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
fill(ax, [0 xr(2) xr(2) 0], [-TH/2 -TH/2 TH/2 TH/2], [0.45 0.45 0.45], 'EdgeColor','k', 'LineWidth', 2);
fill(ax, [xr(1) XL XL xr(1)], [-TH/2 -TH/2 TH/2 TH/2], [0.45 0.45 0.45], 'EdgeColor','k', 'LineWidth', 2);
plot(ax, xr, [0 0], ':', 'Color', [0 0 0], 'LineWidth', 3);
if strcmp(SAMP, 'c500')
    t = linspace(0, 2*pi, 400);
    plot(ax, XC + R*cos(t), R*sin(t), '--', 'Color', CR, 'LineWidth', 5);
end
plot(ax, xk, zeros(size(xk)), 'o', 'Color', [0.05 0.10 0.95], 'MarkerFaceColor', [0.05 0.10 0.95], 'MarkerSize', 8);
% +l_hat definition: arrow from the tip (origin) toward +x, drawn above the plate
ZA = 200;  XA = 300;  HL = 60;  HW = 26;  LWA = 5;
plot(ax, [0 0], [TH/2 ZA+30], '-', 'Color', [0 0 0], 'LineWidth', 2);
plot(ax, [0 XA-HL], [ZA ZA], '-', 'Color', [0 0 0], 'LineWidth', LWA);
fill(ax, [XA-HL XA XA-HL], [ZA-HW ZA ZA+HW], [0 0 0], 'EdgeColor', 'none');
text(ax, XA/2, ZA+12, '$\hat{\ell}$', 'Interpreter','latex', 'FontSize', FSLAB+14, ...
     'HorizontalAlignment','center', 'VerticalAlignment','bottom');
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
