% sp_cost_b_plot.m -- profile cost J(b) (a re-fitted at every b), b = 0..500 um, minimum marked
% (loads sp_cost_b.mat only; computed in sp_cost_b.m).
if ~exist('CASE','var'), CASE = 'single_pole'; end   % 'single_pole' | 'pair_pole' (NTU) | 'zhi_peng' (pair, P1)
if ~exist('SAMP','var'), SAMP = 'tip'; end   % 'tip' | 'c500' (NTU only)
clearvars -except CASE SAMP;  close all;
SFX = '';  if ~strcmp(SAMP, 'tip'), SFX = ['_' SAMP]; end
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
% case -> {data folder model, .mat, figure folder}
CS  = struct('single_pole', {{'NTU_hexapole', 'sp_cost_b.mat', fullfile('NTU_hexapole', 'single_pole')}}, ...
             'pair_pole',   {{'NTU_hexapole', 'pp_cost_b.mat', fullfile('NTU_hexapole', 'pair_pole')}}, ...
             'zhi_peng',    {{'zhi_peng',     'zp_cost_b.mat', fullfile('zhi_peng', 'pair_pole')}});
CS  = CS.(CASE);
S   = load(fullfile(CALROOT, 'utils', 'data', CS{1}, strrep(CS{2}, '.mat', [SFX '.mat'])));
OUT = fullfile(CALROOT, 'figures', CS{3}, ['cost_b' SFX '.png']);
if ~isfolder(fileparts(OUT)), mkdir(fileparts(OUT)); end
FS = 60; FSLAB = 46; LWBOX = 5; CANV = 14.5; LWD = 5; MS = 26;
CR = [0.85 0.10 0.10];
sx = round(max(S.bs) / 4);  tx = (1:3)*sx;  xr = [0, 4*sx];   % 0..500 um ('tip') or 0..2000 um ('c500')
% vertical ticks: 0 .. 4*s with three ticks at s, 2s, 3s; s = smallest nice step covering the data
NS = [1 2 2.5 4 5 10 20 25 40 50 100 200 250 400 500 1000 2000 2500 4000 5000 10000 20000 25000 40000 50000];
nice = @(m) min(NS(4*NS >= 1.02*m));
% large tick numbers crowd out the axis title -> put the power of ten into the axis label
% (cost >= 1e5 -> units of 1e4; cost >= 1e3 -> units of 1e3)
SC = 1;  if max(S.C) >= 1e5, SC = 1e4;  elseif max(S.C) >= 1e3, SC = 1e3; end
C  = S.C / SC;  Cmin = S.Cmin / SC;
sy = nice(max(C));  ty = (1:3)*sy;  yr = [0 4*sy];

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax = axes(fig, 'Position', [0.22 0.15 0.71 0.78]); hold(ax, 'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
plot(ax, S.bs, C, '-', 'Color', [0.05 0.10 0.95], 'LineWidth', LWD);
plot(ax, [S.bmin S.bmin], [yr(1) Cmin], '--', 'Color', CR, 'LineWidth', 3);
plot(ax, S.bmin, Cmin, 'o', 'Color', CR, 'MarkerFaceColor', CR, 'MarkerSize', MS);
% minimum label: latex for the hat (tex has no \hat); concatenated, not sprintf (sprintf eats the backslashes)
% label right of the point, or left of it when the minimum sits in the right part of the axis
HA = 'left';  dxl = 0.03*diff(xr);  if S.bmin > xr(1) + 0.55*diff(xr), HA = 'right';  dxl = -dxl; end
text(ax, S.bmin + dxl, Cmin + 0.30*diff(yr), ['$\hat{\ell} = ' num2str(S.bmin, '%.1f') '$'], 'Interpreter','latex', ...
     'FontSize', FS*0.75, 'Color', CR, 'HorizontalAlignment',HA, 'VerticalAlignment','bottom');
xlim(ax, xr);  ylim(ax, yr);
set(ax, 'XTick', tx, 'YTick', ty);
box(ax, 'on');
set(ax, 'FontSize', FS, 'FontWeight','bold', 'LineWidth', LWBOX, 'TickDir','in', 'Layer','top');
ax.YAxis.Exponent = 0;
pbaspect(ax, [1 1 1]);
yoff = yr(1) - 0.035*diff(yr);
for xv = xr
    text(ax, xv, yoff, sprintf('%g', xv), 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
         'FontSize', FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{\hat{\ell}\;(\mu m)}$', 'Interpreter','latex', 'FontSize', FSLAB);
if SC == 1
    ylabel(ax, '$\mathbf{Cost\;(mT^{2})}$', 'Interpreter','latex', 'FontSize', FSLAB);
else
    ylabel(ax, ['$\mathbf{Cost\;(10^{' num2str(log10(SC)) '}\;mT^{2})}$'], 'Interpreter','latex', 'FontSize', FSLAB);
end
set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
print(fig, OUT, '-dpng', '-r200');
fprintf('saved %s\n', OUT);
