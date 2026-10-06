% sp_noise_hist_plot.m -- histogram of the fitted l_hat (= b) over the noise Monte Carlo (loads sp_noise_mc_N<NT>.mat only).
% Vertical axis = percentage of trials per bin; bin width BW.
if ~exist('NT','var'), NT = 10000; end
if ~exist('CASE','var'), CASE = 'single_pole'; end   % 'single_pole' | 'pair_pole'
clearvars -except NT CASE;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
DAT = fullfile(CALROOT, 'utils', 'data', 'NTU_hexapole');  FIG = fullfile(CALROOT, 'figures', 'NTU_hexapole');
PRE = struct('single_pole','sp', 'pair_pole','pp');  PRE = PRE.(CASE);
S   = load(fullfile(DAT, sprintf('%s_noise_mc_N%d.mat', PRE, NT)));
OUT = fullfile(FIG, CASE, 'noise_hist.png');
FS = 60; FSLAB = 46; LWBOX = 5; CANV = 14.5;
% horizontal: three ticks centred on round(mean), step sx = smallest nice value whose window
% [c-2sx, c+2sx] covers all samples; bin width = sx/12
c  = round(mean(S.B));  SX = [1 2 3 4 5 6 8 10 12 15 20 25 30 40 50];
sx = SX(find(c - 2*SX <= min(S.B) & c + 2*SX >= max(S.B), 1));
tx = c + (-1:1)*sx;  xr = [tx(1)-sx, tx(end)+sx];  BW = sx/12;
assert(min(S.B) >= xr(1) && max(S.B) <= xr(2), 'l_hat outside the plotted range');
edges = xr(1):BW:xr(2);
pc = histcounts(S.B, edges) / numel(S.B) * 100;        % percentage count
NS = [1 2 2.5 4 5 10 20 25];  sy = min(NS(4*NS >= 1.02*max(pc)));  ty = (1:3)*sy;  yr = [0 4*sy];

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax = axes(fig, 'Position', [0.22 0.15 0.71 0.78]); hold(ax, 'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
bar(ax, edges(1:end-1) + BW/2, pc, 1, 'FaceColor', [0.05 0.10 0.95], 'EdgeColor', 'none');
xlim(ax, xr);  ylim(ax, yr);
set(ax, 'XTick', tx, 'YTick', ty);
box(ax, 'on');
set(ax, 'FontSize', FS, 'FontWeight','bold', 'LineWidth', LWBOX, 'TickDir','in', 'Layer','top');
pbaspect(ax, [1 1 1]);
yoff = yr(1) - 0.035*diff(yr);
for xv = xr
    text(ax, xv, yoff, sprintf('%g', xv), 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
         'FontSize', FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{\hat{\ell}\;(\mu m)}$', 'Interpreter','latex', 'FontSize', FSLAB);
ylabel(ax, '$\mathbf{Percentage\;count\;(\%)}$', 'Interpreter','latex', 'FontSize', FSLAB);
set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
print(fig, OUT, '-dpng', '-r200');
fprintf('saved %s\n', OUT);
fprintf('N = %d | l_hat mean %.4f um, variance %.4f um^2 (std %.4f um) | a mean %.6g, variance %.6g (std %.6g) mT*um^2/A\n', ...
        numel(S.B), mean(S.B), var(S.B), std(S.B), mean(S.A), var(S.A), std(S.A));
fprintf('max bin %.2f %%, bins %d x %.1f um\n', max(pc), numel(pc), BW);
