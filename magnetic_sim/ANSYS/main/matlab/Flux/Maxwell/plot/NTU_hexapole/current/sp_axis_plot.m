% sp_axis_plot.m -- |b_x| on the 17 axis samples and the fitted H = g*I/(x - x_c)^2 vs distance from the tip
% (loads sp_axis.mat only; the fit is done in sp_axis_fit.m).
if ~exist('CASE','var'), CASE = 'single_pole'; end   % 'single_pole' | 'pair_pole' (NTU) | 'zhi_peng' (pair, P1)
clearvars -except CASE;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
% case -> {data folder model, .mat, figure folder}
CS  = struct('single_pole', {{'NTU_hexapole', 'sp_axis.mat', fullfile('NTU_hexapole', 'single_pole')}}, ...
             'pair_pole',   {{'NTU_hexapole', 'pp_axis.mat', fullfile('NTU_hexapole', 'pair_pole')}}, ...
             'zhi_peng',    {{'zhi_peng',     'zp_axis.mat', fullfile('zhi_peng', 'pair_pole')}});
CS  = CS.(CASE);
S   = load(fullfile(CALROOT, 'utils', 'data', CS{1}, CS{2}));
OUT = fullfile(CALROOT, 'figures', CS{3}, 'b_axis.png');
if ~isfolder(fileparts(OUT)), mkdir(fileparts(OUT)); end
FS = 60; FSLEG = 45; FSLAB = 46; LWBOX = 5; CANV = 14.5; LWD = 4; MS = 22;
tx = [100 200 300];  sx = tx(2)-tx(1);  xr = [tx(1)-sx, tx(end)+sx];
% vertical ticks: 0 .. 4*s with three ticks at s, 2s, 3s; s = smallest nice step covering the data
NS = [1 2 2.5 4 5 10 20 25 40 50 100 200 250 400 500 1000 2000 2500 4000 5000];
nice = @(m) min(NS(4*NS >= 1.02*m));
sy = nice(max(abs(S.b(:,1))));  ty = (1:3)*sy;  yr = [0 4*sy];

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax = axes(fig, 'Position', [0.22 0.15 0.71 0.78]); hold(ax, 'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
hC = plot(ax, S.dc, S.Hc, '-', 'Color', [0.85 0.10 0.10], 'LineWidth', LWD);
hD = plot(ax, S.d*1e6, abs(S.b(:,1)), 'o', 'Color', [0.05 0.10 0.95], 'MarkerFaceColor', 'none', ...
     'MarkerSize', MS, 'LineWidth', LWD);
xlim(ax, xr);  ylim(ax, yr);
set(ax, 'XTick', tx, 'YTick', ty);
box(ax, 'on');
set(ax, 'FontSize', FS, 'FontWeight','bold', 'LineWidth', LWBOX, 'TickDir','in', 'Layer','top');
pbaspect(ax, [1 1 1]);
lg = legend(ax, [hD hC], {'FEM', 'H(x)'}, 'Location', 'northeast', 'FontSize', FSLEG, 'FontWeight', 'bold');
lg.LineWidth = LWBOX;  lg.ItemTokenSize = [55 25];
yoff = yr(1) - 0.035*diff(yr);
for xv = xr
    text(ax, xv, yoff, sprintf('%g', xv), 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
         'FontSize', FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{Distance\;from\;tip\;(\mu m)}$', 'Interpreter','latex', 'FontSize', FSLAB);
ylabel(ax, '$\mathbf{|b_x|\;(mT)}$', 'Interpreter','latex', 'FontSize', FSLAB);
set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
print(fig, OUT, '-dpng', '-r200');
fprintf('saved %s\n', OUT);
