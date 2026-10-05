%% eps_norm.m -- normalised shell residual eps(r_k)/eps(R) vs r/R, single and eighteen (long2016, current)
%  Reads only utils/data/long2016_hexapole_halfcut/eps_shell.mat (written by utils/scripts/long2016_hexapole_halfcut/eps_shell.m); the normalisation
%  is a plain division done here for presentation.
%  Output: figures/long2016_hexapole_halfcut/current/eps_norm.png

here = fileparts(mfilename('fullpath'));                           % .../plot/<model>/current
CAL  = fileparts(fileparts(fileparts(here)));                      % .../Flux/Maxwell
S    = load(fullfile(CAL, 'utils', 'data', 'long2016_hexapole_halfcut', 'eps_shell.mat'));
OUT  = fullfile(CAL, 'figures', S.MODEL, 'current', 'eps_norm.png');

% ---- style (figure-style.md) ----------------------------------------------
FS = 60;  FSLAB = 36;  FSLEG = 45;  LW = 5.0;  LWBOX = 5.0;  CANV = 14.5;
c1 = [0.05 0.10 0.95];   c2 = [0.85 0.10 0.10];

x = [0, S.rk(:).'] / S.R;                                          % r/R; 0 = centre point
E = [S.e_centre; S.EPS];   E = E ./ E(end, :);                     % eps(r_k)/eps(R)
xt = [0.25 0.5 0.75];                                               % three ticks, equidistant to 0 and 1 (user choice)
YL = [0 1.2];   yt = [0.3 0.6 0.9];

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes(fig);   hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
p1 = plot(ax, x, E(:,1), '-', 'Color',c1, 'LineWidth',LW);
p2 = plot(ax, x, E(:,2), '-', 'Color',c2, 'LineWidth',LW);

box(ax,'on');  grid(ax,'off');
set(ax, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, 'TickLength',[.02 .02]);
xlim(ax, [0 1]);   set(ax, 'XTick', xt, 'XTickLabel', compose('%g', xt));
ylim(ax, YL);      set(ax, 'YTick', yt, 'YTickLabel', compose('%g', yt));

yoff = YL(1) - 0.022*diff(YL);
for xv = [0 1]
    text(ax, xv, yoff, sprintf('%g', xv), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{r/R}$', 'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{\varepsilon(r_k)/\varepsilon(R)}$', 'Interpreter','latex', 'FontSize',FSLAB);

lg = legend(ax, [p1 p2], {'Single parameter', 'Eighteen parameters'}, ...
            'Interpreter','tex', 'Location','northoutside', 'NumColumns',2);
lg.FontSize = FSLEG;   lg.FontWeight = 'bold';   lg.ItemTokenSize = [55 25];
lg.Box = 'on';   lg.EdgeColor = 'k';   lg.LineWidth = LWBOX;
drawnow;
axp = get(ax,'Position');   lgp = get(lg,'Position');
GAPN = 0.022;   newTop = 1 - lgp(4) - GAPN - 0.006;
axp(4) = newTop - axp(2);   set(ax,'Position',axp);
assert(lgp(3) < 0.98, 'legend wider than the canvas');
x0 = min(max(axp(1) + (axp(3)-lgp(3))/2, 0.01), 0.99 - lgp(3));
set(lg, 'Position', [x0, newTop + GAPN, lgp(3), lgp(4)]);

set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
print(fig, OUT, '-dpng', '-r200');
fprintf('wrote %s | eps(0)/eps(R): single %.4f, eighteen %.4f\n', OUT, E(1,1), E(1,2));
