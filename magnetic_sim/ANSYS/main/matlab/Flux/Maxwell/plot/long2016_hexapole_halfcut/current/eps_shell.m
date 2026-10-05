%% eps_shell.m -- shell-averaged residual eps(r_k) vs r, single and eighteen (long2016, current)
%  Reads only utils/data/long2016_hexapole_halfcut/eps_shell.mat (written by utils/scripts/long2016_hexapole_halfcut/eps_shell.m).
%  Output: figures/long2016_hexapole_halfcut/current/eps_shell.png
%  The two models differ by ~100x, so the y axis is logarithmic.

here = fileparts(mfilename('fullpath'));                           % .../plot/<model>/current
CAL  = fileparts(fileparts(fileparts(here)));                      % .../Flux/Maxwell
S    = load(fullfile(CAL, 'utils', 'data', 'long2016_hexapole_halfcut', 'eps_shell.mat'));
OUT  = fullfile(CAL, 'figures', S.MODEL, 'current', 'eps_shell.png');

% ---- style (figure-style.md) ----------------------------------------------
FS = 60;  FSLAB = 36;  FSLEG = 45;  LW = 5.0;  LWBOX = 5.0;  CANV = 14.5;
c1 = [0.05 0.10 0.95];   c2 = [0.85 0.10 0.10];

r  = [0, S.rk(:).'] * 1e6;                                         % [um]; r = 0 is the centre point (its own value)
E  = [S.e_centre; S.EPS];
X0 = 0;   X1 = 150;   xt = X0 + (1:5) * (X1 - X0)/6;               % 25 50 75 100 125
YL = [1e-5 1e1];      yt = 10.^(-4:0);                             % five decades inside, equal spacing to both edges

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes(fig);   hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
p1 = plot(ax, r, E(:,1), '-', 'Color',c1, 'LineWidth',LW);
p2 = plot(ax, r, E(:,2), '-', 'Color',c2, 'LineWidth',LW);

box(ax,'on');  grid(ax,'off');
set(ax, 'YScale','log', 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, 'TickLength',[.02 .02]);
xlim(ax, [X0 X1]);   set(ax, 'XTick', xt, 'XTickLabel', compose('%d', xt));
ylim(ax, YL);   set(ax, 'YTick', yt, 'YTickLabel', compose('10^{%d}', log10(yt)), 'YMinorTick','off');

yoff = YL(1) / 10^(0.022 * log10(YL(2)/YL(1)));                    % same 2.2 % offset, in log space
for xv = [X0 X1]
    text(ax, xv, yoff, sprintf('%d', xv), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{r\;(\mu m)}$', 'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{\varepsilon(r_k)\;(mT^{2})}$', 'Interpreter','latex', 'FontSize',FSLAB);

% legend above the axes (inside it hides tick marks), entries side by side
lg = legend(ax, [p1 p2], {'Single parameter', 'Eighteen parameters'}, ...
            'Interpreter','tex', 'Location','northoutside', 'NumColumns',2);
lg.FontSize = FSLEG;   lg.FontWeight = 'bold';   lg.ItemTokenSize = [55 25];
lg.Box = 'on';   lg.EdgeColor = 'k';   lg.LineWidth = LWBOX;
drawnow;
axp = get(ax,'Position');   lgp = get(lg,'Position');
GAPN = 0.022;   newTop = 1 - lgp(4) - GAPN - 0.006;
axp(4) = newTop - axp(2);   set(ax,'Position',axp);
assert(lgp(3) < 0.98, 'legend wider than the canvas');
x0 = min(max(axp(1) + (axp(3)-lgp(3))/2, 0.01), 0.99 - lgp(3));   % centred on the axes, kept inside the canvas
set(lg, 'Position', [x0, newTop + GAPN, lgp(3), lgp(4)]);

set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
print(fig, OUT, '-dpng', '-r200');
fprintf('wrote %s\n', OUT);
