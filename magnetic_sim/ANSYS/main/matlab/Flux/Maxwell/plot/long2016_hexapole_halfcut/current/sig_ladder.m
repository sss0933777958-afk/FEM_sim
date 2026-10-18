%% sig_ladder.m -- sensitivity eta_l_hat = sigma_l_hat / sigma vs number of calibration points
%  (long2016, single parameter, current; calib_noise_v2.pdf sec. 6, um/mT)
%  Reads only utils/data/long2016_hexapole_halfcut/calib_noise_ladder.mat
%  (written by utils/scripts/long2016_hexapole_halfcut/calib_noise.m with LADDER set).
%  Single parameter only (column 1), discrete points N = 7 ... NMAX (no legend, rule 9).
%  Frame 7 ... 301, round ticks 75 / 150 / 225 (end gaps 68 / 76, user request); the N = 301 point sits on the frame.
%  Vertical: frame 25 ... 125, ticks 50 / 75 / 100 (round numbers, user request).
%  Output: figures/long2016_hexapole_halfcut/hexapole/sig_ladder.png

here = fileparts(mfilename('fullpath'));                           % .../plot/<model>/current
CAL  = fileparts(fileparts(fileparts(here)));                      % .../Flux/Maxwell
S    = load(fullfile(CAL, 'utils', 'data', 'long2016_hexapole_halfcut', 'calib_noise_ladder.mat'));
OUTD = fullfile(CAL, 'figures', S.MODEL, 'hexapole');

% ---- style (figure-style.md) ----------------------------------------------
FS = 60;  FSLAB = 36;  MS = 16;  LWBOX = 5.0;  CANV = 14.5;
CB = [0.05 0.10 0.95];
NMAX = 301;                                  % last plotted level

k  = S.N(:) <= NMAX;
N  = S.N(k).';   s = S.SIG_U(k, 1);          % eta_l_hat [um/mT], single parameter
X0 = 7;   X1 = NMAX;
xt = [75 150 225];                           % round ticks, equally spaced

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes(fig);   hold(ax,'on');
try, ax.Toolbar = []; end                                            %#ok<TRYNC>
try, ax.Interactions = []; end                                       %#ok<TRYNC>
plot(ax, N, s, 'o', 'MarkerSize',MS, 'MarkerFaceColor',CB, 'MarkerEdgeColor',CB, 'Clipping','off');

box(ax,'on');  grid(ax,'off');
set(ax, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, 'TickLength',[.02 .02]);
xlim(ax, [X0 X1]);   set(ax, 'XTick', xt, 'XTickLabel', compose('%g', xt), 'XTickLabelRotation', 0);
yl = [25 125];   yt = [50 75 100];
ylim(ax, yl);   set(ax, 'YTick', yt, 'YTickLabel', compose('%g', yt));

yoff = yl(1) - 0.022*diff(yl);
for xv = [X0 X1]
    text(ax, xv, yoff, sprintf('%d', xv), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold');
end
xlabel(ax, '$\mathbf{Number\;of\;calibration\;points}$', 'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{\eta_{\hat{\ell}}\;(\mu m/mT)}$', 'Interpreter','latex', 'FontSize',FSLAB);

out = fullfile(OUTD, 'sig_ladder.png');
set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');
close(fig);
fprintf('wrote %s | ylim [%g %g] ticks %s\n', out, yl, mat2str(yt));
