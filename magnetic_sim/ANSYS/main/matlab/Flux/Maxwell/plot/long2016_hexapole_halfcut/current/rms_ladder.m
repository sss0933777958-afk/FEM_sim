%% rms_ladder.m -- calibration RMS vs number of calibration points (long2016, current)
%  Reads only utils/data/long2016_hexapole_halfcut/rms_ladder.mat (written by utils/scripts/long2016_hexapole_halfcut/rms_ladder.m).
%  One figure per model (single series -> no legend):
%    figures/long2016_hexapole_halfcut/hexapole/rms_single.png
%    figures/long2016_hexapole_halfcut/hexapole/rms_eighteen.png

here = fileparts(mfilename('fullpath'));                           % .../plot/<model>/current
CAL  = fileparts(fileparts(fileparts(here)));                      % .../Flux/Maxwell
S    = load(fullfile(CAL, 'utils', 'data', 'long2016_hexapole_halfcut', 'rms_ladder.mat'));
OUTD = fullfile(CAL, 'figures', S.MODEL, 'hexapole');

% ---- style (figure-style.md) ----------------------------------------------
FS = 60;  FSLAB = 36;  LW = 10.0;  LWBOX = 5.0;  CANV = 14.5;
MS = 16;                                     % marker size of the points-only version
STYLE = {'line', 'pts'};                     % 'line' -> rms_<m>.png, 'pts' -> rms_<m>_pts.png (one marker per rung)
COL  = {[0.05 0.10 0.95], [0.85 0.10 0.10]};
NAME = {'single', 'eighteen'};

X0 = 7;   X1 = 301;
xt = [75 150 225];                           % round ticks, equally spaced (end gaps 68 / 76, user request)
N  = S.N(:).';
k  = find(N > X1, 1);                        % draw through the first rung past X1; the axes clip it at X1
sel = 1:k;   inb = N <= X1;

for m = 1:2
  for st = STYLE
    r = S.RMS(:,m).';
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes(fig);   hold(ax,'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    if strcmp(st{1}, 'pts')
        plot(ax, N(inb), r(inb), 'o', 'Color',COL{m}, 'MarkerFaceColor',COL{m}, ...
             'MarkerSize',MS, 'LineStyle','none', 'Clipping','off');
        sfx = '_pts';
    else
        plot(ax, N(sel), r(sel), '-', 'Color',COL{m}, 'LineWidth',LW, 'Clipping','on');
        sfx = '';
    end

    box(ax,'on');  grid(ax,'off');
    set(ax, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, 'TickLength',[.02 .02]);
    xlim(ax, [X0 X1]);   set(ax, 'XTick', xt, 'XTickLabel', compose('%d', xt));
    [yl, yt] = ylim_band(min(r(inb)), max(r(inb)));
    ylim(ax, yl);   set(ax, 'YTick', yt, 'YTickLabel', compose('%g', yt));

    yoff = yl(1) - 0.022*diff(yl);
    for xv = [X0 X1]
        text(ax, xv, yoff, sprintf('%d', xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    xlabel(ax, '$\mathbf{Number\;of\;calibration\;points}$', 'Interpreter','latex', 'FontSize',FSLAB);
    ylabel(ax, '$\mathbf{RMS\;(mT)}$', 'Interpreter','latex', 'FontSize',FSLAB);

    out = fullfile(OUTD, ['rms_' NAME{m} sfx '.png']);
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
    print(fig, out, '-dpng', '-r200');
    close(fig);
    fprintf('wrote %s | ylim [%g %g] ticks %s\n', out, yl, mat2str(yt));
  end
end

% ---------------------------------------------------------------------------
function [yl, tk] = ylim_band(lo, hi)
% Three equally spaced ticks, spacing to both frame edges equal to the step; frame edges
% are not labelled. Step 1-5 x 10^e; the frame bottom is any multiple of 10^e, but every
% tick must have at most two significant digits (e.g. 0.17, 0.013; user request).
% Smallest frame wins; among equal frames, the one centred best on the data.
    mg = 0.03 * (hi - lo);   e = floor(log10((hi - lo) / 4));
    best = [inf inf];   yl = [];   tk = [];
    for ee = e-1:e+1
        for c = 1:5
            s = c * 10^ee;   q = 10^ee;
            for y0 = (floor((hi + mg)/q) - 4*c : floor((lo - mg)/q)) * q
                if y0 + 4*s < hi + mg, continue; end
                t = y0 + (1:3)*s;
                if any(abs(t) > 0 & abs(round(t ./ 10.^(floor(log10(abs(t))) - 1)) - ...
                        t ./ 10.^(floor(log10(abs(t))) - 1)) > 1e-6), continue; end
                sc = [4*s, abs(y0 + 2*s - (lo + hi)/2)];
                if sc(1) < best(1) - 1e-12 || (abs(sc(1) - best(1)) < 1e-12 && sc(2) < best(2))
                    best = sc;   yl = round([y0, y0 + 4*s], 12);   tk = round(y0 + (1:3)*s, 12);
                end
            end
        end
    end
end
