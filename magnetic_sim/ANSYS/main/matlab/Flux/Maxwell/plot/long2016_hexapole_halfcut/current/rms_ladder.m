%% rms_ladder.m -- calibration RMS vs number of calibration points (long2016, current)
%  Reads only utils/data/long2016_hexapole_halfcut/rms_ladder.mat (written by utils/scripts/long2016_hexapole_halfcut/rms_ladder.m).
%  One figure per model (single series -> no legend):
%    figures/long2016_hexapole_halfcut/current/rms_single.png
%    figures/long2016_hexapole_halfcut/current/rms_eighteen.png

here = fileparts(mfilename('fullpath'));                           % .../plot/<model>/current
CAL  = fileparts(fileparts(fileparts(here)));                      % .../Flux/Maxwell
S    = load(fullfile(CAL, 'utils', 'data', 'long2016_hexapole_halfcut', 'rms_ladder.mat'));
OUTD = fullfile(CAL, 'figures', S.MODEL, 'current');

% ---- style (figure-style.md) ----------------------------------------------
FS = 60;  FSLAB = 36;  LW = 5.0;  LWBOX = 5.0;  CANV = 14.5;
COL  = {[0.05 0.10 0.95], [0.85 0.10 0.10]};
NAME = {'single', 'eighteen'};

X0 = 7;   X1 = 301;
xt = X0 + (1:5) * (X1 - X0)/6;               % 56 105 154 203 252: five integer ticks, spacing 49 to both ends
N  = S.N(:).';
k  = find(N > X1, 1);                        % draw through the first rung past X1; the axes clip it at X1
sel = 1:k;   inb = N <= X1;

for m = 1:2
    r = S.RMS(:,m).';
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes(fig);   hold(ax,'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    plot(ax, N(sel), r(sel), '-', 'Color',COL{m}, 'LineWidth',LW, 'Clipping','on');

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

    out = fullfile(OUTD, ['rms_' NAME{m} '.png']);
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
    print(fig, out, '-dpng', '-r200');
    close(fig);
    fprintf('wrote %s | ylim [%g %g] ticks %s\n', out, yl, mat2str(yt));
end

% ---------------------------------------------------------------------------
function [yl, tk] = ylim_band(lo, hi)
% Three equally spaced ticks, spacing to both frame edges equal to the step; frame
% edges are not labelled. The frame hugs the data: the bottom is the largest
% multiple of step/2 below the data, and the (step, bottom) pair with the
% tightest frame wins.
    R  = hi - lo;   mg = 0.02 * R;
    e  = floor(log10(R / 4));
    best = inf;   yl = [];   tk = [];
    for ee = e-1:e+1
        for c = [1 1.2 1.5 2 2.5 3 4 5 6 8]
            s  = round(c * 10^ee, 12);
            y0 = floor((lo - mg) / (s/2)) * (s/2);
            if y0 + 4*s < hi + mg, continue; end
            if 4*s < best
                best = 4*s;   yl = round([y0, y0 + 4*s], 12);   tk = round(y0 + (1:3)*s, 12);
            end
        end
    end
end
