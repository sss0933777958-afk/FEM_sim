% cost_c500_plot.m -- profile cost J(l_hat) of the 1-D model for the centre samples (c500, 350..650 um from the
% P1 tip): single pole and pair pole (P1 excitation) overlaid, l_hat = 0 .. 6000 um.  The two costs differ by
% about ten times, so each curve is shown as J / J_min (touches 1 at its own minimum, marked by a dot) on a log axis.
% Loads sp_cost_b_c500.mat and pp_cost_b_c500.mat only (written by utils/scripts/NTU_hexapole/sp_cost_b.m).
% Output (same figure in both folders): figures/NTU_hexapole/{single_pole,pair_pole}/cost_c500.png
% RAW = true: J itself [mT^2] (not normalised), linear axis, one figure per case (single series, no legend)
%   -> figures/NTU_hexapole/<case>/cost_J.png
%   QTY = 'J' (default) | 'dJ': with RAW, plot dJ/dl_hat at a* [mT^2/um] (field dC) instead -> cost_dJ.png
if ~exist('RAW','var'), RAW = false; end
if ~exist('QTY','var'), QTY = 'J'; end
clearvars -except RAW QTY;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
DAT  = fullfile(CALROOT, 'utils', 'data', 'NTU_hexapole');
FIG  = fullfile(CALROOT, 'figures', 'NTU_hexapole');
S = {load(fullfile(DAT, 'sp_cost_b_c500.mat')), load(fullfile(DAT, 'pp_cost_b_c500.mat'))};

% ---- style (figure-style.md) ---------------------------------------------------------
FS = 60;  FSLEG = 45;  FSLAB = 46;  LWBOX = 5.0;  CANV = 14.5;  LW = 5.0;  MS = 30;
CB = [0.05 0.10 0.95];  CR = [0.85 0.10 0.10];  COL = {CB, CR};  NAME = {'Single pole', 'Pair pole'};
XMAX = max(S{1}.bs(end), S{2}.bs(end));
sx = XMAX/4;  tx = (1:3)*sx;  xr = [0 XMAX];                         % 0 .. 6000: ticks 1500 / 3000 / 4500
ty = [1 10 100];  yr = [0.1 1000];                                    % log axis: one decade between ticks and ends

if RAW
    SUB = {'single_pole', 'pair_pole'};
    for k = 1:2
        if strcmp(QTY, 'dJ')
            Y = S{k}.dC;  Ymin = 0;  [yrk, tyk] = band(min(Y), max(Y));
        else
            Y = S{k}.C;   Ymin = S{k}.Cmin;
            s = nice_step(max(Y)/4);  tyk = (1:3)*s;  yrk = [0 4*s];      % 0 .. 4s, ticks s / 2s / 3s
        end
        fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
        ax  = axes(fig, 'Position', [0.20 0.15 0.73 0.80]);  hold(ax, 'on');
        try, ax.Toolbar = []; end
        try, ax.Interactions = []; end
        if strcmp(QTY, 'dJ'), plot(ax, xr, [0 0], '-', 'Color',[0.6 0.6 0.6], 'LineWidth',2); end   % zero
        plot(ax, S{k}.bs, Y, '-', 'Color',COL{k}, 'LineWidth',LW);
        plot(ax, S{k}.bmin, Ymin, '.', 'Color',COL{k}, 'MarkerSize',MS*2.2, 'Clipping','off');
        xlim(ax, xr);  ylim(ax, yrk);
        set(ax, 'XTick', tx, 'XTickLabel', compose('%d', tx), 'YTick', tyk, 'YTickLabel', compose('%g', tyk));
        box(ax, 'on');  grid(ax, 'off');
        set(ax, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, 'TickDir','in', 'Layer','top');
        yoff = yrk(1) - 0.025*diff(yrk);
        for xv = xr
            text(ax, xv, yoff, sprintf('%d', xv), 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
                 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
        end
        xlabel(ax, '$\mathbf{\hat{\ell}\;(\mu m)}$', 'Interpreter','latex', 'FontSize',FSLAB);
        if strcmp(QTY, 'dJ')
            ylabel(ax, '$\mathbf{\partial J/\partial\hat{\ell}\;(mT^2/\mu m)}$', 'Interpreter','latex', 'FontSize',FSLAB);
        else
            ylabel(ax, '$\mathbf{J\;(mT^2)}$', 'Interpreter','latex', 'FontSize',FSLAB);
        end
        set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
        out = fullfile(FIG, SUB{k}, ['cost_' QTY '.png']);
        print(fig, out, '-dpng', '-r200');  close(fig);
        fprintf('wrote %s | l_min %.1f um | %s(0) %.4g, %s(%g) %.4g, range [%.4g %.4g] | y %s\n', out, S{k}.bmin, ...
                QTY, Y(1), QTY, S{k}.bs(end), Y(end), min(Y), max(Y), mat2str(yrk));
    end
    return
end

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes(fig, 'Position', [0.20 0.15 0.73 0.80]);  hold(ax, 'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
h = gobjects(1,2);
for k = 1:2
    h(k) = plot(ax, S{k}.bs, S{k}.C / S{k}.Cmin, '-', 'Color',COL{k}, 'LineWidth',LW);
    plot(ax, S{k}.bmin, 1, '.', 'Color',COL{k}, 'MarkerSize',MS*2.2, 'Clipping','off');
    fprintf('%-11s l_min %.1f um | J_min %.4g mT^2 | J/J_min at 0: %.1f, at %g: %.2f\n', NAME{k}, S{k}.bmin, ...
            S{k}.Cmin, S{k}.C(1)/S{k}.Cmin, S{k}.bs(end), S{k}.C(end)/S{k}.Cmin);
end
set(ax, 'YScale','log');
xlim(ax, xr);  ylim(ax, yr);
set(ax, 'XTick', tx, 'XTickLabel', compose('%d', tx), 'YTick', ty, 'YTickLabel', compose('%g', ty), ...
        'YMinorTick','off');
box(ax, 'on');  grid(ax, 'off');
set(ax, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, 'TickDir','in', 'Layer','top');
yoff = yr(1) / 10^(0.025*log10(yr(2)/yr(1)));                         % just below the frame on the log axis
for xv = xr
    text(ax, xv, yoff, sprintf('%d', xv), 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
         'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{\hat{\ell}\;(\mu m)}$', 'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{J\,/\,J_{min}}$', 'Interpreter','latex', 'FontSize',FSLAB);
lg = legend(ax, h, NAME, 'Location','southeast', 'FontSize',FSLEG, 'FontWeight','bold');
lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.ItemTokenSize = [55 25];

set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
for sub = {'single_pole', 'pair_pole'}
    out = fullfile(FIG, sub{1}, 'cost_c500.png');
    print(fig, out, '-dpng', '-r200');
    fprintf('wrote %s\n', out);
end
close(fig);

% ---------------------------------------------------------------------------------------
function s = nice_step(m)
% Smallest step from {1 1.2 1.5 2 2.5 3 4 5 6 8} x 10^e with s >= 1.02*m.
    e = floor(log10(m));  c = [1 1.2 1.5 2 2.5 3 4 5 6 8 10]*10^e;
    s = c(find(c >= 1.02*m, 1));
end

function [yl, tk] = band(lo, hi)
% Three equally spaced ticks with end gaps equal to the step; frame bottom on a multiple of the step.
    R = hi - lo;  mg = 0.02*R;  e = floor(log10(R/4));  best = inf;  yl = [];  tk = [];
    for ee = e-1:e+1
        for c = [1 1.2 1.5 2 2.5 3 4 5 6 8]
            s  = round(c*10^ee, 12);
            y0 = floor((lo - mg)/s)*s;
            if y0 + 4*s < hi + mg, continue; end
            if 4*s < best, best = 4*s;  yl = round([y0, y0 + 4*s], 12);  tk = round(y0 + (1:3)*s, 12); end
        end
    end
end
