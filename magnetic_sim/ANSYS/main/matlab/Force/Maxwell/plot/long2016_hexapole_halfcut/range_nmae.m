% range_nmae  Test NMAE against the data range, for both field models.
%
%   For every radius R the nodes inside r <= R are split 80 / 20 (rng 0, the same
%   split for both models, training set capped at 4000 nodes from R = 250 on so the
%   RBF stays solvable; the test set is always the full 20 %).  Each model is then
%   tuned in TWO STAGES -- first the smoothness gate, then accuracy:
%       stage 1  d3(b.b)/ds^3 must keep a constant sign on all three actuator axes
%                over |s| <= R; parameters that fail are discarded outright
%       stage 2  among the survivors, take the lowest TEST NMAE
%   Parameters searched:  RBF  a (rho, lambda) grid    sph  the degree L
%   The scan that writes range_scan.mat is temp_code/range_scan.m.
%
%   ⚠ Beyond R ~ 200 um the smoothness test itself stops being the right question:
%   the true field has an inflection about 185 um out on the side AWAY from the
%   excited pole, so a correct model is rejected there and only over-smoothed ones
%   survive (the harmonic degree collapses to L = 6 and then 2; the RBF finds no
%   admissible pair at all at R = 450 and 500).  The rise on the right of these
%   figures therefore measures the criterion, not the models.
%
%   Reads only (plot-scripts-pure); range_scan.mat is written by the scan.
%
% Output: figures/long2016_hexapole_halfcut/range_nmae.png   (both models, log y)
%         figures/long2016_hexapole_halfcut/range_diff.png   (sph - RBF, signed)

HERE = fileparts(mfilename('fullpath'));
FMX  = fileparts(fileparts(HERE));
DAT  = fullfile(FMX,'utils','data');
FIG  = fullfile(FMX,'figures','long2016_hexapole_halfcut','current');

S = load(fullfile(DAT, 'range_scan.mat'));
R  = [S.RES.R].';
nr = arrayfun(@(s) s.rbf.nm, S.RES).';          % Inf where nothing was smooth
ns = arrayfun(@(s) s.sph.nm, S.RES).';
nr(~isfinite(nr)) = NaN;   ns(~isfinite(ns)) = NaN;
dd = ns - nr;                                   % signed, NaN where the RBF has no entry

fprintf('%6s %10s %10s %10s\n','R','RBF','sph','sph-RBF');
for q = 1:numel(R)
    fprintf('%6d %9.4f%% %9.4f%% %+9.4f\n', R(q), nr(q), ns(q), dd(q));
end
fprintf('RBF has no smooth pair at R = %s um\n', mat2str(R(isnan(nr)).'));

FS = 60;  FSLEG = 45;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;  LW = 6;  MS = 20;
cR = [0.85 0.10 0.10];   cS = [0.05 0.10 0.95];
% [MODIFIED 2026-09-16] the frame now ENDS on the data: xlim = [min R, max R], so the
%   first and last samples sit exactly on the left and right box edges (user's call --
%   it overrides figure-style rule 5's "end padding = tick spacing", which cannot hold
%   at the same time).  Inner ticks stay odd in number and evenly spaced.
XR = [min(R) max(R)];   XT = [100 250 400];

% ---------------- figure 1: both curves, log y (three decades apart) -------------
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]); hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
% [MODIFIED 2026-09-16] markers only, no connecting line: each point is an independent
%   scan at that radius, and the line invited reading a trend between the samples.
%   Clipping off so the end markers are drawn whole where they straddle the box edge.
h1 = plot(ax, R, nr, 'o', 'LineStyle','none', 'Color',cR, 'LineWidth',LW*0.6, ...
          'MarkerSize',MS, 'MarkerFaceColor',cR, 'Clipping','off');
h2 = plot(ax, R, ns, 'o', 'LineStyle','none', 'Color',cS, 'LineWidth',LW*0.6, ...
          'MarkerSize',MS, 'MarkerFaceColor',cS, 'Clipping','off');
set(ax,'YScale','log');
set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
       'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
       'XTickLabelRotation',0,'YTickLabelRotation',0);
xlim(ax, XR);   set(ax,'XTick', XT);
YR = [0.01 100];   YT = [0.1 1 10];                    % one decade of margin each side
ylim(ax, YR);   set(ax,'YTick', YT, 'YMinorTick','off');
ax.YAxis.MinorTickValues = [];
for xv = XR
    text(ax, xv, YR(1)*10^(-0.035*log10(YR(2)/YR(1))), sprintf('%g', xv), ...
         'HorizontalAlignment','center','VerticalAlignment','top', ...
         'FontSize',FS,'FontWeight','bold','Clipping','off');
end
xlabel(ax, '$\mathbf{R\;(\mu m)}$',          'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{Test\;NMAE\;(\%)}$',    'Interpreter','latex', 'FontSize',FSLAB);
lg = legend(ax, [h1 h2], {'RBF', 'Spherical harmonics'}, 'Interpreter','tex', ...
            'Location','northwest', 'NumColumns',1);
lg.FontSize = FSLEG;  lg.FontWeight = 'bold';
lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.Color = 'w';
lg.ItemTokenSize = [55 25];
hold(ax,'off');
out = fullfile(FIG,'range_nmae.png');
set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');   close(fig);
fprintf(['wrote %s' newline], out);

% ---------------- figure 2: the signed difference -------------------------------
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]); hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
yline(ax, 0, 'k--', 'LineWidth', LWBOX*0.7);
plot(ax, R, dd, 'o', 'LineStyle','none', 'Color',cS, 'LineWidth',LW*0.6, ...
     'MarkerSize',MS, 'MarkerFaceColor',cS, 'Clipping','off');
set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
       'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
       'XTickLabelRotation',0,'YTickLabelRotation',0);
xlim(ax, XR);   set(ax,'XTick', XT);
[YR2, YT2] = ylim_fit_(min(dd), max(dd), 3);
ylim(ax, YR2);  set(ax,'YTick', YT2);
for xv = XR
    text(ax, xv, YR2(1)-0.022*diff(YR2), sprintf('%g', xv), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{R\;(\mu m)}$', 'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{NMAE_{sph} - NMAE_{RBF}\;(\%)}$', 'Interpreter','latex', 'FontSize',FSLAB);
hold(ax,'off');
out = fullfile(FIG,'range_diff.png');
set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');   close(fig);
fprintf(['wrote %s  (y ticks %s, frame [%g %g])' newline], out, mat2str(YT2), YR2);

% ---- local helpers ------------------------------------------------------------
function [lim, tk] = ylim_fit_(lo, hi, N)
    if nargin < 3 || isempty(N), N = 3; end
    half = (N+1)/2;   need = (hi - lo)/(N+1);
    for kk = floor(log10(need)) : floor(log10(need))+1
        for m = [1 1.5 2 2.5 3 4 5 6 8]
            s = m*10^kk;
            if s < need*(1+1e-9), continue, end
            u = 10^floor(log10(s));
            while abs(s/u - round(s/u)) > 1e-9, u = u/2;  end
            c = ceil((hi - half*s)/u - 1e-9) * u;
            if c - half*s <= lo + 1e-9
                lim = [c-half*s, c+half*s];   tk = c + (-(N-1)/2:(N-1)/2)*s;   return
            end
        end
    end
    error('ylim_fit_:none','no nice frame for [%g %g]', lo, hi);
end
