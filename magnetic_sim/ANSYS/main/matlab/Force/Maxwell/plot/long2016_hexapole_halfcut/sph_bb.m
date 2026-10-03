% sph_bb  b.b on the x_a axis from the solid-harmonic model fitted at R <= 150 um.
%
%   The model is the one the R = 150 sweep selected: degree L = 9 (K = 99), chosen
%   as the HIGHEST degree still passing the smoothness gate (d3(b.b)/ds^3 keeps one
%   sign on all three actuator axes over |s| <= 150), with patience 10.  NMAE plays
%   no part in that choice for the harmonics -- it is reported, not used.
%
%   Drawn on the same 1001 equidistant points over |s| <= 150 um as every other
%   curve in this folder, P1 excited at 1 A, so it overlays the FEM ones directly.
%   Single series, so no legend.
%
%   Reads only (plot-scripts-pure):
%       x150.mat     bbS, s      written by temp_code/x150.m
%       sph_v150.mat RES.nmae    written by temp_code/sph_strat.m (RUN_R = 150)
%
% Output: figures/long2016_hexapole_halfcut/sph_bb_R150.png
%
% [ADDED 2026-09-21]

HERE = fileparts(mfilename('fullpath'));
FMX  = fileparts(fileparts(HERE));
DAT  = fullfile(FMX,'utils','data');
FIG  = fullfile(FMX,'figures','long2016_hexapole_halfcut','current');

% [ADDED 2026-09-29] PLOT_OVR.XFILE draws another record instead of x150.mat: a
%   sphax_R<R>_L<L>.mat written by utils/sph_axes.m (same field names). Such a
%   record is self-contained, so sph_v150.mat is not read for it.
%   PLOT_OVR.YCAP = [lo hi] fixes the vertical window of every panel; without
%   it the window follows the data.
if exist('PLOT_OVR','var'), OVR_ = PLOT_OVR; else, OVR_ = struct(); end
XFILE = 'x150.mat';   if isfield(OVR_,'XFILE'), XFILE = OVR_.XFILE; end
YCAP  = [];           if isfield(OVR_,'YCAP'),  YCAP  = OVR_.YCAP;  end
OWN   = strcmp(XFILE, 'x150.mat');
X = load(fullfile(DAT, XFILE));
s = X.s;   y = X.bbS;
if OWN
V = load(fullfile(DAT, 'sph_v150.mat'));
assert(V.RES.L == X.LDEG, 'sph_bb:degree', ...
       'x150.mat was drawn at L = %d but sph_v150.mat selected L = %d', X.LDEG, V.RES.L);
else
V.RES = struct('R',X.R, 'L',X.LDEG, 'K',X.K, 'Nall',X.Nnodes, 'nmae',X.NMAE, ...
               'npass',NaN, 'Ltried',NaN);
end

fprintf(['solid harmonics, R <= %g um: L = %d (K = %d)' newline], V.RES.R, V.RES.L, V.RES.K);
fprintf(['  NMAE over the %d nodes inside R <= 150 um, six excitations = %.4f %%' newline], ...
        V.RES.Nall, V.RES.nmae);
fprintf(['  %d of %d degrees tried passed the gate; selection = highest passer' newline], ...
        V.RES.npass, V.RES.Ltried);
fprintf(['  b.b on x_a: %.4f .. %.4f mT^2   (at s = 0  %.4f)' newline], ...
        min(y), max(y), interp1(s, y, 0));

% [ADDED 2026-09-28] MODE picks what the figure shows.
%   'bb'    one panel, b.b on the x_a axis                (the original figure)
%   'grad3' three stacked panels, d(b.b)/ds on each of the three actuator axes,
%           which is what the smoothness gate looks at. The three columns of
%           X.GS carry their own excitation, X.POLE = [1 3 6]: x_a is P1, y_a is
%           P3, z_a is P6 -- the same pairing sph_strat uses.
% Rule 8: in the stack only the lowest panel carries the axis title and the tick
% numbers; the ones above it draw the ticks and nothing else.
MODE = 'grad3';

FS = 60;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;  LW = 7;
COL = [0.05 0.10 0.95];                            % blue: the harmonic model

ax_nm = {'x_a','y_a','z_a'};
for q = 1:3
    fprintf(['  d(b.b)/d%s (P%d): %.4f .. %.4f mT^2/um   (at s = 0  %.4f)   negatives %d' newline], ...
            ax_nm{q}, X.POLE(q), min(X.GS(:,q)), max(X.GS(:,q)), ...
            interp1(s, X.GS(:,q), 0), nnz(X.GS(:,q) < 0));
end

[XR, XT] = axis_sym_(X.SPAN);
switch lower(MODE)
    case 'grad3'
        YY  = {X.GS(:,1), X.GS(:,2), X.GS(:,3)};
        % No unit in the label: with (mT^2/um) appended the string is taller than
        % a 0.225-high panel and the three run into one another. All three are
        % mT^2/um; that belongs in the caption.
        YL  = {'$\mathbf{d(b\cdot b)/dx_a}$', '$\mathbf{d(b\cdot b)/dy_a}$', ...
               '$\mathbf{d(b\cdot b)/dz_a}$'};
        % raised off the bottom: at 0.115 the axis title was cut by the canvas
        POS = {[0.2150 0.6850 0.7000 0.2250], ...
               [0.2150 0.4200 0.7000 0.2250], ...
               [0.2150 0.1550 0.7000 0.2250]};
        out = fullfile(FIG, 'sph_du3_R150.png');
        if ~OWN
            ct = '';   if ~isempty(YCAP), ct = '_zoom'; end
            out = fullfile(FIG, sprintf('sph_du3_R%d_L%d%s.png', round(X.R), X.LDEG, ct));
        end
    case 'bb'
        YY  = {y};   YL = {'$\mathbf{b\cdot b\;(mT^{2})}$'};
        POS = {[0.1895 0.1367 0.7155 0.7886]};
        out = fullfile(FIG, 'sph_bb_R150.png');
    otherwise
        error('sph_bb:mode', 'MODE is ''bb'' or ''grad3''');
end

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
AX  = gobjects(1,numel(YY));
for q = 1:numel(YY)
    ax = axes('Parent',fig,'Position',POS{q});  hold(ax,'on');   AX(q) = ax;
    try, ax.Toolbar = []; end                                        %#ok<TRYNC>
    try, ax.Interactions = []; end                                   %#ok<TRYNC>

    plot(ax, s, YY{q}, '-', 'Color', COL, 'LineWidth', LW);

    if ~isempty(YCAP)
        YR = YCAP;   YT = YR(1) + (1:3)*diff(YR)/4;
    elseif min(YY{q}) < 0
        [YR, YT] = ylim_any_(min(YY{q}), max(YY{q}));
    else
        [YR, YT] = ylim_odd_(max(YY{q}));
    end
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
           'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
           'XTickLabelRotation',0,'YTickLabelRotation',0);
    xlim(ax, XR);   set(ax,'XTick', XT);
    ylim(ax, YR);   set(ax,'YTick', YT);
    ylabel(ax, YL{q}, 'Interpreter','latex','FontSize',FSLAB);

    if q < numel(YY)                               % upper panel: ticks only
        set(ax, 'XTickLabel', []);
    else                                           % lower panel: the shared x axis
        for xv = XR                                % rule 4: x ends labelled by hand
            text(ax, xv, YR(1)-0.022*diff(YR), sprintf('%g', xv), ...
                 'HorizontalAlignment','center','VerticalAlignment','top', ...
                 'FontSize',FS,'FontWeight','bold','Clipping','off');
        end
        % 'grad3' walks a different axis in each panel, so the shared abscissa is
        % the along-axis coordinate s, not x_a.
        xlb = '$\mathbf{x_a\;(\mu m)}$';
        if strcmpi(MODE,'grad3'), xlb = '$\mathbf{s\;(\mu m)}$'; end
        xlabel(ax, xlb, 'Interpreter','latex','FontSize',FSLAB);
    end
    hold(ax,'off');
end

drawnow;                                           % keep a long y label off the edge
need = max(arrayfun(@(a) a.TightInset(1), AX)) + 0.010;
for q = 1:numel(AX)
    if AX(q).Position(1) < need
        rt = AX(q).Position(1) + AX(q).Position(3);
        AX(q).Position(1) = need;   AX(q).Position(3) = rt - need;
    end
end
set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');   close(fig);
fprintf(['wrote %s' newline], out);

% ---- local helpers ------------------------------------------------------------
function [lim, tk] = axis_sym_(R)
    s = R/2;   lim = [-R R];   tk = [-s 0 s];
end

function [lim, tk] = ylim_any_(lo, hi)
% three ticks, end gaps equal to the spacing, for data that goes below zero
    cand = sort(reshape(([1 1.5 2 2.5 3 4 5 6 8].') .* 10.^(-4:6), [], 1));
    for st = cand.'
        a = floor(lo/st)*st;
        if a + 4*st >= hi, lim = [a, a + 4*st];   tk = a + (1:3)*st;   return, end
    end
    error('ylim_any_:none', 'no step for [%g %g]', lo, hi);
end

function [lim, tk] = ylim_odd_(maxv)
    cand = [1 1.1 1.2 1.25 1.5 1.6 1.75 2 2.25 2.5 3 4 5 6 7.5 8];
    nice = [1 1.5 2 2.5 3 4 5];   need = maxv/4;
    for kk = floor(log10(need)) : floor(log10(need))+1
        ok = cand(cand*10^kk >= need*(1-1e-12)) * 10^kk;
        if isempty(ok), continue, end
        s = min(ok);
        okn = nice(nice*10^kk >= need*(1-1e-12)) * 10^kk;
        if ~isempty(okn) && min(okn) <= 1.2*s, s = min(okn); end
        lim = [0 4*s];   tk = (1:3)*s;   return
    end
    error('ylim_odd_:none','no step for %g', maxv);
end
