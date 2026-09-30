%% conv_lg.m -- l_hat and F g_I against the point count, flux vs force.
%  Reads conv_R<R>_<tag>.mat written by utils/conv_series.m and draws two
%  figures. No fitting, no field reading: load -> pick -> plot (plot-scripts-pure).
%
%  Style: the nine figure rules. Tick numbers 60, legend 45, box and legend box
%  5.0, horizontal axis labels its endpoints while the vertical axis does not,
%  odd equally spaced ticks with the end gaps equal to the spacing, square
%  canvas via print + PaperPosition (exportgraphics trims the white border and
%  would break the equal sides).
clearvars; clc;

% [ADDED 2026-09-28] What is compared:
%   'ff'  flux vs force, both eighteen      (conv_R<R>_<tag>.mat, conv_series)
%   'se'  force single vs force eighteen    (the two calibration records of
%         main.m, which carry the whole ladder as rec.ladder / rec.SER)
CMP   = 'se';
EXC   = 'singles6';     % CMP='se' only: which excitation set the records used
NEXT  = 2;              % CMP='se' only: ladder levels drawn past the later N_c
TAG   = 'eighteen';
RSEL  = 150;
% Convergence point of the force ladder (Nr=10 -> 61 points), the later of the
% two; the flux ladder stops at 13. The axis ends there, as asked.
NC    = 61;
% Value each curve carries at zero points: l_hat starts from the design length
% l0 = 500 um, the gain from nothing.
L0    = 500;   G0 = 0;
here  = fileparts(mfilename('fullpath'));
FMX   = fileparts(fileparts(here));                       % .../Force/Maxwell
MODEL = 'long2016_hexapole_halfcut';
MATD  = fullfile(FMX,'data',MODEL,'.mat');
if strcmp(CMP, 'ff')
    S   = load(fullfile(MATD, sprintf('conv_R%d_%s.mat',RSEL,TAG)));
    LEG = {'Flux','Force'};   SFX = '';
else
    % Series A = single, series B = eighteen, stored in the l_flux / l_force
    % slots so the drawing loop below is shared. B g is recovered from the
    % stored F g by the gauge relation (a unit conversion, not a refit).
    et = '';   if ~strcmpi(EXC,'pairs21'), et = ['_' lower(EXC)]; end
    ra = pick_(MATD, sprintf('current_R%d_N*_L*%s_single.mat',   RSEL, et), EXC);
    rb = pick_(MATD, sprintf('current_R%d_N*_L*%s_eighteen.mat', RSEL, et), EXC);
    NC = max(ra.npts, rb.npts) + 6*NEXT;       % 6 points per ladder level
    n  = min(size(ra.SER,1), size(rb.SER,1));
    gB = @(r) sqrt(r.SER(1:n,2) .* 2 .* r.SER(1:n,1) ./ (r.UF * r.mgB));
    S  = struct('npts',6*ra.ladder(1:n)+1, ...
                'l_flux',ra.SER(1:n,1), 'l_force',rb.SER(1:n,1), ...
                'g_flux',gB(ra),        'g_force',gB(rb));
    LEG = {'1 parameter','18 parameters'};   SFX = ['_se' et];
    fprintf('N_c: single %d, eighteen %d -> axis ends at %d%s', ...
            ra.npts, rb.npts, NC, newline);
end
OUT   = fullfile(FMX,'figures',MODEL);
if ~exist(OUT,'dir'), mkdir(OUT); end

FS    = 60;        % tick numbers          (rule 1)
FSLEG = 45;        % legend                (rule 2)
FSLAB = 44;        % axis titles
LWBOX = 5.0;       % axes box and legend box (rules 3, 7)
LW    = 6.0;       % data lines            (rule 2)
MS    = 16;        % markers
CANV  = 14.5;      % square canvas [in]    (rule 6)
NTK   = 3;         % odd ticks on both axes (rule 5)
BLU   = [0 0 1];   RED = [0.85 0 0];

keep = S.npts <= NC;                       % the ladder up to the convergence point
x    = [0; S.npts(keep)];

for q = 1:2
    if q == 1
        yF = [L0; S.l_flux(keep)];    yC = [L0; S.l_force(keep)];
        ylab = '$\hat{\ell}~(\mu\mathrm{m})$';   name = sprintf('conv_ell%s_R%d.png', SFX, RSEL);
    else
        yF = [G0; S.g_flux(keep)];    yC = [G0; S.g_force(keep)];
        ylab = '$^{B}\hat{g}_{I}~(\mathrm{mT/A})$';   name = sprintf('conv_g%s_R%d.png', SFX, RSEL);
    end

    % The axis ends at the convergence point (user's call), so the end gaps are
    % not forced to equal the tick spacing here; the ticks stay odd and evenly
    % spaced inside.
    xr = [0 NC];
    xs = nicestep_((NC - 0) / (NTK + 1));
    xt = (1:NTK) * xs;
    lo = min([yF(:); yC(:)]);   hi = max([yF(:); yC(:)]);
    [yr, yt] = odd_axis(lo, hi, NTK);

    % The right edge stops short of the canvas: the last endpoint number is
    % centred on it and would otherwise be cut in half (measured: it was).
    AXPOS = [0.215 0.135 0.700 0.700];
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes(fig, 'Units','normalized', 'Position',AXPOS);
    try, ax.Toolbar = []; end                                        %#ok<TRYNC>
    try, ax.Interactions = []; end                                   %#ok<TRYNC>
    hold(ax,'on');  box(ax,'on');

    % markers thinned to about eight so they do not run into one another
    mi = unique([1 round(linspace(1, numel(x), 8))]);
    hF = plot(ax, x, yF, '-o', 'Color',BLU, 'LineWidth',LW, ...
              'MarkerSize',MS, 'MarkerFaceColor',BLU, 'MarkerIndices',mi);
    hC = plot(ax, x, yC, '-s', 'Color',RED, 'LineWidth',LW, ...
              'MarkerSize',MS, 'MarkerFaceColor','w', 'MarkerIndices',mi);

    set(ax, 'XLim',xr, 'YLim',yr, 'XTick',xt, 'YTick',yt, ...
            'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, ...
            'TickDir','in', 'TickLength',[0.018 0.018]);
    xlabel(ax, '\textbf{Number of points}', 'Interpreter','latex', 'FontSize',FSLAB);
    ylabel(ax, ylab, 'Interpreter','latex', 'FontSize',FSLAB);

    % rule 4: the horizontal axis labels its endpoints, and they are drawn by
    % hand because they are deliberately not ticks.
    yoff = yr(1) - 0.022*diff(yr);
    for xv = xr
        text(ax, xv, yoff, sprintf('%g',xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', ...
             'Clipping','off');
    end

    lg = legend(ax, [hF hC], LEG, 'Interpreter','tex', ...
                'FontSize',FSLEG, 'FontWeight','bold', 'Location','northoutside', ...
                'Orientation','horizontal');
    lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;
    lg.ItemTokenSize = [55 25];
    ax.Position = [0.215 0.135 0.745 0.700];      % northoutside shrinks it back

    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
    print(fig, fullfile(OUT,name), '-dpng', '-r200');
    close(fig);
    fprintf('saved %s\n', fullfile(OUT,name));
end

% ---- the one calibration record matching a pattern --------------------------
function r = pick_(d, pat, exc)
    L = dir(fullfile(d, pat));
    r = [];
    for k = 1:numel(L)
        t = load(fullfile(d, L(k).name));
        if isfield(t,'EXC') && strcmpi(t.EXC, exc) && isfield(t,'SER')
            assert(isempty(r), 'conv_lg:ambiguous', 'more than one record matches %s', pat);
            r = t;
        end
    end
    assert(~isempty(r), 'conv_lg:missing', 'no record matches %s', pat);
end

% ---- odd, equally spaced ticks whose end gaps equal the spacing -------------
function [lim, tk] = odd_axis(lo, hi, n)
    cands = reshape(([1 1.5 2 2.5 3 4 5 6 8].') .* 10.^(-4:6), [], 1);
    cands = sort(cands);
    for s = cands.'
        a = floor(lo/s)*s;                      % a <= lo, a is a multiple of s
        if a + (n+1)*s >= hi
            tk  = a + (1:n)*s;
            lim = [a, a + (n+1)*s];
            return
        end
    end
    error('odd_axis:range', 'no nice step covers [%g %g]', lo, hi);
end

% ---- largest nice step not exceeding t --------------------------------------
function s = nicestep_(t)
    cands = reshape(([1 1.5 2 2.5 3 4 5 6 8].') .* 10.^(-4:6), [], 1);
    cands = sort(cands);
    s = max(cands(cands <= t));
end
