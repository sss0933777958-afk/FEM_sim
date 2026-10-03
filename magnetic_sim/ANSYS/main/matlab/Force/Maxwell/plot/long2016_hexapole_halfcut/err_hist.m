%% err_hist.m -- full-node force residual, single vs eighteen, as a histogram.
%  Reads the two calibration records written by main/main.m (single, eighteen),
%  which carry the full-node residual as res_eval, and draws one figure.
%  No fitting, no field reading: load -> pick -> plot (plot-scripts-pure).
%
%  Style: the nine figure rules. Tick numbers 60, legend 45, box and legend box
%  5.0, the horizontal axis numbers its endpoints while the vertical axis does
%  not, three equally spaced ticks on both axes with the end gaps equal to the
%  spacing, horizontal tick step a multiple of 0.1, square canvas written with
%  print + PaperPosition. Tick numbers and legend are tex + bold (the latex
%  interpreter ignores FontWeight); only the axis titles are latex.
clc;
if exist('PLOT_OVR','var'), OVR_ = PLOT_OVR; else, OVR_ = struct(); end
clearvars -except OVR_;
BASE  = 'current';     % 'current' | 'voltage'  (PLOT_OVR.BASE overrides)
if isfield(OVR_,'BASE'), BASE = OVR_.BASE; end
SOFF  = 3;             % sensor distance [mm], voltage base only  (PLOT_OVR.SOFF overrides)
if isfield(OVR_,'SOFF'), SOFF = OVR_.SOFF; end
% [ADDED 2026-09-29] What the two histograms are (PLOT_OVR.CMP overrides):
%   'se'     single vs eighteen, each over its whole evaluation set
%   'shell'  ONE record (eighteen), its full-node residual split by node radius
%            into r <= RSPLIT and RSPLIT < r <= RSEL. Same model, same
%            evaluation, only grouped by where the node sits.
CMP    = 'se';
RSPLIT = 150;          % [um], CMP = 'shell' only
RSEL  = 150;
EXC   = 'singles6';
BINW  = 0.02;          % shared bin width [pN]; shared so the two heights compare
if isfield(OVR_,'CMP'),  CMP  = OVR_.CMP;  end
if isfield(OVR_,'RSEL'), RSEL = OVR_.RSEL; end
if isfield(OVR_,'BINW'), BINW = OVR_.BINW; end
NTK   = 3;             % ticks per axis
MODEL = 'long2016_hexapole_halfcut';
here  = fileparts(mfilename('fullpath'));
FMX   = fileparts(fileparts(here));                       % .../Force/Maxwell
et    = '';   if ~strcmpi(EXC,'pairs21'), et = ['_' lower(EXC)]; end
MATD  = fullfile(FMX,'data',MODEL,'.mat');
st = '';   if strcmp(BASE,'voltage'), st = ['_soff' strrep(sprintf('%g',SOFF),'.','p') 'mm']; end
S.eighteen = pick_(MATD, sprintf('%s_R%d_N*%s_L*%s_eighteen.mat', BASE, RSEL, st, et), EXC);
OUT   = fullfile(FMX,'figures',MODEL,BASE);

FS = 60;  FSLEG = 45;  FSLAB = 44;  LWBOX = 5.0;  CANV = 14.5;
BLU = [0.05 0.10 0.95];   RED = [0.85 0.10 0.10];

if strcmp(CMP, 'shell')
    % res_eval is ordered (node i, excitation j) -> (j-1)*Np + i, so the node
    % radius is simply repeated once per excitation.
    E   = S.eighteen;
    rad = repmat(vecnorm(E.P_eval, 2, 2), size(E.u,2), 1);
    assert(numel(rad) == numel(E.res_eval), 'err_hist:shape', 'residual and node list disagree');
    in  = rad <= RSPLIT;
    e0  = E.res_eval(in);    e1 = E.res_eval(~in);
    LEG = {sprintf('R \\leq %d \\mum', RSPLIT), sprintf('%d < R \\leq %d \\mum', RSPLIT, RSEL)};
    ct  = '_shell';
    fprintf('inner r <= %d      : %6d values  NMAE %.4f %%  mean %.4f pN  max %.4f pN%s', RSPLIT, ...
            numel(e0), sum(e0)/sum(E.fn_eval(in))*100,  mean(e0), max(e0), newline);
    fprintf('outer %d < r <= %d : %6d values  NMAE %.4f %%  mean %.4f pN  max %.4f pN%s', RSPLIT, RSEL, ...
            numel(e1), sum(e1)/sum(E.fn_eval(~in))*100, mean(e1), max(e1), newline);
else
    S.single = pick_(MATD, sprintf('%s_R%d_N*%s_L*%s_single.mat', BASE, RSEL, st, et), EXC);
    e0 = S.single.res_eval;   e1 = S.eighteen.res_eval;
    LEG = {'Single parameter','Eighteen parameters'};
    ct  = '';
    fprintf('single   : NMAE %.4f %%  mean %.4f pN  max %.4f pN%s', S.single.NMAE,   mean(e0), max(e0), newline);
    fprintf('eighteen : NMAE %.4f %%  mean %.4f pN  max %.4f pN%s', S.eighteen.NMAE, mean(e1), max(e1), newline);
end
fprintf('bin width %.3g pN%s', BINW, newline);

% ---- horizontal axis: [0, (n+1)s], s the smallest multiple of 0.1 that covers
SHELL = strcmp(CMP, 'shell');
hi = max([e0; e1]);
if SHELL
    % Same window rule as paper_fig_plot/plot_err_hist_shell.m: the frame is set
    % by the 99.5th percentile, so a thin tail does not flatten the figure.
    hi = 1.02 * prctile([e0; e1], 99.5);
    fprintf('x window set by the 99.5th percentile (%.3f pN); %.2f %% of the values lie beyond%s', ...
            hi/1.02, 100*mean([e0; e1] > hi), newline);
end
sx = ceil(hi/(NTK+1)/0.1 - 1e-9) * 0.1;
xr = [0, (NTK+1)*sx];   xt = round((1:NTK)*sx*10)/10;

edg = 0 : BINW : (ceil(max([e0; e1; xr(2)])/BINW)*BINW);
ctr = (edg(1:end-1) + edg(2:end))/2;
if SHELL
    % STACKED, ONE DENOMINATOR (user's call, 2026-09-29; the style of
    % err_hist_conv_maxwell_eighteen_R300.png). Every bar is the share of ALL
    % values in that bin, split into its inner and outer part, so the two
    % colours never overlap and the outline is the unsplit distribution.
    nall = numel(e0) + numel(e1);
    p0   = histcounts(e0, edg) / nall * 100;
    p1   = histcounts(e1, edg) / nall * 100;
    ptop = p0 + p1;
else
    p0  = histcounts(e0, edg) / numel(e0) * 100;
    p1  = histcounts(e1, edg) / numel(e1) * 100;
    ptop = max(p0, p1);
end

% ---- vertical axis: same construction, a nice step if it wastes <= 20 %
smin = 1.02*max(ptop)/(NTK+1);
nice = sort(reshape(([1 1.5 2 2.5 3 4 5].') .* 10.^(-3:3), [], 1));
sy   = min(nice(nice >= smin));
if sy > 1.2*smin
    d1 = 10^floor(log10(smin));                           % one significant digit first
    sy = ceil(smin/d1 - 1e-9)*d1;
    if sy > 1.2*smin
        dg = d1/10;                                       % then two
        sy = ceil(smin/dg - 1e-9)*dg;
    end
end
yr = [0, (NTK+1)*sy];   yt = (1:NTK)*sy;

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes(fig, 'Units','normalized', 'Position',[0.215 0.135 0.700 0.760]);
try, ax.Toolbar = []; end                                            %#ok<TRYNC>
try, ax.Interactions = []; end                                       %#ok<TRYNC>
hold(ax,'on');  box(ax,'on');
TDIR = 'in';
if SHELL
    hb = bar(ax, ctr, [p0(:) p1(:)], 1, 'stacked', 'EdgeColor','none');
    hb(1).FaceColor = RED;   hb(1).FaceAlpha = 0.60;      % inner, as in the reference
    hb(2).FaceColor = BLU;   hb(2).FaceAlpha = 0.60;      % outer
    h0 = hb(1);   h1 = hb(2);   TDIR = 'out';
else
    h0 = bar(ax, ctr, p0, 1, 'FaceColor',BLU, 'FaceAlpha',0.60, 'EdgeColor','none');
    h1 = bar(ax, ctr, p1, 1, 'FaceColor',RED, 'FaceAlpha',0.60, 'EdgeColor','none');
end

set(ax, 'XLim',xr, 'YLim',yr, 'XTick',xt, 'YTick',yt, ...
        'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, ...
        'TickDir',TDIR, 'TickLength',[0.018 0.018], 'Layer','top');
xlabel(ax, '$\mathbf{Residual\;(pN)}$',   'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{Percentage\;(\%)}$', 'Interpreter','latex', 'FontSize',FSLAB);

% rule 4: the horizontal endpoints are numbered by hand, they are not ticks
yoff = yr(1) - 0.022*diff(yr);
if strcmp(TDIR, 'out')
    % outward ticks push the tick numbers down by one tick length; the endpoint
    % numbers follow so all five sit on one baseline, and the frame is raised to
    % keep the axis title on the canvas
    yoff = yr(1) - (0.022 + 0.018)*diff(yr);
    ax.Position = [0.215 0.160 0.700 0.735];
end
for xv = xr
    text(ax, xv, yoff, sprintf('%g',xv), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end

lg = legend(ax, [h0 h1], LEG, ...
            'Interpreter','tex', 'FontSize',FSLEG, 'FontWeight','bold', ...
            'Location','northeast', 'NumColumns',1);
lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.Color = 'w';
lg.ItemTokenSize = [55 25];

bt   = '';   if strcmp(BASE,'voltage'), bt = '_V'; end      % current stays untagged
name = fullfile(OUT, sprintf('err_hist%s%s_R%d%s%s.png', bt, ct, RSEL, st, et));
set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
print(fig, name, '-dpng', '-r200');
close(fig);
fprintf('saved %s%s', name, newline);

% ---- the one calibration record matching a pattern --------------------------
function r = pick_(d, pat, exc)
    L = dir(fullfile(d, pat));
    r = [];
    for k = 1:numel(L)
        t = load(fullfile(d, L(k).name));
        if isfield(t,'EXC') && strcmpi(t.EXC, exc) && isfield(t,'res_eval')
            assert(isempty(r), 'err_hist:ambiguous', 'more than one record matches %s', pat);
            r = t;
        end
    end
    assert(~isempty(r), 'err_hist:missing', ...
           'no record with a full-node evaluation matches %s (re-run main.m)', pat);
end
