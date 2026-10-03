% rbf_f18  Force along each actuator axis: eighteen-parameter force model vs RBF.
%
%   One figure per axis, the two models overlaid:
%       x_a  <- P1 excited      y_a  <- P3 excited      z_a  <- P6 excited, drawn as -F_za
%   (P6 sits on -z_a, so its force along +z_a is negative; both curves are flipped
%   together, which changes neither the gap between them nor any sign test.)
%
%   Eighteen-parameter force model -- built from the FLUX eighteen calibration
%   (l_hat + 17 e, with B g_I and K_I bar) done over r <= 300 um:
%       f = base_model(build_L(P, l, e, Pc_base), K_I bar, Fmap, F g_I)
%       F g_I = UF * mgB/(2 l) * (B g_I)^2
%   RBF -- rbf_field(500, struct('rho',240,'lam',1e-4)),  F = 0.5*mgB*UF*d(b.b)/ds.
%   Both are continuous functions of s; they are evaluated densely enough to draw.
%
%   Reads only (plot-scripts-pure).  The .mat was produced in the MATLAB session:
%       cal = load_flux_calib('long2016_hexapole_halfcut','tip40um', ...
%                             'calib_current_maxwell_axshN97_R300_eighteen','current');
%       [~, g] = rbf_field(500, struct('rho',240,'lam',1e-4));
%       for a = 1:3, Pq(:,a) = s;  F18(:,a) = f(a,:,POLE(a));  dU = g(Pq,I);  Frbf(:,a) = ...
%   with mgB = 0.0451, UF = 1e3, POLE = [1 3 6], SG = [1 1 -1].
%
% Output: figures/long2016_hexapole_halfcut/rbf_f18_{xa,ya,za}.png
%         figures/long2016_hexapole_halfcut/rbf_f18_err.png   (|F_RBF - F_18|, three panels stacked)

HERE = fileparts(mfilename('fullpath'));                      % .../plot/<model>
FMX  = fileparts(fileparts(HERE));                            % .../Force/Maxwell
DAT  = fullfile(FMX,'utils','data');
FIG  = fullfile(FMX,'figures','long2016_hexapole_halfcut','current');

MATF = 'rbf_f18_R500';
if exist('MATF_OVR','var') && ~isempty(MATF_OVR), MATF = MATF_OVR; end
S = load(fullfile(DAT, [MATF '.mat']));

fprintf(['eighteen: %s (l_hat %.1f um, B g_I %.4f mT/A)' newline], S.MATNAME, S.l_hat, S.gI);
fprintf(['RBF     : R = %g um, rho = %g um, lam = %.0e' newline newline], S.RR, S.RHO, S.LAM);
% (%) in brackets = sum|dF| / sum|F18| * 100 over that stretch, dF = F_RBF - F_18
fprintf('%5s %22s %22s %16s\n', 'axis', 'max|dF| (rel) |s|<=150', 'max|dF| (rel) |s|<=300', 'RBF d3 changes');
for a = 1:3
    fprintf('%5s %11.4f (%5.2f %%) %11.4f (%5.2f %%) %9d / %d\n', S.AXN{a}, ...
            S.dmax150(a), S.rel150(a), S.dmax300(a), S.rel300(a), S.n3_150(a), S.n3_300(a));
end

% ================================ figures =====================================
FS = 60;  FSLEG = 45;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;  LW = 6;
cM = [0.85 0.10 0.10];   cR = [0.05 0.10 0.95];               % red = eighteen, blue = RBF
XR = [-S.SMAX S.SMAX];   sx = diff(XR)/4;

for a = 1:3
    AN = S.AXN{a};
    y1 = S.F18(:,a);   y2 = S.Frbf(:,a);

    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]); hold(ax,'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    h1 = plot(ax, S.s, y1, '-',  'Color',cM, 'LineWidth',LW);
    h2 = plot(ax, S.s, y2, '-',  'Color',cR, 'LineWidth',LW);   % [MODIFIED 2026-09-14] solid, per user

    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
           'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
           'XTickLabelRotation',0,'YTickLabelRotation',0);
    xlim(ax, XR);   set(ax,'XTick', XR(1)+(1:3)*sx);
    lo = min([y1; y2]);   hi = max([y1; y2]);
    [yr, yt] = ylim_fit_(lo, hi, 3);
    ylim(ax, yr);   set(ax,'YTick', yt);
    for xv = XR
        text(ax, xv, yr(1)-0.022*diff(yr), sprintf('%g', xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    xlabel(ax, ['$\mathbf{' AN '\;(\mu m)}$'], 'Interpreter','latex', 'FontSize',FSLAB);
    if S.SG(a) < 0
        ylabel(ax, ['$\mathbf{-F_{' AN '}\;(pN)}$'], 'Interpreter','latex', 'FontSize',FSLAB);
    else
        ylabel(ax, ['$\mathbf{F_{' AN '}\;(pN)}$'],  'Interpreter','latex', 'FontSize',FSLAB);
    end

    LOC = 'northwest';   if y1(end) < y1(1), LOC = 'northeast'; end   % legend away from the peak
    lg = legend(ax, [h1 h2], {'Eighteen parameters', 'RBF'}, 'Interpreter','tex', ...
                'Location',LOC, 'NumColumns',1);
    lg.FontSize = FSLEG;  lg.FontWeight = 'bold';
    lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.Color = 'w';
    lg.ItemTokenSize = [55 25];
    hold(ax,'off');

    out = fullfile(FIG, ['rbf_f18_' strrep(AN,'_','') '.png']);
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
    print(fig, out, '-dpng', '-r200');   close(fig);
    fprintf(['wrote %s  (y ticks %s, frame [%g %g])' newline], out, mat2str(yt), yr);
end

% ==================== error stack: |F_RBF - F_18| per axis ====================
% [ADDED 2026-09-14] One figure, three panels stacked top to bottom (x_a, y_a, z_a),
%   sharing the position axis.  Rule 8: the upper two panels keep their tick marks but
%   carry neither the x title nor the tick numbers.  The z_a flip multiplies both curves,
%   so it cancels inside the absolute difference.
CANVE = 16;                                             % square canvas [in] (rule 6)
L0 = 0.20;  W0 = 0.72;  B0 = 0.13;  H0 = 0.25;  G0 = 0.045;   % B0 leaves room for the x title;
                                                               % G0 keeps the outward ticks of
                                                               % neighbouring panels apart
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANVE CANVE]);
for a = 1:3
    AN = S.AXN{a};
    ax = axes('Parent',fig,'Position',[L0, B0+(3-a)*(H0+G0), W0, H0]);  hold(ax,'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    e = abs(S.Frbf(:,a) - S.F18(:,a));
    plot(ax, S.s, e, '-', 'Color',cR, 'LineWidth',LW);
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
           'TickLength',[.012 .012],'TickDir','out','Box','on','Layer','top', ...
           'XTickLabelRotation',0,'YTickLabelRotation',0);
    xlim(ax, XR);   set(ax,'XTick', XR(1)+(1:3)*sx);
    [yr, yt] = ylim_zero_(max(e), 3);
    ylim(ax, yr);   set(ax,'YTick', yt);
    ylabel(ax, ['$\mathbf{|\Delta F_{' AN '}|\;(pN)}$'], 'Interpreter','latex', 'FontSize',FSLAB);
    if a < 3
        set(ax, 'XTickLabel', []);                      % rule 8
    else
        for xv = XR
            text(ax, xv, yr(1)-0.03*diff(yr), sprintf('%g', xv), 'HorizontalAlignment','center', ...
                 'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
        end
        xlabel(ax, '$\mathbf{Position\;(\mu m)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    end
    hold(ax,'off');
    [em, im] = max(e);
    fprintf(['err %s: max %.4f pN at s = %+.1f um | frame [0 %g], ticks %s' newline], ...
            AN, em, S.s(im), yr(2), mat2str(yt));
end
out = fullfile(FIG, 'rbf_f18_err.png');
set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANVE CANVE],'PaperSize',[CANVE CANVE]);
print(fig, out, '-dpng', '-r200');   close(fig);
fprintf(['wrote %s' newline], out);

% ---- local helpers ------------------------------------------------------------------
function [lim, tk] = ylim_zero_(maxv, N)
% Non-negative data: frame [0, (N+1)*s], inner ticks (1:N)*s, so neither end is labelled
% and both end gaps equal the spacing.  s is the smallest of {1 1.5 2 2.5 3 4 5 6 8}x10^k
% that covers maxv.
    if nargin < 2 || isempty(N), N = 3; end
    need = maxv / (N+1);
    for kk = floor(log10(need)) : floor(log10(need))+1
        for m = [1 1.5 2 2.5 3 4 5 6 8]
            s = m * 10^kk;
            if s >= need*(1+1e-9), lim = [0 (N+1)*s];  tk = (1:N)*s;  return, end
        end
    end
end

function [lim, tk] = ylim_fit_(lo, hi, N)
% N equally spaced inner ticks, both frame gaps equal to the tick spacing (rules 4/5),
% and the TIGHTEST such frame a nice step allows.  The step is scanned upwards through
% {1 1.5 2 2.5 3 4 5 6 8}x10^k; every tick must sit on a multiple of the step's own
% unit u, so the labels stay round.  The centre tick is then the first multiple of u
% for which the frame covers [lo, hi].
% (rbf_range_gd.m's ylim_sym_ centres the frame on the data midpoint instead, which
%  left ~40 pN of empty frame under these curves.)
    if nargin < 3 || isempty(N), N = 3; end
    half = (N+1)/2;
    need = (hi - lo) / (N+1);
    for kk = floor(log10(need)) : floor(log10(need))+1
        for m = [1 1.5 2 2.5 3 4 5 6 8]
            s = m * 10^kk;
            if s < need*(1+1e-9), continue; end
            u = 10^floor(log10(s));
            while abs(s/u - round(s/u)) > 1e-9, u = u/2; end
            c = ceil((hi - half*s)/u - 1e-9) * u;          % lowest centre whose top covers hi
            if c - half*s <= lo + 1e-9                      % ... and whose bottom covers lo
                lim = [c - half*s, c + half*s];
                tk  = c + (-(N-1)/2:(N-1)/2)*s;
                return
            end
        end
    end
    error('ylim_fit_:none', 'no nice frame found for [%g %g]', lo, hi);
end
