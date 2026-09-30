% sph_d3  Third derivative of b.b along each actuator axis, solid harmonics L = 11, R = 250 um.
%
%   One figure per axis.  The curve is d3(b.b)/ds^3 (z_a is flipped with the excitation,
%   as everywhere else).  The stretch carrying the MINORITY sign is drawn in red: that is
%   the part that fails the "sign constant" smoothness test.
%
%   The point of the figure: on all three axes the single sign change sits at
%   |s| = 195 ... 201 um, i.e. always in the same place, and the minority stretch is the
%   outer 10 % of the segment on the side AWAY from the excited pole.  That is the field's
%   own inflection (the opposite, unexcited iron pole pulls |B| back up out there), not a
%   ripple of the fit -- a fitting artefact would move with L and with the axis.
%
%   Reads only (plot-scripts-pure); sph_d3.mat is written by the fitting run.
%
%   The companion figures sph_grad_{xa,ya,za}.png show the GRADIENT itself,
%   d(b.b)/ds, over the same segment and with the same blue/red split, so the stretch
%   that the d3 test objects to can be seen on the curve that actually matters.
%
% Output: figures/long2016_hexapole_halfcut/sph_d3_{xa,ya,za}.png
%         figures/long2016_hexapole_halfcut/sph_grad_{xa,ya,za}.png

HERE = fileparts(mfilename('fullpath'));
FMX  = fileparts(fileparts(HERE));
DAT  = fullfile(FMX, 'data', 'long2016_hexapole_halfcut', '.mat');
FIG  = fullfile(FMX, 'figures', 'long2016_hexapole_halfcut');

S = load(fullfile(DAT, 'sph_d3.mat'));
fprintf(['solid harmonics L = %d, fitted on r <= %d um' newline], S.Lq, S.R0);

FS = 60;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;  LW = 6;
cM = [0.05 0.10 0.95];   cB = [0.85 0.10 0.10];              % majority blue, minority red
XR = [-S.R0 S.R0];   sx = diff(XR)/4;

for a = 1:3
    AN = S.AXN{a};   y = S.D3(:,a);
    g = sign(y);  g(g==0) = 1;
    bad = (g ~= sign(sum(g)));                                % the minority-sign stretch
    fprintf(['%s: sign change at %s um, minority stretch %.1f %% of the segment' newline], ...
            AN, mat2str(round(S.CHG{a},0)), 100*mean(bad));

    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]); hold(ax,'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    yline(ax, 0, 'k--', 'LineWidth', LWBOX*0.7);
    plot(ax, S.s3, y, '-', 'Color',cM, 'LineWidth',LW);
    yb = nan(size(y));  yb(bad) = y(bad);
    plot(ax, S.s3, yb, '-', 'Color',cB, 'LineWidth',LW*1.4);

    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
           'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
           'XTickLabelRotation',0,'YTickLabelRotation',0);
    xlim(ax, XR);   set(ax,'XTick', XR(1)+(1:3)*sx);
    [yr, yt] = ylim_sym_(min(y), max(y), 3);
    ylim(ax, yr);   set(ax,'YTick', yt);
    for xv = XR
        text(ax, xv, yr(1)-0.022*diff(yr), sprintf('%g', xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    xlabel(ax, ['$\mathbf{' AN '\;(\mu m)}$'], 'Interpreter','latex', 'FontSize',FSLAB);
    if S.SG(a) < 0
        ylabel(ax, ['$\mathbf{-\partial^3(b\cdot b)/\partial ' AN '^3}$'], 'Interpreter','latex', 'FontSize',FSLAB);
    else
        ylabel(ax, ['$\mathbf{\partial^3(b\cdot b)/\partial ' AN '^3}$'],  'Interpreter','latex', 'FontSize',FSLAB);
    end
    hold(ax,'off');

    out = fullfile(FIG, ['sph_d3_' strrep(AN,'_','') '.png']);
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
    print(fig, out, '-dpng', '-r200');   close(fig);
    fprintf(['  wrote %s  (y ticks %s)' newline], out, mat2str(yt,3));
end

% ================== the gradient itself, same segment, same split ==================
for a = 1:3
    AN = S.AXN{a};   y = S.FG(:,a);   d3 = S.D3(:,a);
    g = sign(d3);  g(g==0) = 1;   bad = (g ~= sign(sum(g)));
    fprintf(['%s: d(b.b)/ds from %.4f to %.4f mT^2/um; the d3 minority stretch covers ' ...
             's < %g um' newline], AN, y(1), y(end), max(S.s3(bad)));

    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]); hold(ax,'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    plot(ax, S.s3, y, '-', 'Color',cM, 'LineWidth',LW);
    yb = nan(size(y));  yb(bad) = y(bad);
    plot(ax, S.s3, yb, '-', 'Color',cB, 'LineWidth',LW*1.4);

    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
           'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
           'XTickLabelRotation',0,'YTickLabelRotation',0);
    xlim(ax, XR);   set(ax,'XTick', XR(1)+(1:3)*sx);
    [yr, yt] = ylim_sym_(min(y), max(y), 3);
    ylim(ax, yr);   set(ax,'YTick', yt);
    for xv = XR
        text(ax, xv, yr(1)-0.022*diff(yr), sprintf('%g', xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    xlabel(ax, ['$\mathbf{' AN '\;(\mu m)}$'], 'Interpreter','latex', 'FontSize',FSLAB);
    if S.SG(a) < 0
        ylabel(ax, ['$\mathbf{-\partial(b\cdot b)/\partial ' AN '\;(mT^2/\mu m)}$'], 'Interpreter','latex', 'FontSize',FSLAB);
    else
        ylabel(ax, ['$\mathbf{\partial(b\cdot b)/\partial ' AN '\;(mT^2/\mu m)}$'],  'Interpreter','latex', 'FontSize',FSLAB);
    end
    hold(ax,'off');

    out = fullfile(FIG, ['sph_grad_' strrep(AN,'_','') '.png']);
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
    print(fig, out, '-dpng', '-r200');   close(fig);
    fprintf(['  wrote %s' newline], out);
end

% ---- local helper ---------------------------------------------------------------
function [lim, tk] = ylim_sym_(lo, hi, N)
% Signed data: N inner ticks, equal gaps to both frame edges, zero inside the view.
    if nargin < 3 || isempty(N), N = 3; end
    s  = 1.02*(hi-lo)/(N+1);
    u  = 10^floor(log10(s));   s = ceil(s/u - 1e-12)*u;
    for it = 1:8
        c   = round(((lo+hi)/2)/s)*s;
        lim = [c-((N+1)/2)*s, c+((N+1)/2)*s];
        if lim(1) <= lo && lim(2) >= hi, break, end
        s = s + u;
    end
    tk = c + (-(N-1)/2:(N-1)/2)*s;
end
