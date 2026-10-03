% xa_grad  d(b.b)/dx_a on the x_a axis: RBF at three data ranges, plus the
%          18-parameter charge model from the flux calibration.
%
%   Four curves, one unit (mT^2/um), P1 excited at 1 A, drawn over |s| <= 150 um
%   so that every model stays inside its own fitted range.
%
%       METHOD = 'rbf' | 'sph' picks which run to draw.  Four curves at
%       R = 150 / 250 / 350 / 450 um, each with that range's winning setting from
%       range_scan.m, fitted on ALL nodes inside that radius.
%       Charge model                 grad(b.b) = (gB^2/l_hat) g' L_k g, 18 par.
%
%   The eighteen-parameter curve is dashed and black: it is the only one not fitted
%   to the FEM nodes directly, so it reads as the reference the RBF curves are being
%   compared against.
%
%   This is the GRADIENT, not the force; the force is 22.55 times it (0.5*mgB*UF).
%
%   Reads only (plot-scripts-pure); xa_grad.mat is written by temp_code/xa_grad.m.
%
% Output: figures/long2016_hexapole_halfcut/xa_grad.png

HERE = fileparts(mfilename('fullpath'));
FMX  = fileparts(fileparts(HERE));
DAT  = fullfile(FMX,'utils','data');
FIG  = fullfile(FMX,'figures','long2016_hexapole_halfcut','current');

METHOD = 'rbf';                 % 'rbf' | 'sph' -- which run to draw
S  = load(fullfile(DAT, sprintf('xa_grad_%s.mat', METHOD)));
CV = S.CV;   s = S.s;
% the flux-derived curve is labelled "Eighteen parameters" (user's call); accept
% .mat files still carrying the earlier "Charge model" wording
for q = 1:numel(CV)
    if contains(CV(q).name,'Charge'), CV(q).name = 'Eighteen parameters';  end
end
for q = 1:numel(CV)
    fprintf(['%-24s %.4f .. %.4f mT^2/um' newline], CV(q).name, min(CV(q).g), max(CV(q).g));
end

FS = 60;  FSLEG = 45;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;  LW = 7;
COL = {[0.85 0.10 0.10], [0.00 0.60 0.25], [0.05 0.10 0.95], [0.55 0.20 0.75], [0 0 0]};
STY = {'-','-','-','-','--'};
LWS = {LW, LW, LW, LW, LW-1};

[XR, XT] = axis_sym_(S.SPAN);
gall = cell2mat(arrayfun(@(c) c.g(:), CV, 'UniformOutput', false));
[YR, YT] = ylim_odd_(max(gall(:)));

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]);  hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end

h = gobjects(1,numel(CV));
for q = 1:numel(CV)
    h(q) = plot(ax, s, CV(q).g, STY{q}, 'Color', COL{q}, 'LineWidth', LWS{q});
end

set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
       'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
       'XTickLabelRotation',0,'YTickLabelRotation',0);
xlim(ax, XR);   set(ax,'XTick', XT);
ylim(ax, YR);   set(ax,'YTick', YT);
for xv = XR                                       % rule 4: x ends labelled by hand
    text(ax, xv, YR(1)-0.022*diff(YR), sprintf('%g', xv), ...
         'HorizontalAlignment','center','VerticalAlignment','top', ...
         'FontSize',FS,'FontWeight','bold','Clipping','off');
end
xlabel(ax, '$\mathbf{x_a\;(\mu m)}$',                        'Interpreter','latex','FontSize',FSLAB);
ylabel(ax, '$\mathbf{d(b\cdot b)/dx_a\;(mT^{2}/\mu m)}$',     'Interpreter','latex','FontSize',FSLAB);
lg = legend(ax, h, {CV.name}, 'Interpreter','tex', 'Location','northwest', 'NumColumns',1);
lg.FontSize = FSLEG;  lg.FontWeight = 'bold';
lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.Color = 'w';
lg.ItemTokenSize = [55 25];
hold(ax,'off');

out = fullfile(FIG, sprintf('xa_grad_%s.png', METHOD));
set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');   close(fig);
fprintf(['wrote %s' newline], out);

% ---- local helpers ------------------------------------------------------------
function [lim, tk] = axis_sym_(R)
    s = R/2;   lim = [-R R];   tk = [-s 0 s];
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
