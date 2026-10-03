% sph_range  d(b.b)/dx_a on the x_a axis, one curve per data range R = 150:10:500.
%
%   Each curve is a solid-harmonic fit on ALL nodes inside its own R, at the degree
%   that (1) passes the smoothness gate over |s| <= R and (2) has the lowest NMAE on
%   the fixed R <= 150 working region.  Colour encodes R.
%
%   The bundle stays together from R = 150 to about 440 and then peels away: the gate
%   forces the degree down to L = 4 (R >= 450) and L = 2 (R >= 480), and those low
%   orders cannot hold the curvature.  The stragglers therefore measure the gate, not
%   the harmonic basis.
%
%   P1 excited at 1 A, |x_a| <= 150 um.  Gradient, not force (force = 22.55 x this).
%
%   Reads only (plot-scripts-pure); sph_range.mat is written by temp_code/sph_range.m.
%
% Output: figures/long2016_hexapole_halfcut/sph_range.png

HERE = fileparts(mfilename('fullpath'));
FMX  = fileparts(fileparts(HERE));
DAT  = fullfile(FMX,'utils','data');
FIG  = fullfile(FMX,'figures','long2016_hexapole_halfcut','current');

S   = load(fullfile(DAT,'sph_range.mat'));
RES = S.RES;   s = S.s;
R   = [RES.R].';   nR = numel(R);
fprintf(['%d ranges, R = %d .. %d um | NMAE(R<=%d) %.4f .. %.4f %%' newline], ...
        nR, min(R), max(R), S.REVAL, min([RES.nmae]), max([RES.nmae]));

FS = 60;  FSLEG = 45;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;  LW = 7;
% [MODIFIED 2026-09-17] a legend replaces the colour bar (user's call).  36 entries
%   will not fit at 45 pt, so only every fifth range is drawn -- 150:50:500 still
%   spans all three regimes the sweep found (L = 6..9, then L = 4, then L = 2).
%   RSHOW / OUTTAG are the per-run knobs: which ranges to draw, and the suffix that
%   keeps each selection in its own file.  '' + 150:50:500 is the overview figure.
RSHOW  = [470 480];
OUTTAG = '_470_480';
ish    = arrayfun(@(rr) find(R == rr, 1), RSHOW);
if numel(ish) <= 3, CM = [0.05 0.10 0.95; 0.85 0.10 0.10; 0 0 0];  CM = CM(1:numel(ish),:);
else,               CM = turbo(numel(ish));  end

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]);  hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end

h = gobjects(1,numel(ish));   lab = cell(1,numel(ish));
for q = 1:numel(ish)
    k = ish(q);
    h(q)   = plot(ax, s, RES(k).g, '-', 'Color', CM(q,:), 'LineWidth', LW);
    lab{q} = sprintf('R = %d, L = %d  (NMAE %.3f%%)', RES(k).R, RES(k).L, RES(k).nmae);
end

[XR, XT] = axis_sym_(S.SPAN);
gall = cell2mat(arrayfun(@(r) r.g(:), RES(ish), 'UniformOutput', false));
[YR, YT] = ylim_odd_(max(gall(:)));
set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
       'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
       'XTickLabelRotation',0,'YTickLabelRotation',0);
xlim(ax, XR);   set(ax,'XTick', XT);
ylim(ax, YR);   set(ax,'YTick', YT);
for xv = XR                                       % rule 4: x ends laid by hand
    text(ax, xv, YR(1)-0.022*diff(YR), sprintf('%g', xv), ...
         'HorizontalAlignment','center','VerticalAlignment','top', ...
         'FontSize',FS,'FontWeight','bold','Clipping','off');
end
xlabel(ax, '$\mathbf{x_a\;(\mu m)}$',                    'Interpreter','latex','FontSize',FSLAB);
ylabel(ax, '$\mathbf{d(b\cdot b)/dx_a\;(mT^{2}/\mu m)}$','Interpreter','latex','FontSize',FSLAB);

lg = legend(ax, h, lab, 'Interpreter','tex', 'Location','northwest', 'NumColumns',1);
lg.FontSize = FSLEG*0.78;  lg.FontWeight = 'bold';
lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.Color = 'w';
lg.ItemTokenSize = [55 25];
hold(ax,'off');

out = fullfile(FIG, ['sph_range' OUTTAG '.png']);
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
