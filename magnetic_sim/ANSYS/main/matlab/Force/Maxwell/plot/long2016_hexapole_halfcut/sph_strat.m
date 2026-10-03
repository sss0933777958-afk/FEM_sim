% sph_strat  d(b.b)/dx_a on the x_a axis: solid harmonics at four data ranges,
%            stratified sampling, plus the eighteen-parameter model.
%
%   Each range is fitted on 5000 nodes -- the 1771 inside R <= 150 kept in full plus
%   3229 drawn uniformly from outside (R = 150 uses its 1771 alone) -- at the degree
%   that (1) keeps d3(b.b) of one sign on all three actuator axes over |s| <= 150 and
%   (2) has the lowest NMAE on those same 1771 nodes.  The sampling and the criterion
%   are identical to the RBF run, so the two figures can be read side by side.
%
%   P1 excited at 1 A, |x_a| <= 150 um.  Gradient, not force (force = 22.55 x this).
%
%   Reads only (plot-scripts-pure); sph_strat.mat is written by temp_code/sph_strat.m.
%
% Output: figures/long2016_hexapole_halfcut/sph_strat.png

HERE = fileparts(mfilename('fullpath'));
FMX  = fileparts(fileparts(HERE));
DAT  = fullfile(FMX,'utils','data');
FIG  = fullfile(FMX,'figures','long2016_hexapole_halfcut','current');

S   = load(fullfile(DAT,'sph_strat.mat'));
RES = S.RES;   s = S.sA;   g18 = S.g18;

FS = 60;  FSLEG = 45;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;  LW = 7;
COL = {[0.85 0.10 0.10], [0.00 0.60 0.25], [0.05 0.10 0.95], [0.55 0.20 0.75], [0 0 0]};
STY = {'-','-','-','-','--'};
LWS = {LW, LW, LW, LW, LW-1};

Y   = [arrayfun(@(r) {r.g}, RES), {g18}];
% legend carries the range only; the selected L is in the printout and in the .mat
LAB = [arrayfun(@(r) {sprintf('R = %d \\mum', r.R)}, RES), {'Eighteen parameters'}];

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]);  hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end

h = gobjects(1,numel(Y));
for q = 1:numel(Y)
    h(q) = plot(ax, s, Y{q}, STY{q}, 'Color', COL{q}, 'LineWidth', LWS{q});
end

[XR, XT] = axis_sym_(S.RWORK);
[YR, YT] = ylim_odd_(max(cell2mat(cellfun(@(v) v(:), Y, 'UniformOutput', false)), [], 'all'));
set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
       'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
       'XTickLabelRotation',0,'YTickLabelRotation',0);
xlim(ax, XR);   set(ax,'XTick', XT);
ylim(ax, YR);   set(ax,'YTick', YT);
for xv = XR
    text(ax, xv, YR(1)-0.022*diff(YR), sprintf('%g', xv), ...
         'HorizontalAlignment','center','VerticalAlignment','top', ...
         'FontSize',FS,'FontWeight','bold','Clipping','off');
end
xlabel(ax, '$\mathbf{x_a\;(\mu m)}$',                    'Interpreter','latex','FontSize',FSLAB);
ylabel(ax, '$\mathbf{d(b\cdot b)/dx_a\;(mT^{2}/\mu m)}$','Interpreter','latex','FontSize',FSLAB);
lg = legend(ax, h, LAB, 'Interpreter','tex', 'Location','northwest', 'NumColumns',1);
lg.FontSize = FSLEG*0.8;  lg.FontWeight = 'bold';
lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.Color = 'w';
lg.ItemTokenSize = [55 25];
hold(ax,'off');

out = fullfile(FIG,'sph_strat.png');
set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');   close(fig);

fprintf([newline '%5s %5s %6s %14s %12s %14s' newline], 'R','L','K','NMAE(R<=150)%','grad(+150)','rel.err vs 18');
for q = 1:numel(RES)
    gg = RES(q).g;
    fprintf('%5d %5d %6d %14.4f %12.4f %13.2f %%\n', RES(q).R, RES(q).L, RES(q).K, ...
            RES(q).nmae, gg(end), (gg(end)-g18(end))/g18(end)*100);
end
fprintf(['eighteen parameters at x_a=+150: %.4f mT^2/um' newline], g18(end));
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
