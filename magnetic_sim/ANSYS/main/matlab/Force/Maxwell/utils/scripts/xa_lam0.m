%% xa_lam0.m -- one unregularised narrow kernel, two data ranges
%
%  Diagnostic (temp_code, so fit and plot live together).
%
%  SAME kernel for both curves: rho = 20 um, lam = 0.  That is about one node
%  spacing wide (1771 nodes inside a 150 um ball sit ~20 um apart) with no
%  regularisation at all, so each fit is forced to pass exactly through every node.
%
%  The only thing that differs is the DATA RANGE, and both use ALL nodes inside it:
%      R = 150   1771 nodes
%      R = 250   8225 nodes
%
%  One figure, on the x_a axis (P1 excited, 1 A) over |x_a| <= 150: b.b itself, in
%  mT^2.  The gradient is still computed and its gate verdict printed (d3 sign
%  changes per actuator axis), but it is no longer plotted -- the field alone shows
%  what the unregularised fit does.
%
% Output: figures/long2016_hexapole_halfcut/xa_lam0_bb.png

clear;  clc;

RS   = [150 250];                        % the two data ranges, all nodes in each
RHO  = 20;   LAM = 0;                    % the shared kernel
SPAN = 150;                              % plotted segment |x_a| <= SPAN
NS   = 601;   NQ = 1201;
POLE = [1 3 6];   SG = [1 1 -1];

here=fileparts(mfilename('fullpath'));  FMX=fileparts(fileparts(here));
MAIN=fileparts(fileparts(fileparts(FMX)));  FLX=fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
FIG = fullfile(FMX,'figures','long2016_hexapole_halfcut','current');

s  = linspace(-SPAN, SPAN, NS).';   Pq = [s zeros(NS,2)];
sG = linspace(-SPAN, SPAN, NQ).';
Ix = zeros(6,1);  Ix(1) = 1;

BB = cell(1,numel(RS));   G = cell(1,numel(RS));   LAB = cell(1,numel(RS));
for q = 1:numel(RS)
    R0 = RS(q);
    [bq, gq, iq] = rbf_field(R0, struct('rho',RHO,'lam',LAM));
    BB{q} = sum(bq(Pq, Ix).^2, 2);              % b.b            [mT^2]
    Gx    = gq(Pq, Ix, 'du');   G{q} = Gx(:,1); % d(b.b)/dx_a    [mT^2/um]
    nch   = zeros(1,3);
    for a = 1:3
        Pv = zeros(NQ,3);  Pv(:,a) = sG;
        I  = zeros(6,1);   I(POLE(a)) = 1;
        d3 = gq(Pv, I, 'd3');
        v  = SG(a)*d3(:,a);   g = sign(v);   nz = find(g~=0);
        nch(a) = sum(diff(g(nz)) ~= 0);
    end
    if all(nch==0), verdict = 'PASS'; else, verdict = 'FAIL'; end
    LAB{q} = sprintf('R = %d \\mum', R0);
    fprintf(['R = %3d : %5d nodes | b.b %8.2f .. %8.2f mT^2 | grad %+8.3f .. %+8.3f | d3 [%d %d %d] %s' newline], ...
            R0, size(iq.P,1), min(BB{q}), max(BB{q}), min(G{q}), max(G{q}), nch, verdict);
    clear bq gq iq
end

draw_(s, BB, LAB, SPAN, '$\mathbf{b\cdot b\;(mT^{2})}$',              false, fullfile(FIG,'xa_lam0_bb.png'));

% ---- local functions ----------------------------------------------------------
function draw_(s, Y, LAB, SPAN, ylab, sym, out)
    FS = 60;  FSLEG = 45;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;  LW = 7;
    COL = {[0.85 0.10 0.10], [0.05 0.10 0.95]};
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]);  hold(ax,'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end

    h = gobjects(1,numel(Y));
    for q = 1:numel(Y), h(q) = plot(ax, s, Y{q}, '-', 'Color', COL{q}, 'LineWidth', LW); end

    lo = min(cellfun(@min, Y));   hi = max(cellfun(@max, Y));
    if sym || lo < 0, [YR, YT] = ylim_sym_(max(abs([lo hi])));
    else,             [YR, YT] = ylim_odd_(hi);  end
    XR = [-SPAN SPAN];   XT = [-SPAN/2 0 SPAN/2];
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
    xlabel(ax, '$\mathbf{x_a\;(\mu m)}$', 'Interpreter','latex','FontSize',FSLAB);
    ylabel(ax, ylab,                      'Interpreter','latex','FontSize',FSLAB);
    lg = legend(ax, h, LAB, 'Interpreter','tex', 'Location','northwest', 'NumColumns',1);
    lg.FontSize = FSLEG;  lg.FontWeight = 'bold';
    lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.Color = 'w';
    lg.ItemTokenSize = [55 25];
    hold(ax,'off');
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
    print(fig, out, '-dpng', '-r200');   close(fig);
    fprintf(['wrote %s' newline], out);
end

function [lim, tk] = ylim_sym_(maxabs)
    s = nice_(maxabs/2);   lim = [-2*s 2*s];   tk = [-s 0 s];
end

function [lim, tk] = ylim_odd_(maxv)
    s = nice_(maxv/4);     lim = [0 4*s];      tk = (1:3)*s;
end

function s = nice_(need)
    cand = [1 1.1 1.2 1.25 1.5 1.6 1.75 2 2.25 2.5 3 4 5 6 7.5 8];
    for kk = floor(log10(need)) : floor(log10(need))+1
        ok = cand(cand*10^kk >= need*(1-1e-12)) * 10^kk;
        if ~isempty(ok), s = min(ok);  return, end
    end
    error('nice_:none','no step for %g', need);
end
