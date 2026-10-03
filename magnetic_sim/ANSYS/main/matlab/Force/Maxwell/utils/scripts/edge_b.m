%% edge_b.m -- B_x, B_y, B_z along ONE tetrahedron edge, straight from the Maxwell export.
%
%  The edge joins two mesh nodes of the 2026-09-20 22:15 solve (current.ngmesh):
%     pid 292093  ( 0.0035010,  0.0574076, -12.7129041) mm
%     pid 277693  ( 0.0147323,  0.0551877, -12.6965238) mm
%  length 19.9846 um, shared by 5 tets, all of them in body Sphere3 (air), neither end node
%  on a surface facet.  The 201 samples run from 2 % to 98 % of the edge so the shared nodes
%  themselves are never evaluated.
%
%  COMPONENTS, NOT |b|.  The basis is linear in the components; |b| = sqrt(B.B) is nonlinear
%  and would show curvature even for a linear B.
%
%  x axis = arc length from node A [um].  One panel per component, own y scale -- the three
%  ranges differ by an order of magnitude and a shared axis would flatten two of them.
%
%  [ADDED 2026-09-21]  temp_code one-off (compute + draw in one file, per the exception in
%  plot-scripts-pure.md).  Output: figures/long2016_hexapole_halfcut/edge_b.png
clear; clc;

HERE = fileparts(mfilename('fullpath'));
FMX  = fileparts(fileparts(HERE));
FIG  = fullfile(FMX,'figures','long2016_hexapole_halfcut','current');
DAT  = fullfile(FMX,'utils','data');
FLD  = 'D:\Maxwell_sim\long2016_hexapole_halfcut\export\No gap\B_p1_0.02_pts.fld';
PTS  = 'D:\Maxwell_sim\long2016_hexapole_halfcut\export\No gap\edge_201.pts';

A = [ 3.50100731847553e-06,  5.74076454151674e-05, -1.27129041335509e-02];   % pid 292093 [m]
Bn = [ 1.47322629061608e-05,  5.51876593987401e-05, -1.26965238087826e-02];  % pid 277693 [m]

M = readmatrix(FLD, 'FileType','text', 'NumHeaderLines',2);
Q = readmatrix(PTS, 'FileType','text', 'NumHeaderLines',1);
assert(size(M,1) == size(Q,1), 'edge_b:n', 'export %d points, list %d', size(M,1), size(Q,1));
assert(max(abs(M(:,1:3) - Q), [], 'all') < 1e-12, 'edge_b:xyz', 'export points differ from the list');

P  = M(:,1:3)*1e-3;                       % mm -> m
B  = M(:,4:6)*1e3;                        % T  -> mT
d  = (Bn - A) / norm(Bn - A);
t  = ((P - A) * d.') * 1e6;               % arc length from node A [um]
off = max(vecnorm((P - A) - (t*1e-6)*d, 2, 2));
fprintf('%d points, %.4f .. %.4f um along a %.4f um edge, off-edge %.2e um\n', ...
        numel(t), t(1), t(end), norm(Bn-A)*1e6, off*1e6);
for c = 1:3
    fprintf('  B%c  %+10.6f .. %+10.6f mT   (span %.2e mT)\n', 'x'+c-1, ...
            min(B(:,c)), max(B(:,c)), max(B(:,c))-min(B(:,c)));
end
save(fullfile(DAT,'edge_b.mat'), 't','B','P','A','Bn','FLD','PTS');

%% ---- draw -----------------------------------------------------------------------
FS = 44;  FSLAB = 36;  LWBOX = 5.0;  LW = 7;  CANV = 14.5;     % 3-panel stack: ticks 44
COL = [0.05 0.10 0.95];
LBL = {'$\mathbf{B_x\;(mT)}$', '$\mathbf{B_y\;(mT)}$', '$\mathbf{B_z\;(mT)}$'};

XR = [0 20];   XT = [5 10 15];
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
H = 0.255;   BOT = [0.700 0.420 0.140];
for c = 1:3
    ax = axes('Parent',fig,'Position',[0.200 BOT(c) 0.740 H]);  hold(ax,'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    plot(ax, t, B(:,c), '-', 'Color', COL, 'LineWidth', LW);
    [YR, YT] = ylim_win_(min(B(:,c)), max(B(:,c)));
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
           'TickLength',[.02 .02],'TickDir','out','Box','on','Layer','top');
    xlim(ax, XR);   set(ax,'XTick', XT);
    ylim(ax, YR);   set(ax,'YTick', YT);
    ylabel(ax, LBL{c}, 'Interpreter','latex', 'FontSize',FSLAB);
    if c < 3
        set(ax,'XTickLabel',[]);                       % rule 8: only the bottom panel labels x
    else
        for xv = XR
            text(ax, xv, YR(1)-0.055*diff(YR), sprintf('%g', xv), ...
                 'HorizontalAlignment','center','VerticalAlignment','top', ...
                 'FontSize',FS,'FontWeight','bold','Clipping','off');
        end
        xlabel(ax, '$\mathbf{arc\;length\;from\;node\;A\;(\mu m)}$', ...
               'Interpreter','latex', 'FontSize',FSLAB);
    end
    hold(ax,'off');
end
out = fullfile(FIG, 'edge_b.png');
set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');   close(fig);
fprintf('wrote %s\n', out);

% a window that does not contain zero: centre on the data, 4 steps across, 3 interior ticks
function [lim, tk] = ylim_win_(lo, hi)
    c = 0.5*(lo+hi);   need = max(hi-lo, eps)/4 * 1.15;
    cand = [1 1.1 1.2 1.25 1.5 1.6 1.75 2 2.25 2.5 3 4 5 6 7.5 8];
    for kk = floor(log10(need)) : floor(log10(need))+1
        ok = cand(cand*10^kk >= need*(1-1e-12)) * 10^kk;
        if isempty(ok), continue, end
        st = min(ok);   c = round(c/st)*st;
        lim = [c-2*st c+2*st];   tk = c + [-st 0 st];   return
    end
    error('ylim_win_:none','no step for %g..%g', lo, hi);
end
