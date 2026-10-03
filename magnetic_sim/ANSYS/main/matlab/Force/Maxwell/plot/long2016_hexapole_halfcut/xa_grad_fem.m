% xa_grad_fem  b.b and d(b.b)/dx_a on the x_a axis, straight from the Maxwell
%              solve for P1 -- no fitted model anywhere.
%
%   MESH picks the solve: '0.1' -> tag h0p1, '0.02' -> tag h0p02.  Both were
%   exported on the same 0.02 mm Cartesian grid, so the curves overlay directly.
%
%   QTY = 'bb'    b.b                [mT^2]     -> xa_bb_<tag>.png
%   QTY = 'grad'  d(b.b)/dx_a        [mT^2/um]  -> xa_grad_<tag>.png
%   QTY = 'gradp' the same, but evaluated BY AEDT at the points
%                                    [mT^2/um]  -> xa_gradp_<tag>.png
%   QTY = 'gradsm' the same again, but from AEDT's SMOOTHED stack
%                 Smooth(Grad(Smooth(Dot(B,B))))
%                                    [mT^2/um]  -> xa_gradsm_<tag>.png
%   QTY = 'bbp'   b.b evaluated BY AEDT at the points
%                                    [mT^2]     -> xa_bbp_<tag>.png
%   QTY = 'res'   p(x_a) - b.b       [mT^2]     -> xa_res_<tag>.png
%   QTY = 'resp'  the same residual, for the 'bbp' curve
%                                    [mT^2]     -> xa_resp_<tag>.png
%   QTY = 'cmp'   'grad' and 'gradp' overlaid   -> xa_cmp_<tag>.png
%
%   'bb' and 'res' need the B export of that solve (only the 0.1 mm run has one);
%   'gradp' and 'cmp' need the point-list export (only the 0.02 mm run has one).
%
%   TWO ROUTES TO THE SAME QUANTITY.  'grad' goes AEDT -> Cartesian grid -> our
%   trilinear interpolation; 'gradp' hands AEDT the same 1001 points and lets it
%   interpolate in its own tetrahedral mesh.  Neither is uninterpolated -- the
%   comparison is one interpolation against the other.
%
%   Over |s| <= 150 um, 1001 equidistant points on the x_a axis of the actuator
%   frame.  A single series gets no legend (the caption carries it).
%
%   'res' is the degree-DEG least-squares polynomial minus the curve it was fitted
%   to, so it is the LS RESIDUAL: it averages to zero by construction and changes
%   sign.  Its y axis is therefore symmetric about zero and a dashed zero line is
%   drawn; KINKS may mark the trilinear cell-face crossings on it.
%
%   The gradient is the GRADIENT, not the force; the force is 22.55 times it
%   (0.5*mgB*UF).
%
%   Reads only (plot-scripts-pure); xa_fem_<tag>.mat is written by
%   temp_code/xa_grad_fem.m.
%
% Output: figures/long2016_hexapole_halfcut/xa_{bb,grad,gradp,res,cmp}_<tag>.png
%
% [ADDED 2026-09-20]

MESH  = '0.02';             % '0.1' | '0.02' -- which solve's .mat to draw
QTY   = 'bpts';               % 'bb' | 'grad' | 'gradp' | 'gradsm' | 'bpts' | 'bbp' | 'res' | 'resp' | 'cmp'
ZOOM  = [-72 -50];            % [] = the whole +-SPAN; [s0 s1] draws that window only
                              %   (file gets a _zoom suffix).  Use it to see whether the
                              %   curve is straight between element faces.
KINKS = 1;                  % 'res' only: which cell-face families to mark, [] = none
                            %   1 = x faces (every 24.5 um, the dominant kinks)
                            %   3 = z faces (every 34.6 um, weaker)  |  [1 3] = both

HERE = fileparts(mfilename('fullpath'));
FMX  = fileparts(fileparts(HERE));
DAT  = fullfile(FMX,'utils','data');
FIG  = fullfile(FMX,'figures','long2016_hexapole_halfcut','current');

QTY = validatestring(QTY, {'bb','grad','gradp','gradsm','bpts','bbp','res','resp','cmp'}, mfilename, 'QTY');
switch MESH
    case '0.1',  TAG = 'h0p1';
    case '0.02', TAG = 'h0p02';
    otherwise,   error('xa_grad_fem:mesh','MESH must be ''0.1'' or ''0.02''');
end
S = load(fullfile(DAT, sprintf('xa_fem_%s.mat', TAG)));
fprintf(['%s' newline], S.SOLVE);

GLAB = '$\mathbf{d(b\cdot b)/dx_a\;(mT^{2}/\mu m)}$';
SYM  = false;
switch QTY
    case 'bb'
        Y = {S.bb};      NAME = {''};   ylab = '$\mathbf{b\cdot b\;(mT^{2})}$';   unit = S.UNIT_BB;
    case 'grad'
        Y = {S.g_xa};    NAME = {''};   ylab = GLAB;   unit = S.UNIT_G;
    case 'gradp'
        Y = {S.g_pts};   NAME = {''};   ylab = GLAB;   unit = S.UNIT_G;
    case 'bpts'
        % [ADDED 2026-09-21] |b| from AEDT's B_Vector at the points
        Y = {S.bmag_pts};  NAME = {''};   ylab = '$\mathbf{|b|\;(mT)}$';   unit = 'mT';
    case 'gradsm'
        % [ADDED 2026-09-21] SIGNED here: this export comes out negative over the whole axis
        % (the unsmoothed one is positive), so it gets the symmetric axis + zero line.
        Y = {S.g_pts_sm};  NAME = {''};   ylab = GLAB;   unit = S.UNIT_G;   SYM = true;
    case 'bbp'
        Y = {S.bb_pts};  NAME = {''};   ylab = '$\mathbf{b\cdot b\;(mT^{2})}$';   unit = S.UNIT_BB;
    case 'resp'
        Y = {S.res_pts}; NAME = {''};   ylab = '$\mathbf{p(x_a)-b\cdot b\;(mT^{2})}$';
        unit = S.UNIT_BB;   SYM = true;
    case 'res'
        Y = {S.res};     NAME = {''};   ylab = '$\mathbf{p(x_a)-b\cdot b\;(mT^{2})}$';
        unit = S.UNIT_BB;   SYM = true;
    case 'cmp'
        Y = {S.g_xa, S.g_pts};
        NAME = {'grid + trilinear', 'AEDT at points'};
        ylab = GLAB;   unit = S.UNIT_G;
end
for q = 1:numel(Y)
    assert(~isempty(Y{q}), 'xa_grad_fem:noqty', ...
           ['''%s'' is not available for mesh %s mm -- ''bb''/''res'' need that solve''s B ' ...
            'export, ''gradp''/''cmp'' need its point-list export.'], QTY, S.MESH);
end
s   = S.s;
if strcmp(QTY,'bpts') && isfield(S,'s_b') && ~isempty(S.s_b)
    s = S.s_b(:);          % that export has its own (denser) point list
end
out = fullfile(FIG, sprintf('xa_%s_%s.png', QTY, TAG));
for q = 1:numel(Y)
    fprintf(['  %-18s %+.4f .. %+.4f %s' newline], NAME{q}, min(Y{q}), max(Y{q}), unit);
end
if strcmp(QTY,'res')
    fprintf(['degree-%d fit, max |rel| %.4f %%, rms |rel| %.4f %%' newline], ...
            S.DEG, 100*max(abs(S.res./S.bb)), 100*rms(S.res./S.bb));
elseif strcmp(QTY,'resp')
    fprintf(['degree-%d fit, max |rel| %.4f %%, rms |rel| %.4f %%' newline], ...
            S.DEG, S.DCHK.max_rel_pct, S.DCHK.rms_rel_pct);
end

FS = 60;  FSLEG = 45;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;  LW = 7;
COL = {[0.85 0.10 0.10], [0.05 0.10 0.95]};   % 'cmp' overlay: red then blue
COL1 = [0.05 0.10 0.95];                      % [MODIFIED 2026-09-21 user] a single series
                                              % is drawn blue; the overlay keeps red+blue
if numel(Y) == 1, COL{1} = COL1; end
STY = {'-', '-'};

if isempty(ZOOM)
    [XR, XT] = axis_sym_(S.SPAN);
    vis = true(size(s));
else
    % [ADDED 2026-09-21] zoom: x from the window, y from what is visible in it, so the
    % element-scale structure is not flattened by the full-range limits.
    XR = [min(ZOOM) max(ZOOM)];   XT = XR(1) + (1:3)*diff(XR)/4;
    vis = s >= XR(1) & s <= XR(2);
    out = strrep(out, '.png', '_zoom.png');
end
ally = cell2mat(cellfun(@(v) v(vis), Y, 'UniformOutput', false));
if ~isempty(ZOOM)
    [YR, YT] = ylim_win_(min(ally(:)), max(ally(:)));
elseif SYM
    [YR, YT] = ylim_sym_(max(abs(ally(:))));
else
    [YR, YT] = ylim_odd_(max(ally(:)));
end

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]);  hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end

hL = gobjects(0);   lbl = {};
if strcmp(QTY,'res') && ~isempty(KINKS)   % only the grid route has Cartesian kinks;
                                          % the point-wise one is cut by tet faces instead
    % Where the trilinear interpolant kinks: one vertical line per cell-face
    % crossing, drawn first so the residual stays on top.  Each family goes in as a
    % single NaN-separated line object, so it takes exactly one legend entry.
    cF = {[0.05 0.10 0.95], [], [0.00 0.55 0.20]};   sk = {'--', '', ':'};
    for d = KINKS(:).'
        xk = S.KINK.s{d}(:).';
        if isempty(xk), continue, end
        XX = [xk; xk; nan(1,numel(xk))];
        YY = repmat([YR(1); YR(2); NaN], 1, numel(xk));
        hL(end+1) = plot(ax, XX(:), YY(:), sk{d}, 'Color', cF{d}, 'LineWidth', LW-4); %#ok<SAGROW>
        lbl{end+1} = sprintf('%s faces (%.1f \\mum)', S.KINK.name{d}, S.KINK.period(d)); %#ok<SAGROW>
    end
end
if SYM                                            % zero reference for a signed quantity
    plot(ax, XR, [0 0], '--', 'Color', [0.45 0.45 0.45], 'LineWidth', LW-3);
end
for q = 1:numel(Y)
    h = plot(ax, s, Y{q}, STY{q}, 'Color', COL{q}, 'LineWidth', LW);
    if numel(Y) > 1, hL(end+1) = h;  lbl{end+1} = NAME{q};  end %#ok<SAGROW>
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
xlabel(ax, '$\mathbf{x_a\;(\mu m)}$', 'Interpreter','latex','FontSize',FSLAB);
ylabel(ax, ylab,                      'Interpreter','latex','FontSize',FSLAB);

if ~isempty(hL)                                   % tex, not latex: latex ignores bold
    loc = 'south';   if ~SYM, loc = 'northwest';  end
    lg = legend(ax, hL, lbl, 'Interpreter','tex', 'Location',loc, 'NumColumns',1);
    lg.FontSize = FSLEG;   lg.FontWeight = 'bold';
    lg.Box = 'on';   lg.EdgeColor = 'k';   lg.LineWidth = LWBOX;   lg.Color = 'w';
    lg.ItemTokenSize = [55 25];
end

% A long y label does not fit the fixed left margin and is silently clipped at the
% canvas edge ('res' needs 0.198 against the 0.190 on offer).  Give it the room it
% asks for, keeping the right edge put; the short labels ask for less than the
% margin already is, so those figures come out byte-for-byte as before.
drawnow;
need = ax.TightInset(1) + 0.010;
if ax.Position(1) < need
    rt = ax.Position(1) + ax.Position(3);
    ax.Position(1) = need;   ax.Position(3) = rt - need;
    fprintf(['  [layout] y label needed more room: left margin %.4f -> %.4f' newline], ...
            0.1895, need);
end
hold(ax,'off');

set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');   close(fig);
fprintf(['wrote %s' newline], out);

% ---- local helpers ------------------------------------------------------------
% A window that does not contain zero: centre on the data, 4 steps across, 3 interior
% ticks -- same tick/label rules as the other two helpers.
function [lim, tk] = ylim_win_(lo, hi)
    c = 0.5*(lo+hi);   need = max(hi-lo, eps)/4 * 1.15;
    cand = [1 1.1 1.2 1.25 1.5 1.6 1.75 2 2.25 2.5 3 4 5 6 7.5 8];
    for kk = floor(log10(need)) : floor(log10(need))+1
        ok = cand(cand*10^kk >= need*(1-1e-12)) * 10^kk;
        if isempty(ok), continue, end
        st = min(ok);
        c  = round(c/st)*st;
        lim = [c-2*st c+2*st];   tk = c + [-st 0 st];   return
    end
    error('ylim_win_:none','no step for %g..%g', lo, hi);
end
function [lim, tk] = axis_sym_(R)
    s = R/2;   lim = [-R R];   tk = [-s 0 s];
end

% Signed quantity: symmetric about zero.  lim = [-2s, 2s], ticks = [-s 0 s], so the
% gap to each end equals the tick spacing (rule 5) and neither end is labelled
% (rule 4).  Step chosen as in ylim_odd_: smallest two-significant-figure step that
% fits, unless a round one wastes no more than 20 % more.
function [lim, tk] = ylim_sym_(maxabs)
    cand = [1 1.1 1.2 1.25 1.5 1.6 1.75 2 2.25 2.5 3 4 5 6 7.5 8];
    nice = [1 1.5 2 2.5 3 4 5];   need = maxabs/2;
    for kk = floor(log10(need)) : floor(log10(need))+1
        ok = cand(cand*10^kk >= need*(1-1e-12)) * 10^kk;
        if isempty(ok), continue, end
        s = min(ok);
        okn = nice(nice*10^kk >= need*(1-1e-12)) * 10^kk;
        if ~isempty(okn) && min(okn) <= 1.2*s, s = min(okn); end
        lim = [-2*s 2*s];   tk = [-s 0 s];   return
    end
    error('ylim_sym_:none','no step for %g', maxabs);
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
