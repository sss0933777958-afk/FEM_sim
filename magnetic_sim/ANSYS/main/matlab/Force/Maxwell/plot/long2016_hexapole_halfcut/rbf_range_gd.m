% rbf_range_gd  F_xa on the xa axis for a ladder of DATA RANGES, one figure.
%
%   Companion to rbf_range_bb.m, reading the same .mat.  The kernel is held fixed
%   at the settled pair (rho = 120 um, lam = 1e-4) and only the node radius R
%   moves.  b.b turned out to be almost blind to R -- the five curves there agree
%   to 0.29 % -- so this is the figure that actually answers whether the pair
%   still holds up: the gradient is where a poorly conditioned weight vector
%   shows itself.
%
%       F_xa = 0.5 * mgB * UF * d(b.b)/dx_a        [pN]
%
%   Smoothness is measured the way it has been throughout: fit a 5th-order
%   polynomial to F_xa and take the peak-to-peak of what is left over.  The
%   settled R = 250 um reference is 0.0410 pN.
%
%   Reads only; see rbf_range_bb.m for the call that produces the .mat.
%
%   [ADDED 2026-09-13] PICK selects which levels to draw, so one level can be
%   inspected on its own without a second script.  A single series gets no
%   legend (figure-style: name it in the caption instead).
%
% Output: figures/long2016_hexapole_halfcut/rbf_range_gd.png
%         ... or rbf_range_gd<R>.png when PICK names a single level

HERE = fileparts(mfilename('fullpath'));                      % .../plot/<model>
FMX  = fileparts(fileparts(HERE));                            % .../Force/Maxwell
DAT  = fullfile(FMX,'utils','data');
FIG  = fullfile(FMX,'figures','long2016_hexapole_halfcut','current');

PICK = [];                        % [] = every level; e.g. 350 or [300 350]
QTY  = 'F';                       % 'F' | 'U' | 'd3' | 'err' (F 減掉第一個保留欄)
XZOOM = [];                       % [] = the whole axis, or [lo hi] in um
SRC  = 'range';                   % 'range' = the R ladder | 'axes' = three axes at one R
if exist('PICK_OVR','var'),  PICK  = PICK_OVR;  end
if exist('QTY_OVR','var'),   QTY   = QTY_OVR;   end
if exist('XZOOM_OVR','var'), XZOOM = XZOOM_OVR; end
if exist('SRC_OVR','var'),   SRC   = SRC_OVR;   end

if strcmpi(SRC,'axes'), MATF = 'rbf_axes.mat'; else, MATF = 'rbf_range.mat'; end
if exist('MATF_OVR','var') && ~isempty(MATF_OVR), MATF = MATF_OVR; end
S = load(fullfile(DAT, MATF));
xa = S.xa;   Fx = S.Fx;   RLIST = S.RLIST;   ok = S.okv(:).';
keep = find(ok);
if ~isempty(PICK), keep = keep(ismember(RLIST(keep), PICK)); end
AXN = {'x_a','y_a','z_a'};        % SRC='axes' 時三欄依序對應這三根軸
PRE = 'rbf_range';   if strcmpi(SRC,'axes'), PRE = 'rbf_axes'; end
STEM = [PRE '_gd'];
if strcmpi(QTY,'d3'), STEM = [PRE '_d3']; end
if strcmpi(QTY,'U'),  STEM = [PRE '_u'];  end
if strcmpi(QTY,'err'),  STEM = [PRE '_err'];  end
if strcmpi(QTY,'aerr'), STEM = [PRE '_aerr']; end
if numel(keep) == 1
    if strcmpi(SRC,'axes')
        STEM = sprintf('%s_%s', STEM, strrep(AXN{keep},'_',''));   % -> _xa / _ya / _za
    else
        STEM = sprintf('%s%d', STEM, RLIST(keep));
    end
end
if ~isempty(XZOOM),   STEM = [STEM 'z']; end
if strcmpi(QTY,'d3')
    assert(isfield(S,'U3'), 'this .mat has no U3; rerun the d3 pass first');
    YY = S.U3;   YLAB = '$\mathbf{d^3(b\cdot b)/dx_a^3}$';
elseif strcmpi(QTY,'U')
    YY = S.Ubb;  YLAB = '$\mathbf{b\cdot b\;(mT^2)}$';
elseif strcmpi(QTY,'err')          % 相對於第一個保留欄的差值
    YY = Fx - Fx(:, find(ok,1));
    YLAB = '$\mathbf{\Delta F_{x_a}\;(pN)}$';
elseif strcmpi(QTY,'aerr')         % 同上，取絕對值
    YY = abs(Fx - Fx(:, find(ok,1)));
    YLAB = '$\mathbf{|\Delta F_{x_a}|\;(pN)}$';
else
    YY = Fx;     YLAB = '$\mathbf{F_{x_a}\;(pN)}$';
end
XLAB = '$\mathbf{x_a\;(\mu m)}$';
if strcmpi(SRC,'axes')
    if numel(keep) == 1                % 單軸：標籤用該軸自己的名字
        AN   = AXN{keep};
        XLAB = ['$\mathbf{' AN '\;(\mu m)}$'];
        if strcmpi(QTY,'F'),  YLAB = ['$\mathbf{F_{' AN '}\;(pN)}$'];  end
        if strcmpi(QTY,'d3'), YLAB = ['$\mathbf{d^3(b\cdot b)/d' AN '^3}$']; end
    else                               % 疊圖：共用一條 s 軸
        XLAB = '$\mathbf{s\;(\mu m)}$';
        if strcmpi(QTY,'F'),  YLAB = '$\mathbf{F_s\;(pN)}$';  end
        if strcmpi(QTY,'d3'), YLAB = '$\mathbf{d^3(b\cdot b)/ds^3}$'; end
    end
end
sel = true(size(xa));
if ~isempty(XZOOM), sel = xa >= XZOOM(1) & xa <= XZOOM(2); end
assert(~isempty(keep), 'no level in rbf_range.mat completed');

fprintf('rho = %g um, lam = %.0e held fixed (R = 250 um reference: wobble 0.0410 pN)\n\n', ...
        S.kern.rho, S.kern.lam);
fprintf('%8s %9s %11s %12s %12s %12s\n', ...
        'R [um]','Np','F(0) [pN]','wobble [pN]','max|W|','resid [mT]');
for k = keep
    fprintf('%8d %9d %11.4f %12.4f %12.3e %12.3e\n', ...
            RLIST(k), S.Np(k), Fx(xa==0,k), S.ppk(k), S.wmax(k), S.resid(k));
end
fprintf('\nspread across the ladder: max over xa of (max_k - min_k) = %.4f pN\n', ...
        max(max(Fx(:,keep),[],2) - min(Fx(:,keep),[],2)));
if isfield(S,'f3neg')                       % the analytic third-derivative test
    fprintf(['%8s %14s %12s %16s' newline], 'R [um]','d3<0 [%]','sign changes','min d3');
    for k = keep
        fprintf(['%8d %14.2f %12d %16.4e' newline], RLIST(k), S.f3neg(k), S.n3sgn(k), S.u3min(k));
    end
end

% ================================ figure ======================================
FS = 60;  FSLEG = 45;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;  LW = 6;
XR = [-150 150];   if ~isempty(XZOOM), XR = XZOOM; end
% 專案標準色：單條藍、兩條藍紅、三條以上才用 turbo
if     numel(keep) == 1, CM = [0.05 0.10 0.95];
elseif numel(keep) == 2, CM = [0.05 0.10 0.95; 0.85 0.10 0.10];
else,                    CM = turbo(numel(keep));  end
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]); hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
hh = gobjects(1,numel(keep));  lbl = cell(1,numel(keep));
for i = 1:numel(keep)
    k = keep(i);
    hh(i)  = plot(ax, xa(sel), YY(sel,k), '-', 'Color',CM(i,:), 'LineWidth',LW);
    if strcmpi(QTY,'d3')            % mark the stretch that fails the sign test
        yk = YY(:,k);   sg = sign(yk);  sg(sg==0) = 1;
        bad = (sg ~= sign(sum(sg))) & sel;      % minority sign = the sign test's failure
        if any(bad)
            yb = nan(size(yk));  yb(bad) = yk(bad);
            plot(ax, xa, yb, '-', 'Color',[0.85 0.10 0.10], 'LineWidth',LW*1.6);
        end
    end
    if strcmpi(SRC,'axes'), lbl{i} = S.LBL{k};
    else,                   lbl{i} = ['R = ' num2str(RLIST(k)) ' \mum'];  end
end

set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
       'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
       'XTickLabelRotation',0,'YTickLabelRotation',0);
sx = diff(XR)/4;
xlim(ax, XR);  set(ax,'XTick', XR(1)+(1:3)*sx);
if strcmpi(QTY,'d3')
    yline(ax, 0, 'k--', 'LineWidth', LWBOX*0.7);      % the sign test lives on this line
    [yr, yt] = ylim_sym_(min(min(YY(sel,keep),[],'all'),0), ...
                         max(max(YY(sel,keep),[],'all'),0), 3);
elseif ~isempty(XZOOM)                                % a zoom fits the data, not zero
    [yr, yt] = ylim_sym_(min(YY(sel,keep),[],'all'), max(YY(sel,keep),[],'all'), 3);
elseif min(YY(sel,keep), [], 'all') < 0
    [yr, yt] = ylim_sym_(min(YY(sel,keep),[],'all'), max(YY(sel,keep),[],'all'), 3);
else
    [yr, yt] = ylim_nice_(max(YY(sel,keep), [], 'all'), 3);
end
ylim(ax, yr);  set(ax,'YTick', yt);  ytop = yr(2);
for xv = XR
    text(ax, xv, yr(1)-0.022*diff(yr), sprintf('%g', xv), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, XLAB, 'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, YLAB, 'Interpreter','latex', 'FontSize',FSLAB);
if numel(keep) > 1                  % a single series is named in the caption
    lg = legend(ax, hh, lbl, 'Interpreter','tex', 'Location','northwest', 'NumColumns',1);
    lg.FontSize = FSLEG*0.8;  lg.FontWeight = 'bold';
    lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.Color = 'w';
    lg.ItemTokenSize = [50 22];
end
hold(ax,'off');

out = fullfile(FIG, [STEM '.png']);
set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');  close(fig);
fprintf(['wrote %s  (y ticks %s, top %g, fill %.0f%%)' newline], ...
        out, mat2str(yt), ytop, 100*max(YY(sel,keep),[],'all')/ytop);

function [lim, tk] = ylim_sym_(lo, hi, N)
% Signed data: N equally spaced inner ticks, the gaps to both frame edges equal
% to the tick spacing, and zero guaranteed inside the view.
    if nargin < 3 || isempty(N), N = 3; end
    s  = 1.02*(hi-lo)/(N+1);
    u  = 10^floor(log10(s));   s = ceil(s/u - 1e-12)*u;
    for it = 1:8
        c   = round(((lo+hi)/2)/s)*s;
        lim = [c-((N+1)/2)*s, c+((N+1)/2)*s];
        if lim(1) <= lo && lim(2) >= hi, break; end
        s = s + u;
    end
    tk = c + (-(N-1)/2:(N-1)/2)*s;
end

function [lim, tk] = ylim_nice_(maxv, N)
    if nargin < 2 || isempty(N), N = 3; end
    smin = 1.08*maxv/(N+1);
    u    = 10^floor(log10(smin));
    s    = ceil(smin/u - 1e-12) * u;
    lim  = [0 (N+1)*s];   tk = (1:N)*s;
end
