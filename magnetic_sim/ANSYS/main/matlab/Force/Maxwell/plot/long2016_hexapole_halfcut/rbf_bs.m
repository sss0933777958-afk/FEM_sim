% rbf_bs  U = b.b along the xa axis when the FITTING node set is restricted.
%
%   RFIT = 150 um  ->  only the 1771 nodes inside the working region are used as RBF
%   centres, instead of the settled r <= 250 um (8225).  rho = 20 um, lam = 0 (exact
%   interpolation).  Evaluated along the whole xa axis, +-150 um, 1001 points, P1.
%
%   The evaluation line therefore REACHES the edge of the node cloud, which is the point
%   of the figure: a Gaussian RBF has no centres beyond the boundary to hold the surface
%   up, so it collapses towards zero there.
%
% Output: figures/long2016_hexapole_halfcut/rbf_bs.png + temp_code/data/rbf_bs.mat

RFIT = 150;                      % [um] radius of the FITTING node set
% [ADDED 2026-09-11] override + tagged output, so the r<=150 and r<=250 node sets
% each keep their own figure.
if exist('RFIT_OVR','var'), RFIT = RFIT_OVR; end
STEM = sprintf('rbf_bs%d', RFIT);
% [ADDED 2026-09-27] rho sweep. One entry keeps the original single-curve figure
% untouched; several entries overlay them and write a separate, tagged file, so
% the question "does a wider kernel cure the boundary collapse at lam = 0?" can
% be read straight off one picture.
RHO_LIST = 20;                   % [um]  e.g. [20 40 60 80 120]
if exist('RHO_OVR','var'), RHO_LIST = RHO_OVR; end
SWEEP = numel(RHO_LIST) > 1;
if SWEEP, STEM = sprintf('rbf_bs%d_rho', RFIT); end
RHO  = RHO_LIST(1);              % [um]
LAM  = 0;
% [ADDED 2026-09-27] lam override + a stem tag whenever rho/lam leave the
% defaults, so a sweep entry can never overwrite the rho = 20, lam = 0 figure.
if exist('LAM_OVR','var'), LAM = LAM_OVR; end
if ~SWEEP && (RHO ~= 20 || LAM ~= 0)
    STEM = sprintf('%s_r%g', STEM, RHO);
    if LAM ~= 0, STEM = [STEM 'l' strrep(sprintf('%g',LAM),'-','m')]; end
end
% [ADDED 2026-09-27] which quantity to draw along the axis.
%   'bb' : U = b.b            'du' : dU/dx_a, the quantity the force model uses
QTY = 'bb';
if exist('QTY_OVR','var'), QTY = QTY_OVR; end
if strcmpi(QTY,'du'), STEM = [STEM '_du']; end
COIL = 1;                        % P1
NQ   = 1001;  XR = [-150 150];

% [MODIFIED 2026-09-11] moved into matlab/Force/Maxwell/plot/long2016_hexapole_halfcut/.
%   HERE is now .../plot/<model>, so MAIN sits five levels up, and the figure goes to
%   the package's own figures/<model>/ instead of temp_figures.  The .mat inputs are
%   still read from temp_code/data via MAIN.
HERE = fileparts(mfilename('fullpath'));                      % .../plot/<model>
FMX  = fileparts(fileparts(HERE));                            % .../matlab/Force/Maxwell
MAIN = fileparts(fileparts(fileparts(FMX)));                  % .../main
FIG  = fullfile(FMX,'figures','long2016_hexapole_halfcut','current'); % figure output dir
DAT  = fullfile(FMX,'utils','data');  % .mat home
%   [MODIFIED 2026-09-11] the rbf_*.mat now live with the Force package, not temp_code.
S = load(fullfile(DAT,'rbf_w6_r20l0.mat'));
sel = vecnorm(S.P,2,2) <= RFIT;
P   = S.P(sel,:);   B = S.B6(sel,:,COIL);   Np = size(P,1);

t = linspace(XR(1), XR(2), NQ).';
Q = [t, zeros(NQ,2)];

s2  = sum(P.^2, 2);
D2n = max(s2 + s2.' - 2*(P*P.'), 0);  D2n(1:Np+1:end) = 0;
D2q = (Q(:,1)-P(:,1).').^2 + (Q(:,2)-P(:,2).').^2 + (Q(:,3)-P(:,3).').^2;

fprintf('fitting nodes r <= %g um : Np = %d (settled run uses %d at r <= 250)\n', ...
        RFIT, Np, S.Np);
Uall = nan(NQ, numel(RHO_LIST));
fprintf('%7s %11s %11s %13s %11s\n','rho','rcond(Phi)','max|W|','nodal resid','U(+150)');
for k = 1:numel(RHO_LIST)
    rk  = RHO_LIST(k);
    Phi = exp(-D2n / rk^2);
    ws  = warning('off','MATLAB:nearlySingularMatrix');
    W   = (Phi + LAM*eye(Np)) \ B;
    warning(ws);
    K   = exp(-D2q / rk^2);
    b   = K * W;
    if strcmpi(QTY,'du')
        % d/dx_a of exp(-r^2/rho^2) = -(2/rho^2)(x_a - p_a) exp(...), so
        % dU/dx_a = 2 sum_c b_c db_c/dx_a.
        dbx = (-(2/rk^2) * ((Q(:,1) - P(:,1).') .* K)) * W;
        Uall(:,k) = 2 * sum(b .* dbx, 2);
    else
        Uall(:,k) = sum(b.^2, 2);
    end
    fprintf('%7g %11.3e %11.3e %13.3e %11.4f\n', rk, rcond(Phi), ...
            max(abs(W),[],'all'), max(abs(Phi*W - B),[],'all'), Uall(end,k));
end
U = Uall(:,1);                                  % the single-curve path is k = 1
fprintf('\n');

% [MODIFIED 2026-09-11] the trilinear reference moved out of temp_code (cleared) into
% the Force package's own data folder, next to the rbf_*.mat.
% The stored reference is b.b, so it only means anything for QTY = 'bb'.
if ~SWEEP && strcmpi(QTY,'bb')
T = load(fullfile(DAT,'bb_xa_raw.mat'));     % trilinear reference
fprintf('%8s %12s %12s %10s\n','xa[um]','b.b RBF','b.b trilin','diff');
for tt = [-150 -120 -75 0 75 120 150]
    i = find(abs(t-tt) < 1e-9, 1);   j = find(abs(T.t-tt) < 1e-9, 1);
    fprintf('%8.0f %12.4f %12.4f %10.4f\n', t(i), U(i), T.s(j), U(i)-T.s(j));
end
if isequal(T.t(:), t)
    e = U - T.s(:);
    fprintf('\nvs trilinear: mean %+.4f | mean|.| %.4f | max %.4f mT^2 | NMAE %.4f %%\n', ...
            mean(e), mean(abs(e)), max(abs(e)), sum(abs(e))/sum(abs(T.s))*100);
    in = abs(t) <= 120;
    fprintf('             |xa| <= 120 only : mean|.| %.4f | max %.4f mT^2 | NMAE %.4f %%\n', ...
            mean(abs(e(in))), max(abs(e(in))), sum(abs(e(in)))/sum(abs(T.s(in)))*100);
end
end
save(fullfile(DAT, [STEM '.mat']), 't','U','Uall','RHO_LIST','RFIT','RHO','LAM','Np','COIL');

% ================================ figure ======================================
FS = 60;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;  LW = 8;  CBLU = [0.05 0.10 0.95];
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]); hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
% [MODIFIED 2026-09-11] user: draw it as a continuous curve, not 1001 markers.
if SWEEP
    CM = lines(numel(RHO_LIST));   hh = gobjects(1,numel(RHO_LIST));   lb = cell(1,numel(RHO_LIST));
    for k = 1:numel(RHO_LIST)
        hh(k) = plot(ax, t, Uall(:,k), '-', 'Color',CM(k,:), 'LineWidth',LW-2);
        lb{k} = sprintf('\\rho = %g', RHO_LIST(k));
    end
    lg = legend(ax, hh, lb, 'Interpreter','tex', 'FontSize',45, 'FontWeight','bold', ...
                'Location','northwest');
    lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.ItemTokenSize = [55 25];
else
    plot(ax, t, U, '-', 'Color',CBLU, 'LineWidth',LW);
end

set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
       'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
       'XTickLabelRotation',0,'YTickLabelRotation',0);
sx = diff(XR)/4;
xlim(ax, XR);  set(ax,'XTick', XR(1)+(1:3)*sx);
if min(Uall(:)) < 0                     % the gradient changes sign; ylim_nice_ assumes [0, top]
    [yr, yt] = ylim_sym_(min(Uall(:)), max(Uall(:)), 3);
else
    [yr, yt] = ylim_nice_(max(Uall(:)), 3);
end
ylim(ax, yr);  set(ax,'YTick', yt);  ytop = yr(2);
for xv = XR
    text(ax, xv, -0.022*ytop, sprintf('%g', xv), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{x_a\;(\mu m)}$',       'Interpreter','latex', 'FontSize',FSLAB);
if strcmpi(QTY,'du')
    ylabel(ax, '$\mathbf{d(b\cdot b)/dx_a\;(mT^2/\mu m)}$', 'Interpreter','latex', 'FontSize',FSLAB);
else
    ylabel(ax, '$\mathbf{b\cdot b\;(mT^2)}$',   'Interpreter','latex', 'FontSize',FSLAB);
end
hold(ax,'off');

out = fullfile(FIG,[STEM '.png']);
set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');  close(fig);
fprintf(['wrote %s  (y ticks %s, top %g)' newline], out, mat2str(yt), ytop);

% odd ticks on a range that straddles zero, end gaps equal to the spacing
function [lim, tk] = ylim_sym_(lo, hi, N)
    cands = reshape(([1 1.5 2 2.5 3 4 5 6 8].') .* 10.^(-4:6), [], 1);
    cands = sort(cands);
    for s = cands.'
        a = floor(lo/s)*s;
        if a + (N+1)*s >= hi, tk = a + (1:N)*s;  lim = [a, a+(N+1)*s];  return; end
    end
    error('ylim_sym_:range','no nice step covers [%g %g]', lo, hi);
end

function [lim, tk] = ylim_nice_(maxv, N)
    if nargin < 2 || isempty(N), N = 3; end
    smin = 1.08*maxv/(N+1);
    u    = 10^floor(log10(smin));
    s    = ceil(smin/u - 1e-12) * u;
    lim  = [0 (N+1)*s];   tk = (1:N)*s;
end
