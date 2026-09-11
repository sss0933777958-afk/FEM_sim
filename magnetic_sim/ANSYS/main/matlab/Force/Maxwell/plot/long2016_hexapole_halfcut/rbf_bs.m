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
RHO  = 20;                       % [um]
LAM  = 0;
COIL = 1;                        % P1
NQ   = 1001;  XR = [-150 150];

% [MODIFIED 2026-09-11] moved into matlab/Force/Maxwell/plot/long2016_hexapole_halfcut/.
%   HERE is now .../plot/<model>, so MAIN sits five levels up, and the figure goes to
%   the package's own figures/<model>/ instead of temp_figures.  The .mat inputs are
%   still read from temp_code/data via MAIN.
HERE = fileparts(mfilename('fullpath'));                      % .../plot/<model>
FMX  = fileparts(fileparts(HERE));                            % .../matlab/Force/Maxwell
MAIN = fileparts(fileparts(fileparts(FMX)));                  % .../main
FIG  = fullfile(FMX, 'figures', 'long2016_hexapole_halfcut'); % figure output dir
DAT  = fullfile(FMX, 'data', 'long2016_hexapole_halfcut', '.mat');  % .mat home
%   [MODIFIED 2026-09-11] the rbf_*.mat now live with the Force package, not temp_code.
S = load(fullfile(DAT,'rbf_w6_r20l0.mat'));
sel = vecnorm(S.P,2,2) <= RFIT;
P   = S.P(sel,:);   B = S.B6(sel,:,COIL);   Np = size(P,1);

t = linspace(XR(1), XR(2), NQ).';
Q = [t, zeros(NQ,2)];

s2  = sum(P.^2, 2);
D2n = max(s2 + s2.' - 2*(P*P.'), 0);  D2n(1:Np+1:end) = 0;
Phi = exp(-D2n / RHO^2);
W   = (Phi + LAM*eye(Np)) \ B;
D2q = (Q(:,1)-P(:,1).').^2 + (Q(:,2)-P(:,2).').^2 + (Q(:,3)-P(:,3).').^2;
b   = exp(-D2q / RHO^2) * W;
U   = sum(b.^2, 2);

fprintf('fitting nodes r <= %g um : Np = %d (settled run uses %d at r <= 250)\n', ...
        RFIT, Np, S.Np);
fprintf('rho = %g um, lam = %g | cond(Phi) = %.3e | max|W| = %.3e\n', ...
        RHO, LAM, cond(Phi), max(abs(W),[],'all'));
fprintf('nodal residual max |Phi W - B| = %.3e mT\n\n', max(abs(Phi*W - B),[],'all'));

% [MODIFIED 2026-09-11] the trilinear reference moved out of temp_code (cleared) into
% the Force package's own data folder, next to the rbf_*.mat.
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
save(fullfile(DAT, [STEM '.mat']), 't','U','b','RFIT','RHO','LAM','Np','COIL');

% ================================ figure ======================================
FS = 60;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;  LW = 8;  CBLU = [0.05 0.10 0.95];
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]); hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
% [MODIFIED 2026-09-11] user: draw it as a continuous curve, not 1001 markers.
plot(ax, t, U, '-', 'Color',CBLU, 'LineWidth',LW);

set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
       'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
       'XTickLabelRotation',0,'YTickLabelRotation',0);
sx = diff(XR)/4;
xlim(ax, XR);  set(ax,'XTick', XR(1)+(1:3)*sx);
[yr, yt] = ylim_nice_(max(U), 3);
ylim(ax, yr);  set(ax,'YTick', yt);  ytop = yr(2);
for xv = XR
    text(ax, xv, -0.022*ytop, sprintf('%g', xv), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{x_a\;(\mu m)}$',       'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{b\cdot b\;(mT^2)}$',   'Interpreter','latex', 'FontSize',FSLAB);
hold(ax,'off');

out = fullfile(FIG,[STEM '.png']);
set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');  close(fig);
fprintf(['wrote %s  (y ticks %s, top %g)' newline], out, mat2str(yt), ytop);

function [lim, tk] = ylim_nice_(maxv, N)
    if nargin < 2 || isempty(N), N = 3; end
    smin = 1.08*maxv/(N+1);
    u    = 10^floor(log10(smin));
    s    = ceil(smin/u - 1e-12) * u;
    lim  = [0 (N+1)*s];   tk = (1:N)*s;
end
