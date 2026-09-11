% rbf_rsw  rho sweep at fixed lambda: NMAE of the predicted field vs rho.
%
%   lambda = 1e-4 fixed, rho = 20 : 10 : 100 um (node spacing h = 20 um).
%   Nodes r <= 250 um (Np = 8225), six excitations solved together (18 RHS).
%   Scored at the 1771 nodes with r <= RIN, all six coils:
%       NMAE = sum_i || b_pred,i - B_i ||  /  sum_i || B_i ||  * 100
%
%   Companion to rbf_lsw.m (lambda sweep at fixed rho).  Each rho needs its own Phi,
%   so this is nine builds and nine factorisations -- no eigendecomposition shortcut.
%
% Output: figures/long2016_hexapole_halfcut/rbf_rsw.png + temp_code/data/rbf_rsw.mat

LAM  = 1e-4;
RHOS = 20:10:200;                % [um]
RIN  = 150;                      % [um] scoring radius
% [ADDED 2026-09-11] reuse the saved sweep so a figure tweak does not redo 19 solves
REUSE = true;

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
MATF = fullfile(DAT,'rbf_rsw.mat');
if REUSE && isfile(MATF)
    Q = load(MATF);
    if isequal(Q.RHOS(:).', RHOS(:).') && Q.LAM == LAM
        nmae = Q.nmae;  wmax = Q.wmax;  rcn = Q.rcn;  Ni = Q.Ni;  Np = Q.Np;
        REUSE = 'done';
        fprintf(['reused %s' newline], MATF);
    end
end
if ~ischar(REUSE)
S = load(fullfile(DAT,'rbf_w6_r20l0.mat'));
P = S.P;  Np = S.Np;  B6 = S.B6;

in = vecnorm(P,2,2) <= RIN;  Pi = P(in,:);  Ni = nnz(in);
Bi = reshape(B6(in,:,:), Ni, 18);
nb = sum(vecnorm(reshape(Bi,[],3), 2, 2));
RHS = reshape(B6, Np, 18);

s2  = sum(P.^2, 2);
D2n = max(s2 + s2.' - 2*(P*P.'), 0);  D2n(1:Np+1:end) = 0;
D2q = (Pi(:,1)-P(:,1).').^2 + (Pi(:,2)-P(:,2).').^2 + (Pi(:,3)-P(:,3).').^2;

nR = numel(RHOS);  nmae = zeros(nR,1);  wmax = zeros(nR,1);  rcn = zeros(nR,1);
fprintf('lambda = %.0e fixed, Np = %d, scored at %d nodes x 6 coils\n\n', LAM, Np, Ni);
fprintf('%8s %12s %14s %12s\n','rho [um]','max|W|','1/rcond(Phi)','NMAE %');
for k = 1:nR
    rho = RHOS(k);
    Phi = exp(-D2n / rho^2);
    rcn(k) = 1/rcond(Phi + LAM*eye(Np));
    W  = (Phi + LAM*eye(Np)) \ RHS;
    bp = exp(-D2q / rho^2) * W;                       % Ni x 18
    nmae(k) = sum(vecnorm(reshape(bp - Bi, [], 3), 2, 2)) / nb * 100;
    wmax(k) = max(abs(W), [], 'all');
    fprintf('%8d %12.3e %14.3e %12.5f\n', rho, wmax(k), rcn(k), nmae(k));
    clear Phi W bp
end
save(fullfile(DAT,'rbf_rsw.mat'), 'RHOS','nmae','wmax','rcn','LAM','RIN','Ni','Np');
end

% ================================ figure ======================================
FS = 60;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;  CRED = [0.85 0.10 0.10];
XR = [RHOS(1) RHOS(end)];   xt = 65:45:155;   % ends 20 / 200, inner 65/110/155
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]); hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
% [MODIFIED 2026-09-11] markers only: these are 19 discrete evaluations, not a curve.
plot(ax, RHOS, nmae, 'o', 'MarkerFaceColor',CRED, 'MarkerEdgeColor',CRED, ...
     'MarkerSize',18, 'LineStyle','none');

set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
       'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
       'XTickLabelRotation',0,'YTickLabelRotation',0);
xlim(ax, XR);  set(ax,'XTick', xt);
% [MODIFIED 2026-09-11] user: do not force one decimal on the y axis (rule 9 binds the
% horizontal axis only).  ylim_nice_ rounds the step up to the next multiple of its own
% decade unit, which fills the frame far better than the 0.1-multiple rule would.
[yr, yt] = ylim_nice_(max(nmae), 3);
ylim(ax, yr);  set(ax,'YTick', yt);  ytop = yr(2);
for xv = XR
    text(ax, xv, -0.022*ytop, sprintf('%g', xv), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{\rho\;(\mu m)}$', 'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{NMAE\;(\%)}$', 'Interpreter','latex', 'FontSize',FSLAB);
hold(ax,'off');

out = fullfile(FIG,'rbf_rsw.png');
set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');  close(fig);
fprintf(['wrote %s  (y ticks %s, top %g, fill %.0f%%)' newline], out, ...
        mat2str(yt), ytop, 100*max(nmae)/ytop);

function [lim, tk] = ylim_nice_(maxv, N)
% Same shape as ylim_from_zero_ (0 start, N inner ticks, [0,T] split into N+1, the top not
% labelled) but the step is a nice number at ANY decade -- the err_hist version only allows
% multiples of 0.1, which for maxv ~ 0.07 would give 0.1/0.2/0.3 and a 17 % fill.
% Rule 9 (one decimal) binds the horizontal axis only, so 0.02/0.04/0.06 is allowed here.
    if nargin < 2 || isempty(N), N = 3; end
    smin = 1.08*maxv/(N+1);
    u    = 10^floor(log10(smin));           % decade unit of the minimum step
    s    = ceil(smin/u - 1e-12) * u;        % round up to the next multiple of it
    lim  = [0 (N+1)*s];   tk = (1:N)*s;
end
