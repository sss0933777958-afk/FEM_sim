% rbf_lsw  lambda sweep at fixed rho: NMAE of the predicted field vs lambda.
%
%   rho = 20 um (node spacing), nodes r <= 250 um (Np = 8225, six excitations).
%   lambda = 0 : DLAM : LMAX, linear, so the sweep starts at the exact interpolant.
%   Scored at the 1771 nodes with r <= 150 um, all six coils:
%       NMAE = sum_i || b_pred,i - B_i ||  /  sum_i || B_i ||  * 100
%
% All lambdas share ONE eigendecomposition of Phi:
%       W(lam) = V * ( (V' * B) ./ (d + lam) )        <=>   (Phi + lam I) W = B
% so 10001 values cost one eig plus one small matrix product each.
% (Verified against backslash and inv() to 1e-10 relative, 2026-09-11.)
%
% rho = 20 is the right place to start a sweep at lam = 0: min eig(Phi) = 0.029 > 0,
% so the exact interpolant is well posed there (at rho = 120 it is numerically singular).
%
% Output: figures/long2016_hexapole_halfcut/rbf_lsw.png + temp_code/data/rbf_lsw.mat

RHO  = 20;                       % [um]
LMAX = 10;    DLAM = 0.001;      % lambda grid
RIN  = 150;                      % [um] scoring radius
XMAX = 1;                        % [--] plot range; the curve fills the frame
% [ADDED 2026-09-11] reuse the saved sweep instead of redoing the eig + 10001 solves
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
MATF = fullfile(DAT,'rbf_lsw.mat');
if REUSE && isfile(MATF)
    Q = load(MATF);
    if Q.RHO == RHO && abs(Q.DLAM - DLAM) < 1e-12 && Q.LMAX >= XMAX
        lam = Q.lam;  nmae = Q.nmae;  Ni = Q.Ni;  Np = Q.Np;
        fprintf(['reused %s (%d lambdas, rho = %g)' newline], MATF, numel(lam), RHO);
        REUSE = 'done';
    end
end
if ~ischar(REUSE)
S = load(fullfile(DAT,'rbf_w6_r20l0.mat'));
P = S.P;  Np = S.Np;  B6 = S.B6;
assert(S.RHO == RHO, 'weights file was built at rho = %g', S.RHO);

in = vecnorm(P,2,2) <= RIN;  Pi = P(in,:);  Ni = nnz(in);
Bi = reshape(B6(in,:,:), Ni, 18);
nb = sum(vecnorm(reshape(Bi,[],3), 2, 2));          % denominator of NMAE

% ---- Phi and its eigendecomposition (once) --------------------------------------
s2  = sum(P.^2, 2);
D2n = max(s2 + s2.' - 2*(P*P.'), 0);  D2n(1:Np+1:end) = 0;
Phi = exp(-D2n / RHO^2);  Phi = (Phi + Phi.')/2;
D2q = (Pi(:,1)-P(:,1).').^2 + (Pi(:,2)-P(:,2).').^2 + (Pi(:,3)-P(:,3).').^2;
Gs  = exp(-D2q / RHO^2);                             % Ni x Np evaluation kernel
clear D2n D2q
tic;  [V, Dg] = eig(Phi);  d = diag(Dg);  te = toc;
fprintf('eig(Phi): %.1f s | eigenvalues %.4e .. %.4e (cond %.3e)\n', te, min(d), max(d), max(d)/min(d));
clear Phi Dg
A = Gs * V;                                          % Ni x Np
Cm = V.' * reshape(B6, Np, 18);                      % Np x 18
clear V Gs

% ---- sweep ----------------------------------------------------------------------
lam  = (0:DLAM:LMAX).';   nL = numel(lam);
nmae = zeros(nL,1);
tic;
for j = 1:nL
    bp = A * (Cm ./ (d + lam(j)));                   % Ni x 18
    nmae(j) = sum(vecnorm(reshape(bp - Bi, [], 3), 2, 2)) / nb * 100;
end
fprintf('swept %d lambdas in %.1f s\n\n', nL, toc);

save(fullfile(DAT,'rbf_lsw.mat'), 'lam','nmae','RHO','RIN','Ni','Np','DLAM','LMAX');

end
fprintf('%10s %12s\n','lambda','NMAE %');
for v = [0 0.001 0.01 0.05 0.1 0.25 0.5 0.75 1]
    [~,j] = min(abs(lam - v));
    fprintf('%10.3f %12.4f\n', lam(j), nmae(j));
end
fprintf('\nmonotone increasing: %d (largest step down %.3e)\n', all(diff(nmae) >= 0), min(diff(nmae)));
for t = [1 5 10 25 50]
    j = find(nmae >= t, 1);
    if ~isempty(j), fprintf('NMAE first reaches %5.1f %% at lambda = %.3f\n', t, lam(j)); end
end

% ================================ figure ======================================
FS = 60;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;  CBLU = [0.05 0.10 0.95];
% [MODIFIED 2026-09-11] the frame now ENDS where the data ends, so the curve spans
% the whole view.  Four inner ticks 0.2/0.4/0.6/0.8: equally spaced, the gaps to both
% end points are the same 0.2, and every label has one decimal (rule 9).
kp = lam <= XMAX;
XR = [0 XMAX];   xt = (0.2:0.2:0.8) * XMAX;
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]); hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
% [MODIFIED 2026-09-11] markers, not a line: these are swept evaluations.  The grid is
% 0.001, i.e. 2 px apart on this canvas, so it is decimated to ~41 markers (lambda every
% 0.025) -- otherwise the markers merge into a solid band.
iv = find(kp);   DEC = max(1, round(numel(iv)/41));   iv = iv(1:DEC:end);
plot(ax, lam(iv), nmae(iv), 'o', 'MarkerFaceColor',CBLU, 'MarkerEdgeColor',CBLU, ...
     'MarkerSize',16, 'LineStyle','none');

set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
       'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
       'XTickLabelRotation',0,'YTickLabelRotation',0);
xlim(ax, XR);  set(ax,'XTick', xt);
[yr, yt] = ylim_from_zero_(max(nmae(kp)), 5);   % five ticks -> tighter fill
ylim(ax, yr);  set(ax,'YTick', yt);  ytop = yr(2);
for xv = XR
    text(ax, xv, -0.022*ytop, sprintf('%g', xv), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{\lambda}$',            'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{NMAE\;(\%)}$',         'Interpreter','latex', 'FontSize',FSLAB);
hold(ax,'off');

out = fullfile(FIG,'rbf_lsw.png');
set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');  close(fig);
fprintf('\nwrote %s  (y ticks %s, top %g, fill %.0f%%)\n', out, mat2str(yt), ytop, 100*max(nmae(kp))/ytop);

function [lim, tk] = ylim_from_zero_(maxv, N)
    if nargin < 2 || isempty(N), N = 3; end
    smin = 1.08*maxv/(N+1);
    if smin >= 1
        k0 = ceil(smin - 1e-9);   s = k0;
        for k = k0:ceil(1.25*smin)
            if mod(k,10) == 0 || mod(k,5) == 0, s = k;  break; end
        end
    else
        k0 = max(1, ceil(smin/0.1 - 1e-9));
        k1 = max(k0, ceil(1.15*smin/0.1));
        s  = k0*0.1;
        for k = (k0+1):k1
            if (mod(k,5) == 0 || mod(k,10) == 0) && (k0/k) >= 0.95, s = k*0.1;  break; end
        end
    end
    lim = [0 (N+1)*s];   tk = (1:N)*s;
end
