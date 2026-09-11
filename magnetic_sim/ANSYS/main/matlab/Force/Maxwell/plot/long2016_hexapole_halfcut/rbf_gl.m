% rbf_gl  Gradient along the xa axis for a ladder of (rho, lambda) pairs.
%
% PAIRS lists the combinations to overlay, so one script serves both sweeps:
%   lambda sweep at fixed rho : PAIRS = [20 1e-10; 20 1e-8; 20 1e-6; 20 1e-4]  -> rbf_gl.png
%   rho sweep at fixed lambda : PAIRS = [20 1e-4; 40 1e-4; ...; 200 1e-4]       -> rbf_gr.png
% Phi and the evaluation kernels are rebuilt inside the loop because rho changes them.
%
%   nodes = r <= 250 um (Np = 8225), six excitations solved together (18 RHS); the curve
%   drawn is the P1 (coil 1) excitation, as in every earlier gradient figure.
%
%   analytic gradient, the document's form simplified for the Gaussian:
%       phi'(r)/r = -(2/rho^2) phi(r)  ->  db_c/dx = -(2/rho^2) sum_n W_cn phi(r_n)(x-x_n)
%       d(b.b)/dx = 2 (b_xa db_xa/dx + b_ya db_ya/dx + b_za db_za/dx)
%       F_xa = 0.5 * mgB * UF * d(b.b)/dx        [pN]
%
% Output: figures/long2016_hexapole_halfcut/rbf_gl.png or rbf_gr.png (+ the matching .mat)

PAIRS = [20 1e-4; 40 1e-4; 60 1e-4; 80 1e-4; 120 1e-4; 200 1e-4];
if exist('PAIRS_OVR','var'), PAIRS = PAIRS_OVR; end
COIL = 1;                                    % excitation drawn
NQ   = 1001;  XR = [-150 150];
mgB  = 0.0451;  UF = 1e3;

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
P = S.P;  Np = S.Np;  B6 = S.B6;

t   = linspace(XR(1), XR(2), NQ).';
Q   = [t, zeros(NQ,2)];
s2  = sum(P.^2, 2);
D2n = max(s2 + s2.' - 2*(P*P.'), 0);  D2n(1:Np+1:end) = 0;
D2q = (Q(:,1)-P(:,1).').^2 + (Q(:,2)-P(:,2).').^2 + (Q(:,3)-P(:,3).').^2;
dxn = Q(:,1) - P(:,1).';
RHS = reshape(B6, Np, 18);

nL = size(PAIRS,1);  FF = zeros(NQ,nL);  PP = zeros(nL,1);
fprintf('Np = %d, coil %d (P%d)\n\n', Np, COIL, S.map(COIL));
fprintf('%8s %10s %12s %14s %12s %10s %12s\n', ...
        'rho','lambda','max|W|','max|PhiW-B|','p-p [pN]','% of F0','F(0) [pN]');
for k = 1:nL
    rho = PAIRS(k,1);   lam = PAIRS(k,2);
    Phi = exp(-D2n / rho^2);
    G   = exp(-D2q / rho^2);
    Gd  = -(2/rho^2) * G .* dxn;
    W   = (Phi + lam*eye(Np)) \ RHS;
    Wk  = reshape(W, Np, 3, 6);  Wk = Wk(:,:,COIL);
    b   = G  * Wk;
    db  = Gd * Wk;
    FF(:,k) = mgB * UF * sum(b.*db, 2);          % 0.5*mgB*UF * 2*(b.db)
    r5 = FF(:,k) - polyval(polyfit(t,FF(:,k),5), t);
    PP(k) = max(r5) - min(r5);
    fprintf('%8g %10.0e %12.3e %14.3e %12.4f %10.2f %12.4f\n', rho, lam, ...
            max(abs(Wk),[],'all'), max(abs(Phi*Wk - B6(:,:,COIL)),[],'all'), ...
            PP(k), PP(k)/FF(t==0,k)*100, FF(t==0,k));
    clear Phi G Gd W Wk b db
end
fprintf('\nspread between the curves: max |F - F(row 1)| = %.4f pN\n', ...
        max(abs(FF - FF(:,1)), [], 'all'));

if     all(PAIRS(:,1) == PAIRS(1,1)), STEM = 'rbf_gl';  VARY = 'lam';
elseif all(PAIRS(:,2) == PAIRS(1,2)), STEM = 'rbf_gr';  VARY = 'rho';
else,                                 STEM = 'rbf_gp';  VARY = 'both';  end
save(fullfile(DAT,[STEM '.mat']), 't','FF','PP','PAIRS','COIL','VARY');

% ================================ figure ======================================
FS = 60;  FSLEG = 45;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;  LW = 6;
CM = turbo(nL);
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]); hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
hh = gobjects(1,nL);  lbl = cell(1,nL);
for k = 1:nL
    hh(k) = plot(ax, t, FF(:,k), '-', 'Color',CM(k,:), 'LineWidth',LW);
    switch VARY
    case 'lam', lbl{k} = ['\lambda = 10^{' num2str(round(log10(PAIRS(k,2)))) '}'];
    case 'rho', lbl{k} = ['\rho = ' num2str(PAIRS(k,1)) ' \mum'];
    otherwise,  lbl{k} = ['\rho=' num2str(PAIRS(k,1)) ', \lambda=10^{' ...
                          num2str(round(log10(PAIRS(k,2)))) '}'];
    end
end
set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
       'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
       'XTickLabelRotation',0,'YTickLabelRotation',0);
sx = diff(XR)/4;
xlim(ax, XR);  set(ax,'XTick', XR(1)+(1:3)*sx);
[yr, yt] = ylim_nice_(max(FF(:)), 3);
ylim(ax, yr);  set(ax,'YTick', yt);  ytop = yr(2);
for xv = XR
    text(ax, xv, -0.022*ytop, sprintf('%g', xv), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{x_a\;(\mu m)}$',  'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{F_{x_a}\;(pN)}$', 'Interpreter','latex', 'FontSize',FSLAB);
lg = legend(ax, hh, lbl, 'Interpreter','tex', 'Location','northwest', 'NumColumns',1);
lg.FontSize = FSLEG*0.85;  lg.FontWeight = 'bold';
lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.Color = 'w';
lg.ItemTokenSize = [50 22];
hold(ax,'off');

out = fullfile(FIG,[STEM '.png']);
set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');  close(fig);
fprintf(['wrote %s  (y ticks %s, top %g, fill %.0f%%)' newline], out, ...
        mat2str(yt), ytop, 100*max(FF(:))/ytop);

function [lim, tk] = ylim_nice_(maxv, N)
    if nargin < 2 || isempty(N), N = 3; end
    smin = 1.08*maxv/(N+1);
    u    = 10^floor(log10(smin));
    s    = ceil(smin/u - 1e-12) * u;
    lim  = [0 (N+1)*s];   tk = (1:N)*s;
end
