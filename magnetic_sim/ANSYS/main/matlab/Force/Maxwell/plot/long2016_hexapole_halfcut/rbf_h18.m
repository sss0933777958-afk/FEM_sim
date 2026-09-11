% rbf_h18  Residual histogram: smoothing RBF vs the eighteen-parameter charge model.
%
% Both models predict B at the SAME points -- the 1771 FEM nodes with r <= RIN inside the
% r <= 250 um node cloud -- for all six excitations, so 1771 x 6 = 10626 residuals each.
% Residual per point = || b_pred - B_FEM ||  (vector norm, the project's NMAE convention).
%
%   RBF        : rho = 120 um, lam = 1e-4, weights from rbf_w6.mat (8225 nodes, 18 RHS)
%   eighteen   : b = S(l_hat, Pc) * G,  Pc = Pc_base + e(17-vector), from the Flux
%                calibration calib_current_maxwell_axshN13_R150_eighteen.mat
%
% G's columns are in COIL order (pole-coil-numbering.md), matching B6.
%
% Style copied verbatim from figures/paper_fig_plot/plot/plot_err_hist_shell.m, which
% produced the reference figure err_hist_conv_maxwell_eighteen_R300_zhi.png:
%   red [0.85 0.10 0.10] / blue [0.05 0.10 0.95], FaceAlpha 0.60, bars with no edges,
%   each series self-normalised to 100 %, shared 2.8 uT bins (the err_hist family's one
%   ruler), latex bold axis titles, three inner ticks on both axes with the end points
%   labelled by text only, legend top right.
%
% Output: figures/long2016_hexapole_halfcut/rbf_h18.png + temp_code/data/rbf_h18.mat

RIN   = 150;                                                  % [um]
MATN  = 'calib_current_maxwell_axshN13_R150_eighteen';
BINWUT = 2.8;                                                 % [uT] shared with the err_hist family
WFILE = 'rbf_w6';   % weights to score; the diagnostic run is rbf_w6_r20l0
if exist('WFILE_OVR','var'), WFILE = WFILE_OVR; end
FTAG  = '';  if ~strcmp(WFILE,'rbf_w6'), FTAG = WFILE(7:end); end

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
CAL  = fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(CAL,'function'), fullfile(CAL,'utils'), fullfile(CAL,'common_path'));

S = load(fullfile(DAT,[WFILE '.mat']));
fprintf(['weights ' WFILE ' : rho = %g um, lam = %.0e' newline], S.RHO, S.LAM);
P = S.P;  Np = S.Np;
in = vecnorm(P,2,2) <= RIN;  Pi = P(in,:);  Ni = nnz(in);
Bi = S.B6(in,:,:);                                            % Ni x 3 x 6, mT

% ---- model 1: smoothing RBF ----------------------------------------------------
D2 = (Pi(:,1)-P(:,1).').^2 + (Pi(:,2)-P(:,2).').^2 + (Pi(:,3)-P(:,3).').^2;
GK = exp(-D2 / S.RHO^2);  clear D2
Drbf = zeros(Ni,6);
for k = 1:6
    Drbf(:,k) = vecnorm(GK * S.W6(:,:,k) - Bi(:,:,k), 2, 2);
end

% ---- model 2: eighteen-parameter charge model ----------------------------------
C   = load(fullfile(CAL,'data','long2016_hexapole_halfcut','.mat',[MATN '.mat']));
cfg = model_config('long2016_hexapole_halfcut', C.GEOM);
Pc  = make_Pc_(C.e, cfg.Pc_base);
Sk  = build_S_(C.l_hat, Pc, Pi*1e-6);                         % l_hat is in metres
D18 = zeros(Ni,6);
for k = 1:6
    D18(:,k) = vecnorm(reshape(Sk * C.G(:,k), 3, Ni).' - Bi(:,:,k), 2, 2);
end

nb = sum(vecnorm(reshape(Bi,[],3), 2, 2));
fprintf('points %d x 6 coils = %d | l_hat = %.4f um | NMAE stored in the .mat = %.4f %%\n\n', ...
        Ni, Ni*6, C.l_hat*1e6, C.NMAE);
fprintf('%14s %10s %12s %12s %12s\n','model','NMAE%','mean [mT]','median [mT]','max [mT]');
fprintf('%14s %10.4f %12.5f %12.5f %12.5f\n','RBF', sum(Drbf(:))/nb*100, mean(Drbf(:)), median(Drbf(:)), max(Drbf(:)));
fprintf('%14s %10.4f %12.5f %12.5f %12.5f\n','Eighteen', sum(D18(:))/nb*100, mean(D18(:)), median(D18(:)), max(D18(:)));
fprintf('\nper coil NMAE %% :\n%6s %10s %10s %8s\n','coil','RBF','Eighteen','ratio');
for k = 1:6
    n1 = sum(Drbf(:,k))/sum(vecnorm(Bi(:,:,k),2,2))*100;
    n2 = sum(D18(:,k)) /sum(vecnorm(Bi(:,:,k),2,2))*100;
    fprintf('%6d %10.4f %10.4f %8.2f\n', k, n1, n2, n2/n1);
end
save(fullfile(DAT,['rbf_h18' FTAG '.mat']), 'Drbf','D18','RIN','Ni','MATN','BINWUT','WFILE');

% ================================ figure ======================================
e1 = D18(:);   e2 = Drbf(:);                    % red = eighteen, blue = RBF
cIN  = [0.85 0.10 0.10];   cOUT = [0.05 0.10 0.95];   ALPH = 0.60;
FS = 60;  FSLEG = 45;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;

BINW = BINWUT*1e-3;
eAll = [e1; e2];   maxE = max(eAll);
edg  = 0 : BINW : (ceil(maxE/BINW)*BINW);
ctr  = (edg(1:end-1) + edg(2:end)) / 2;
p1   = histcounts(e1, edg) / numel(e1) * 100;   % each series self-normalised
p2   = histcounts(e2, edg) / numel(e2) * 100;
fprintf('\nbin width %.2f uT -> %d bins (data max %.4f mT)\n', BINWUT, numel(edg)-1, maxE);

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]); hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
h1 = bar(ax, ctr, p1, 1, 'FaceColor',cIN,  'FaceAlpha',ALPH, 'EdgeColor','none');
h2 = bar(ax, ctr, p2, 1, 'FaceColor',cOUT, 'FaceAlpha',ALPH, 'EdgeColor','none');

set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
       'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
       'XTickLabelRotation',0,'YTickLabelRotation',0);

% x: three inner ticks, [0 T] split into four, step prefers multiples of 0.1 (rule 9)
need = prctile_(eAll, 99.5);   got = [];
for mult = [0.1 0.05 0.02]
    k = ceil(1.02*need/(4*mult) - 1e-9);
    T = 4*k*mult;   fill = need/T;
    if isempty(got) || fill > got(3), got = [k*mult T fill mult]; end
    if abs(mult-0.1) < 1e-9 && fill >= 0.60, break; end
end
sX = got(1);   xr = [0 got(2)];   xt = (1:3)*sX;
fprintf('x: three ticks %s, frame [0 %g] mT, fill %.0f%%%s\n', mat2str(round(xt,3)), xr(2), ...
        100*got(3), repmat(' (two decimals, breaks rule 9)', 1, double(got(4) < 0.1-1e-9)));
xlim(ax, xr);   set(ax,'XTick', xt);

[yr, yt] = ylim_from_zero_(max([p1 p2]), 3);
ylim(ax, yr);   set(ax,'YTick', yt);   ytop = yr(2);
for xv = xr
    text(ax, xv, -0.022*ytop, sprintf('%g', xv), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{Residual\;(mT)}$',      'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{Percentage\;(\%)}$',    'Interpreter','latex', 'FontSize',FSLAB);
lg = legend(ax, [h1 h2], {'Eighteen parameters', 'RBF'}, 'Interpreter','tex', ...
            'Location','northeast', 'NumColumns',1);
lg.FontSize = FSLEG;  lg.FontWeight = 'bold';
lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.Color = 'w';
lg.ItemTokenSize = [55 25];
hold(ax,'off');

out = fullfile(FIG,['rbf_h18' FTAG '.png']);
set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');  close(fig);
fprintf('peak: eighteen %.4f %% | RBF %.4f %%\nwrote %s\n', max(p1), max(p2), out);

% ---- local helpers ----------------------------------------------------------------
function Pc = make_Pc_(e17, Pc_base)     % pipeline's local make_Pc
    E = zeros(3,6);
    E(:,1) = e17(1:3);   E(:,2) = e17(4:6);   E(:,3) = e17(7:9);
    E(:,4) = e17(10:12); E(:,5) = e17(13:15);
    E(1,6) = e17(16);    E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);
    Pc = Pc_base + E;
end
function S = build_S_(l_hat, Pc, P)      % pipeline's local build_S
    Np = size(P,1);  pbar = P / l_hat;  S = zeros(3*Np, 6);
    for k = 1:6
        d  = pbar - Pc(:,k).';
        r3 = sum(d.^2, 2).^1.5;
        S(:,k) = reshape((d ./ r3).', 3*Np, 1);
    end
end
function [lim, tk] = ylim_from_zero_(maxv, N)
% From plot_err_hist_shell.m, with one change: when the peak is large the step is
% forced to an INTEGER (user, 2026-09-11) and a multiple of 5 is preferred inside a
% 25 % window, so ticks read 30/60/90 rather than 27.5/55/82.5.  Small peaks keep the
% original 0.1-multiple rule so that 0.4/0.8/1.2 style ticks still appear.
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
        best = k0;
        for k = (k0+1):k1
            if (mod(k,5) == 0 || mod(k,10) == 0) && (k0/k) >= 0.95, best = k;  break; end
        end
        s = best*0.1;
    end
    lim = [0 (N+1)*s];   tk = (1:N)*s;
end
function q = prctile_(x,p)
    x = sort(x(:));  n = numel(x);
    q = interp1((0.5:n-0.5)/n, x, p/100, 'linear', 'extrap');
end
