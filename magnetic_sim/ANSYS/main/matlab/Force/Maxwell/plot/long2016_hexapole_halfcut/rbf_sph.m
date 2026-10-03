% rbf_sph  Test-set residual histogram: RBF vs the solid-harmonic model.
%
%   Both models are fitted on the SAME 80 % of the FEM nodes inside r <= 150 um and
%   scored on the same held-out 20 % (354 nodes x 6 excitations = 2124 residuals).
%   Residual per point = || b_pred - B_FEM ||, the project's vector-norm convention.
%
%       RBF   rho = 400 um, lam = 4.64e-7 -- the parameter pair that minimises the
%             test NMAE AMONG those whose gradient passes the smoothness test
%             (d3(b.b) does not change sign on any of the three actuator axes).
%             The unconstrained NMAE optimum (rho = 80) is more accurate but its
%             d3 changes sign 7 to 9 times per axis, so it cannot be used for force.
%       sph   L = 7 (K = 63 coefficients) -- the lowest test NMAE of the L sweep,
%             and smooth on all three axes.
%
%   Reads only (plot-scripts-pure); rbf_sph.mat is written by the fitting run.
%
% Output: figures/long2016_hexapole_halfcut/rbf_sph.png

HERE = fileparts(mfilename('fullpath'));                      % .../plot/<model>
FMX  = fileparts(fileparts(HERE));                            % .../Force/Maxwell
DAT  = fullfile(FMX,'utils','data');
FIG  = fullfile(FMX,'figures','long2016_hexapole_halfcut','current');

S = load(fullfile(DAT, 'rbf_sph.mat'));
e1 = S.D_rbf(:);   e2 = S.D_sph(:);                           % red = RBF, blue = harmonics
fprintf(['RBF  rho = %g um, lam = %.2e : NMAE %.4f %% | mean %.5f mT' newline], S.rho, S.lam, S.nmr, mean(e1));
fprintf(['sph  L = %d (K = %d)          : NMAE %.4f %% | mean %.5f mT' newline], S.Lsph, (S.Lsph+1)^2-1, S.nms, mean(e2));

% ================================ figure ======================================
FS = 60;  FSLEG = 45;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;
cR = [0.85 0.10 0.10];   cS = [0.05 0.10 0.95];   ALPH = 0.60;

BINW = S.BINW;
edg = 0 : BINW : (ceil(max([e1;e2])/BINW)*BINW);
ctr = (edg(1:end-1) + edg(2:end)) / 2;
p1  = histcounts(e1, edg) / numel(e1) * 100;                  % each series self-normalised
p2  = histcounts(e2, edg) / numel(e2) * 100;

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]); hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
h1 = bar(ax, ctr, p1, 1, 'FaceColor',cR, 'FaceAlpha',ALPH, 'EdgeColor','none');
h2 = bar(ax, ctr, p2, 1, 'FaceColor',cS, 'FaceAlpha',ALPH, 'EdgeColor','none');

set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
       'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
       'XTickLabelRotation',0,'YTickLabelRotation',0);

[xr, xt] = ticks3_(0, max([e1;e2]));                          % three inner ticks, x
xlim(ax, xr);   set(ax,'XTick', xt);
[yr, yt] = ticks3_(0, max([p1 p2]));                          % three inner ticks, y
ylim(ax, yr);   set(ax,'YTick', yt);
for xv = xr                                                   % ends: numbers only (rule 4)
    text(ax, xv, yr(1)-0.022*diff(yr), sprintf('%g', xv), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{Residual\;(mT)}$',   'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{Percentage\;(\%)}$', 'Interpreter','latex', 'FontSize',FSLAB);

% Same legend styling as the rbf_f18_* figures: bold tex, framed, box line as thick
% as the axes.  It sits top right because the bars occupy the left of the frame.
lg = legend(ax, [h1 h2], {'RBF', 'Spherical harmonics'}, 'Interpreter','tex', ...
            'Location','northeast', 'NumColumns',1);
lg.FontSize = FSLEG;  lg.FontWeight = 'bold';
lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.Color = 'w';
lg.ItemTokenSize = [55 25];
hold(ax,'off');

out = fullfile(FIG, 'rbf_sph.png');
set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');   close(fig);
fprintf(['x ticks %s (frame [%g %g]) | y ticks %s (frame [%g %g])' newline], ...
        mat2str(xt), xr, mat2str(yt), yr);
fprintf(['peak: RBF %.2f %% | harmonics %.2f %%' newline 'wrote %s' newline], max(p1), max(p2), out);

% ---- local helper ------------------------------------------------------------
function [lim, tk] = ticks3_(lo, hi)
% Three inner ticks, equally spaced, with the gaps to both frame edges equal to the
% tick spacing (rules 4 and 5): frame = [lo, lo+4s], ticks = lo + (1:3)*s, with s the
% smallest nice step that covers the data.
    need = (hi - lo)/4;
    for kk = floor(log10(need)) : floor(log10(need))+1
        for m = [1 1.5 2 2.5 3 4 5 6 8]
            s = m*10^kk;
            if s >= need*(1+1e-9)
                lim = [lo, lo + 4*s];   tk = lo + (1:3)*s;   return
            end
        end
    end
    error('ticks3_:none', 'no nice step for [%g %g]', lo, hi);
end
