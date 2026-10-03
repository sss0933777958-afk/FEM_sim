%% pred3.m -- the two R <= 150 um models predicting b at the trilinear sample points
%%             on all three actuator axes, all six excitations.
%
%  REFERENCE.  xa_interp.mat: 1001 equidistant points per axis (3 x 1001 = 3003
%  points), b taken from the FEM export by the project's regular-grid trilinear, in
%  the actuator frame, six excitations.  None of the 3003 points sits on a lattice
%  node (checked there), so every reference value is a genuine interpolation.
%
%  MODELS.  Both fitted on the 1771 nodes inside R = 150 um, at the setting that
%  passes the smoothness gate with the lowest working-region NMAE:
%      RBF   rho = 490 um, lam = 1e-8      (rbf_v150.mat)
%      sph   L = 9, K = 99                 (sph_v150.mat)
%
%  So each model is asked for 6 x 3003 = 18018 field vectors.
%
%      NMAE = sum ||b_model - b_ref|| / sum ||b_ref|| * 100
%
%  pooled over points and excitations, plus the per-axis and per-excitation splits.
%  The histogram shows the pointwise relative error ||db|| / ||b_ref|| * 100 over all
%  18018 (point, excitation) pairs, one curve per model.
%
%  Output: data/long2016_hexapole_halfcut/.mat/pred3.mat
%          figures/long2016_hexapole_halfcut/pred3_hist.png     [ADDED 2026-09-18]

clearvars;  clc;
HERE = fileparts(mfilename('fullpath'));  FMX = fileparts(fileparts(HERE));
MAIN = fileparts(fileparts(fileparts(FMX)));
FLX  = fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT = fullfile(FMX,'utils','data');
FIG = fullfile(FMX,'figures','long2016_hexapole_halfcut','current');

R0 = 150;   ANM = {'x_a','y_a','z_a'};

%% ---- the reference and the two models ---------------------------------------
X  = load(fullfile(DAT,'xa_interp.mat'));
NQ = X.NQ;
SR = load(fullfile(DAT,'rbf_v150.mat'));
[~,iw] = min(SR.TAB.nmae);   RHO = SR.TAB.rho(iw);   LAM = SR.TAB.lam(iw);
SS = load(fullfile(DAT,'sph_v150.mat'));   r = SS.RES([SS.RES.R]==R0);   LDEG = r.L;
fprintf(['RBF : rho = %g um, lam = %.0e   (fit NMAE on the 1771 nodes %.4f %%)' newline], ...
        RHO, LAM, SR.TAB.nmae(iw));
fprintf(['sph : L = %d (K = %d)             (fit NMAE on the 1771 nodes %.4f %%)' newline], ...
        LDEG, r.K, r.nmae);

[bR,~,iN] = rbf_field(R0, struct('rho',RHO,'lam',LAM));
assert(iN.Np == 1771, 'pred3:nodes', 'expected 1771 nodes, got %d', iN.Np);
bS = sph_field(R0, struct('P',iN.P,'B',iN.B6), ...
               struct('L',LDEG,'Rn',R0,'quiet',true,'condmax',1500));

%% ---- predict, and stack every (axis, point, excitation) ----------------------
NC = 6;
REF = [];   PRD_R = [];   PRD_S = [];   AXID = [];   EXID = [];
for a = 1:3
    Pq = X.P_act{a};                                    % NQ x 3, actuator frame [um]
    for j = 1:NC
        I = zeros(6,1);  I(j) = 1;
        REF   = [REF;   X.B{a}(:,:,j)];                 %#ok<AGROW>
        PRD_R = [PRD_R; bR(Pq, I)];                     %#ok<AGROW>
        PRD_S = [PRD_S; bS(Pq, I)];                     %#ok<AGROW>
        AXID  = [AXID;  repmat(a, NQ, 1)];              %#ok<AGROW>
        EXID  = [EXID;  repmat(j, NQ, 1)];              %#ok<AGROW>
    end
end
fprintf([newline 'stacked %d field vectors per model (%d points x %d excitations)' newline], ...
        size(REF,1), 3*NQ, NC);
assert(size(REF,1) == 3*NQ*NC, 'pred3:count', 'expected %d rows', 3*NQ*NC);

nref = vecnorm(REF, 2, 2);
eR   = vecnorm(PRD_R - REF, 2, 2);
eS   = vecnorm(PRD_S - REF, 2, 2);
NMAE_R = sum(eR)/sum(nref)*100;
NMAE_S = sum(eS)/sum(nref)*100;
relR = eR ./ nref * 100;   relS = eS ./ nref * 100;     % pointwise, per cent

fprintf([newline 'NMAE over all %d vectors :  RBF %.4f %%   sph %.4f %%' newline], ...
        numel(nref), NMAE_R, NMAE_S);
fprintf([newline '%6s %12s %12s' newline], 'axis', 'RBF %', 'sph %');
NMA = zeros(2,3);
for a = 1:3
    k = AXID == a;
    NMA(1,a) = sum(eR(k))/sum(nref(k))*100;   NMA(2,a) = sum(eS(k))/sum(nref(k))*100;
    fprintf('%6s %12.4f %12.4f\n', ANM{a}, NMA(1,a), NMA(2,a));
end
fprintf([newline '%6s %12s %12s' newline], 'exc', 'RBF %', 'sph %');
NME = zeros(2,NC);
for j = 1:NC
    k = EXID == j;
    NME(1,j) = sum(eR(k))/sum(nref(k))*100;   NME(2,j) = sum(eS(k))/sum(nref(k))*100;
    fprintf('%6s %12.4f %12.4f\n', sprintf('P%d',j), NME(1,j), NME(2,j));
end
fprintf([newline 'pointwise relative error: RBF median %.4f, p95 %.4f, max %.4f %%' newline], ...
        median(relR), prctile_(relR,95), max(relR));
fprintf(['                          sph median %.4f, p95 %.4f, max %.4f %%' newline], ...
        median(relS), prctile_(relS,95), max(relS));

save(fullfile(DAT,'pred3.mat'), 'REF','PRD_R','PRD_S','AXID','EXID','relR','relS', ...
     'NMAE_R','NMAE_S','NMA','NME','RHO','LAM','LDEG','R0','NQ','ANM');

%% ---- the histogram, styled exactly like plot/long2016_hexapole_halfcut/rbf_sph.m
%  Filled translucent bars (the two overlap into purple where they agree), x = the
%  residual in mT, not a percentage -- same figure the earlier test-set histogram
%  used, so the two can be read side by side.
FS = 60;  FSLEG = 45;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;
cR = [0.85 0.10 0.10];   cS = [0.05 0.10 0.95];   ALPH = 0.60;

% BINW matches the reference figure (rbf_sph.mat: 0.001 mT).  Freedman-Diaconis on
% these 18018 residuals asks for 0.0003 mT, but 149 bars over the same [0, 0.06]
% frame read as noise next to the 55 of the reference, and the two figures are meant
% to be compared side by side -- so the bin is pinned, not derived.
BINW = 0.001;
edg = 0 : BINW : (ceil(max([eR;eS])/BINW)*BINW);
ctr = (edg(1:end-1) + edg(2:end)) / 2;
p1  = histcounts(eR, edg) / numel(eR) * 100;             % each series self-normalised
p2  = histcounts(eS, edg) / numel(eS) * 100;
fprintf([newline 'histogram: bin %.4f mT, %d bins, peak RBF %.2f %% / sph %.2f %%' newline], ...
        BINW, numel(p1), max(p1), max(p2));
fprintf(['residual: RBF mean %.5f, median %.5f, max %.5f mT' newline], ...
        mean(eR), median(eR), max(eR));
fprintf(['          sph mean %.5f, median %.5f, max %.5f mT' newline], ...
        mean(eS), median(eS), max(eS));

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]); hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
h1 = bar(ax, ctr, p1, 1, 'FaceColor',cR, 'FaceAlpha',ALPH, 'EdgeColor','none');
h2 = bar(ax, ctr, p2, 1, 'FaceColor',cS, 'FaceAlpha',ALPH, 'EdgeColor','none');

set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
       'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
       'XTickLabelRotation',0,'YTickLabelRotation',0);
[xr, xt] = ticks3_(0, max([eR;eS]));
xlim(ax, xr);   set(ax,'XTick', xt);
[yr, yt] = ticks3_(0, max([p1 p2]));
ylim(ax, yr);   set(ax,'YTick', yt);
for xv = xr                                              % ends: numbers only (rule 4)
    text(ax, xv, yr(1)-0.022*diff(yr), sprintf('%g', xv), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{Residual\;(mT)}$',   'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{Percentage\;(\%)}$', 'Interpreter','latex', 'FontSize',FSLAB);
lg = legend(ax, [h1 h2], {'RBF', 'Spherical harmonics'}, 'Interpreter','tex', ...
            'Location','northeast', 'NumColumns',1);
lg.FontSize = FSLEG;  lg.FontWeight = 'bold';
lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.Color = 'w';
lg.ItemTokenSize = [55 25];
hold(ax,'off');
out = fullfile(FIG,'pred3_hist.png');
set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');   close(fig);
fprintf(['x ticks %s (frame [%g %g]) | y ticks %s (frame [%g %g])' newline], ...
        mat2str(xt), xr, mat2str(yt), yr);
fprintf(['wrote %s' newline], out);

%% ---- local helpers ----------------------------------------------------------
function y = prctile_(v, p)
    v = sort(v(:));   y = v(max(1, round(p/100*numel(v))));
end

function y = iqr_(v)
    y = prctile_(v,75) - prctile_(v,25);
end

function [lim, tk] = ticks3_(lo, hi)
% Three inner ticks, equally spaced, gaps to both frame edges equal to the spacing
% (rules 4 and 5) -- copied from plot/long2016_hexapole_halfcut/rbf_sph.m.
    need = (hi - lo)/4;
    for kk = floor(log10(need)) : floor(log10(need))+1
        for m = [1 1.5 2 2.5 3 4 5 6 8]
            st = m*10^kk;
            if st >= need*(1+1e-9)
                lim = [lo, lo + 4*st];   tk = lo + (1:3)*st;   return
            end
        end
    end
    error('ticks3_:none', 'no nice step for [%g %g]', lo, hi);
end
