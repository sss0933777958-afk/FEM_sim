%% x200.m -- RBF fitted on R <= 200 um against the R <= 150 um spherical harmonics:
%%            the gradient on the three actuator axes, and their difference.
%
%  The two models are deliberately NOT on the same data range (user's call):
%      RBF   fitted on every node inside R = 200 um, at the (rho,lam) that passes
%            the smoothness gate there  -- read from rbf_v200.mat
%      sph   unchanged from the R = 150 um study, L = 9  -- read from sph_v150.mat
%  Both are drawn over |s| <= 150 um, the working region, so the RBF is being asked
%  what its extra ring of data does to the gradient INSIDE the region.
%
%  The eighteen-parameter curve is left out of both figures (user's call).
%
%  SIGNS.  z_a is the P5 direction and the excited pole is P6, which sits at
%  NEGATIVE z_a, so that panel's gradient is negated -- all three panels then read
%  "rises towards the excited pole".  The DIFFERENCE panels carry the same negation,
%  so RBF - sph keeps its meaning across the three.
%
%  Output: data/long2016_hexapole_halfcut/.mat/x200.mat
%          figures/long2016_hexapole_halfcut/x200_g3.png    gradient, both models
%          figures/long2016_hexapole_halfcut/x200_dif.png   RBF - sph
%
%  [ADDED 2026-09-17]

clearvars;  clc;
HERE = fileparts(mfilename('fullpath'));
FMX  = fileparts(fileparts(HERE));
MAIN = fileparts(fileparts(fileparts(FMX)));
FLX  = fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT = fullfile(FMX,'utils','data');
FIG = fullfile(FMX,'figures','long2016_hexapole_halfcut','current');

% [ADDED 2026-09-18 user's call] FIXRHO pins the RBF to the R = 150 setting instead
% of the gate's own pick at R = 200 (rho 340), so the only difference between the two
% RBF fits is the data range.  Checked against rbf_v200.mat: (490, 1e-8) DOES pass the
% smoothness gate at R = 200 -- the gate is in fact looser there (662 of 1188 pairs
% pass vs 501 at R = 150, and the smallest passing rho drops from 340 to 170).
% Outputs take a _fix suffix so the per-range-pick versions survive.
FIXRHO = true;   FIX_RHO = 490;   FIX_LAM = 1e-8;
SUF = '';   if FIXRHO, SUF = '_fix';  end

R_RBF = 200;   R_SPH = 150;
SPAN  = 150;   NQ = 1001;
POLE  = [1 3 6];   DSGN = [1 1 -1];   ANM = {'x_a','y_a','z_a'};

%% ---- the two settings --------------------------------------------------------
SR = load(fullfile(DAT,'rbf_v200.mat'));
if FIXRHO
    RHO = FIX_RHO;   LAM = FIX_LAM;
    k = abs(SR.TAB.rho-RHO) < 1e-9 & abs(SR.TAB.lam-LAM) < 1e-18;
    assert(any(k) && SR.TAB.pass(k), 'x200:gate', ...
           'rho = %g, lam = %.0e does not pass the gate at R = %d', RHO, LAM, R_RBF);
    fprintf(['RBF @ R = %d : rho = %g um, lam = %.0e (PINNED, gate passes), NMAE = %.4f %%' newline], ...
            R_RBF, RHO, LAM, SR.TAB.nmae(k));
else
    [~,iw] = min(SR.TAB.nmae);   RHO = SR.TAB.rho(iw);   LAM = SR.TAB.lam(iw);
    fprintf(['RBF @ R = %d : %d of %d pairs pass;  rho = %g um, lam = %.0e, NMAE = %.4f %%' newline], ...
            R_RBF, nnz(SR.TAB.pass), numel(SR.TAB.rho), RHO, LAM, SR.TAB.nmae(iw));
end
SS = load(fullfile(DAT,'sph_v150.mat'));
r  = SS.RES([SS.RES.R] == R_SPH);   LDEG = r.L;
fprintf(['sph @ R = %d : L = %d (K = %d), NMAE = %.4f %%' newline], R_SPH, LDEG, r.K, r.nmae);

%% ---- fit ---------------------------------------------------------------------
[~, gR, iR] = rbf_field(R_RBF, struct('rho',RHO, 'lam',LAM));
% sph_field is handed its nodes; take them from a plain R = 150 build so it keeps
% the identical 1771-node set it was selected on.
[~, ~, i150] = rbf_field(R_SPH, struct('rho',RHO, 'lam',LAM));
[~, gS, iS]  = sph_field(R_SPH, struct('P',i150.P, 'B',i150.B6), ...
                         struct('L',LDEG, 'Rn',R_SPH, 'quiet',true));
fprintf(['RBF on %d nodes, sph on %d nodes' newline], iR.Np, iS.Np);

%% ---- the three axes ----------------------------------------------------------
s  = linspace(-SPAN, SPAN, NQ).';
GR = zeros(NQ,3);   GS = zeros(NQ,3);
for a = 1:3
    Pq = zeros(NQ,3);   Pq(:,a) = s;
    I  = zeros(6,1);    I(POLE(a)) = 1;
    G = gR(Pq, I);        GR(:,a) = DSGN(a) * G(:,a);
    G = gS(Pq, I, 'du');  GS(:,a) = DSGN(a) * G(:,a);
end
D = GR - GS;                                   % the difference that gets its own figure

%% ---- mean error --------------------------------------------------------------
%     PCT(a) = sum_i |grad_RBF - grad_SPH| / sum_i grad_SPH * 100
% a ratio of totals, the same form as the range-comparison table, with the harmonics
% as the reference.  PCTPT is the mean of the plotted pointwise curve -- a different
% statistic, kept alongside because the two are easy to confuse.
fprintf([newline '%6s %14s %14s %12s %12s' newline], ...
        'axis','mean|RBF-sph|','max|RBF-sph|','sum/sum %','mean ptwise %');
mabs = zeros(1,3);   PCT = zeros(1,3);   PCTPT = zeros(1,3);
for a = 1:3
    mabs(a)  = mean(abs(D(:,a)));
    PCT(a)   = sum(abs(D(:,a))) / sum(GS(:,a)) * 100;
    PCTPT(a) = mean(abs(D(:,a) ./ GS(:,a))) * 100;
    fprintf('%6s %14.3e %14.3e %11.3f %% %11.3f %%\n', ANM{a}, mabs(a), ...
            max(abs(D(:,a))), PCT(a), PCTPT(a));
end
mean_error = mean(abs(D(:)));
PCTA = sum(abs(D(:))) / sum(GS(:)) * 100;
fprintf([newline 'all three axes : mean|RBF-sph| = %.4e mT^2/um,  sum/sum = %.3f %%' newline], ...
        mean_error, PCTA);
msig = arrayfun(@(a) mean(D(:,a)), 1:3);   mrel = PCT;

save(fullfile(DAT,['x200' SUF '.mat']), 's','GR','GS','D','mean_error','mabs','msig','mrel', ...
     'RHO','LAM','LDEG','R_RBF','R_SPH','SPAN','NQ','POLE','DSGN','ANM');

%% ---- figures -----------------------------------------------------------------
FS = 44;  FSLAB = 32;  FSLEG = 45;  CANV = 14.5;  LWBOX = 5.0;  LW = 7;
CR = [0.85 0.10 0.10];   CS = [0.05 0.10 0.95];   CD = [0 0 0];
XR = [-SPAN SPAN];   XT = [-SPAN/2 0 SPAN/2];
% [MODIFIED 2026-09-17 user's call] every panel carries its own x tick numbers and
% its own axis title (x_a / y_a / z_a), overriding rule 8.  The gap grows to 0.09 to
% hold both, so FSLAB drops to 26 to keep the two-line y label inside its frame.
% GAP is 0.12, not the 0.10 that looked sufficient on paper.  TickLength is measured
% against the LONGER axis, so on a panel 10 in wide and 2.4 in tall the outward ticks
% on the frame stick out about 0.16 in -- the ticks along the TOP edge of the panel
% below reach up into the gap and were being struck by the title sitting above them.
H  = 0.1617;   GAP = 0.1200;   Y0 = 0.1550;   YLG = 0.9750;   FSLAB = 26;

% ---- (a) the two gradients ---------------------------------------------------
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
for a = 1:3
    yb = Y0 + (3-a)*(H+GAP);
    ax = axes('Parent',fig,'Position',[0.1895 yb 0.7155 H]);  hold(ax,'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    h1 = plot(ax, s, GR(:,a), '-',  'Color', CR, 'LineWidth', LW);
    h2 = plot(ax, s, GS(:,a), '--', 'Color', CS, 'LineWidth', LW-3);
    style_(ax, FS, LWBOX);
    [YRa, YTa] = yaxis3_(min([GR(:,a);GS(:,a)]), max([GR(:,a);GS(:,a)]));
    xlim(ax, XR);  set(ax,'XTick',XT);   ylim(ax, YRa);  set(ax,'YTick',YTa);
    ylabel(ax, {['$\mathbf{d(b\cdot b)/d' ANM{a} '}$'], '$\mathbf{(mT^{2}/\mu m)}$'}, ...
           'Interpreter','latex','FontSize',FSLAB);
    if a == 1
        lg = leg_(ax, [h1 h2], {'RBF','Spherical harmonics'}, FSLEG, LWBOX);
        ax.Position = [0.1895 yb 0.7155 H];
    end
    xends_(ax, XR, YRa, FS);
    xttl_(ax, YRa, FSLAB, H*CANV, ANM{a});
    hold(ax,'off');
end
print_(fig, fullfile(FIG,['x200_g3' SUF '.png']), CANV);

% ---- (b) the difference ------------------------------------------------------
% One series per panel, so no legend (figure-style.md): the y label says what it is.
H2 = 0.1817;   GAP2 = 0.1200;   Y02 = Y0;   % no legend here: taller panels
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
for a = 1:3
    yb = Y02 + (3-a)*(H2+GAP2);
    ax = axes('Parent',fig,'Position',[0.1895 yb 0.7155 H2]);  hold(ax,'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    % [MODIFIED 2026-09-18 user's call] plotted as a PERCENTAGE, the pointwise
    % relative difference against the harmonics' own gradient:
    %     100 * ( grad_RBF(p) - grad_SPH(p) ) / grad_SPH(p)
    % Sign kept, so the shape reads as before.  The harmonics are the denominator
    % only because something has to be -- neither model is truth, so this is a
    % relative DIFFERENCE, not RBF's error.  grad_SPH never approaches zero here.
    Dp = 100 * D(:,a) ./ GS(:,a);
    plot(ax, s, Dp, '-', 'Color', CD, 'LineWidth', LW);
    style_(ax, FS, LWBOX);
    [YRa, YTa] = yaxis3_(min(Dp), max(Dp));
    xlim(ax, XR);  set(ax,'XTick',XT);   ylim(ax, YRa);  set(ax,'YTick',YTa);
    ax.YAxis.Exponent = 0;                     % a percentage needs no shared factor
    ylabel(ax, {['$\mathbf{(\nabla_{RBF}-\nabla_{SPH})/\nabla_{SPH}}$'], ...
                ['$\mathbf{on\;' ANM{a} '\;(\%)}$']}, ...
           'Interpreter','latex','FontSize',FSLAB);
    xends_(ax, XR, YRa, FS);
    xttl_(ax, YRa, FSLAB, H2*CANV, ANM{a});
    hold(ax,'off');
end
print_(fig, fullfile(FIG,['x200_dif' SUF '.png']), CANV);

%% ---- local helpers -----------------------------------------------------------
function style_(ax, FS, LWBOX)
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
           'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
           'XTickLabelRotation',0,'YTickLabelRotation',0);
end

function lg = leg_(ax, h, names, FSLEG, LWBOX)
% Built at a normal location and then MOVED: the 'outside' locations resize the axes
% they belong to, and the three panels would stop matching.
    lg = legend(ax, h, names, 'Interpreter','tex', 'Location','northwest', 'NumColumns',2);
    lg.FontSize = FSLEG;   lg.FontWeight = 'bold';
    lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.Color = 'w';
    lg.ItemTokenSize = [60 25];    % 90 leaves a visible blank gap at two entries
    drawnow;
    lg.Units = 'normalized';
    lg.Position(1) = 0.1895 + (0.7155 - lg.Position(3))/2;
    lg.Position(2) = 0.9750 - lg.Position(4);
end

function xttl_(ax, YR, FSLAB, Hin, nm) %#ok<INUSD>
% The axis title, placed from the axes' OWN measured extent instead of a hand-tuned
% offset.  Two things went wrong with the obvious approaches:
%   xlabel()      MATLAB positions it from the tick-label extent, which in a stacked
%                 layout pushes it into the NEXT panel's rectangle, where that
%                 panel's opaque background hides it -- the label vanished entirely.
%   fixed offset  every change of tick font size or panel height re-broke it, and it
%                 kept landing on the tick numbers.
% TightInset(2) is the margin the axes needs below itself for its tick labels, so
% placing the box just under that clears them whatever the font size is.
    fig = ancestor(ax, 'figure');
    drawnow;                                   % TightInset is stale until drawn
    ti = ax.TightInset;   px = ax.Position;    hgt = 0.026;
    annotation(fig, 'textbox', [px(1), px(2)-ti(2)-hgt, px(3), hgt], ...
               'String', ['$\mathbf{' nm '\;(\mu m)}$'], 'Interpreter','latex', ...
               'FontSize', FSLAB, 'EdgeColor','none', 'FitBoxToText','off', ...
               'HorizontalAlignment','center', 'VerticalAlignment','middle');
end

function expo_(ax, XR, YR, FS)
% The shared factor, above the frame and RIGHT-aligned at its left edge -- which is
% where MATLAB right-aligns the y tick numbers, so the two line up.
    text(ax, XR(1), YR(2) + 0.06*diff(YR), '\times10^{-3}', ...
         'HorizontalAlignment','right','VerticalAlignment','bottom', ...
         'FontSize',FS,'FontWeight','bold','Clipping','off');
end

function xends_(ax, XR, YR, FS)
    for xv = XR
        text(ax, xv, YR(1)-0.022*diff(YR), sprintf('%g', xv), ...
             'HorizontalAlignment','center','VerticalAlignment','top', ...
             'FontSize',FS,'FontWeight','bold','Clipping','off');
    end
end

function print_(fig, out, CANV)
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
    print(fig, out, '-dpng', '-r200');   close(fig);
    fprintf(['wrote %s' newline], out);
end

function [lim, tk] = yaxis3_(a, bmax)
    nice = [1 1.25 1.5 2 2.5 3 4 5 6 7.5 8];
    span = bmax - a;   ctr = (a + bmax)/2;   need = span/4;
    for kk = floor(log10(need)) : floor(log10(need))+2
        ok = nice(nice*10^kk >= need*(1-1e-12)) * 10^kk;
        if isempty(ok), continue, end
        for st = ok
            t2  = round(ctr/st)*st;
            tk  = t2 + [-st 0 st];
            lim = [tk(1)-st, tk(end)+st];
            if lim(1) <= a + 1e-12 && lim(2) >= bmax - 1e-12, return, end
        end
    end
    error('yaxis3_:none','no step covers [%g %g]', a, bmax);
end
