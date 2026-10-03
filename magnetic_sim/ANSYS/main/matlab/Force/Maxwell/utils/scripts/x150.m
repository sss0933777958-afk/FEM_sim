%% x150.m -- both field models fitted on the R <= 150 um nodes only, each taken at
%%            the setting that passes the smoothness gate, then drawn side by side.
%
%  THE GATE (user's wording, 2026-09-17): put NQ equidistant points on each of the
%  three actuator axes, evaluate the model's gradient there, and ask whether the
%  curve those points trace is smooth.  Operationally that is
%
%       d3(b.b)/ds^3 keeps ONE SIGN over |s| <= 150 um, on all three axes
%
%  because d3(b.b)/ds^3 IS the second derivative of the gradient curve: no sign
%  change means the gradient is convex (or concave) throughout and cannot carry a
%  ripple.  The derivative is analytic, so the verdict does not depend on how the
%  points are spaced -- the sampling only decides how fine a ripple can hide
%  between two points.  NQ = 1001 per axis (0.3 um spacing).
%
%      x_a <- P1 excited      y_a <- P3      z_a <- P6      each at 1 A
%
%  BOTH models see exactly the same 1771 nodes: rbf_field selects them and hands
%  info.P / info.B6 straight to sph_field, so nothing but the basis differs.
%
%  SIGNS.  z_a is the P5 direction, so the excited P6 sits at NEGATIVE z_a and its
%  gradient is negative there.  The gradient is a vector component, so the z_a panel
%  is drawn NEGATED and all three panels then read "rises towards the excited pole".
%  b.b is a scalar and is never negated.
%
%  Inputs   rbf_v150.mat  (rbf_strat.m, RUN_R=150 RUN_NCAP=Inf RUN_NQ=1001)
%           sph_v150.mat  (sph_strat.m, same, RUN_LMAX=Inf)
%  Outputs  data/long2016_hexapole_halfcut/.mat/x150.mat
%           figures/long2016_hexapole_halfcut/x150_bb.png    b.b on x_a, both models
%           figures/long2016_hexapole_halfcut/x150_g3.png    gradient, 3 axes, both
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

MODEL = 'long2016_hexapole_halfcut';   GEOM = 'tip40um';
CALNAME = 'calib_current_maxwell_R150_eighteen.mat';   % the flux calibration to compare against
R0   = 150;
SPAN = 150;                 % the segment drawn and judged, |s| <= SPAN
NQ   = 1001;                % equidistant points per axis (user's call)
POLE = [1 3 6];             % x_a <- P1, y_a <- P3, z_a <- P6
DSGN = [1 1 -1];            % display sign of the gradient (see SIGNS above)
ANM  = {'x_a','y_a','z_a'};

%% ---- 1. the settings that pass the gate --------------------------------------
SR  = load(fullfile(DAT,'rbf_v150.mat'));
TAB = SR.TAB;
[~, iw] = min(TAB.nmae);                       % NaN wherever the gate failed
RHO = TAB.rho(iw);   LAM = TAB.lam(iw);
fprintf(['RBF : %d of %d (rho,lam) pairs pass;  smallest passing rho = %g um' newline], ...
        nnz(TAB.pass), numel(TAB.rho), min(TAB.rho(TAB.pass)));
fprintf(['      chosen rho = %g um, lam = %.0e, NMAE = %.4f %%' newline], ...
        RHO, LAM, TAB.nmae(iw));

SS  = load(fullfile(DAT,'sph_v150.mat'));
r   = SS.RES([SS.RES.R] == R0);
LDEG = r.L;
fprintf(['sph : %d of %d degrees pass;  chosen L = %d (K = %d), NMAE = %.4f %%' newline], ...
        r.npass, r.Ltried, LDEG, r.K, r.nmae);

%% ---- 2. fit both on the same 1771 nodes --------------------------------------
[bR, gR, iR] = rbf_field(R0, struct('rho',RHO, 'lam',LAM));
assert(iR.Np == 1771, 'x150:nodes', 'expected 1771 nodes, got %d', iR.Np);
[bS, gS, iS] = sph_field(R0, struct('P',iR.P, 'B',iR.B6), ...
                         struct('L',LDEG, 'Rn',R0, 'quiet',true));
assert(iS.Np == iR.Np, 'x150:mismatch', 'the two models were given different nodes');
fprintf(['both models fitted on the same %d nodes' newline], iR.Np);

%% ---- 3. the NQ points on each axis, and the gate re-checked there ------------
% The eighteen-parameter charge model from the flux calibration, for comparison.
% Only its GRADIENT is available here: build_L is the gradient kernel, and this
% package has no evaluator for the charge model's b.b itself, so that curve is
% absent from the b.b figure rather than faked.
cfg = model_config(MODEL, GEOM);
C   = load(fullfile(FLX,'data',MODEL,'.mat',CALNAME));
assert(C.USE_BIAS == 1, 'x150:notEighteen', 'the calibration loaded is not the 18-parameter one');
lum = C.l_hat*1e6;
fprintf(['18-par: l_hat = %.2f um, gI_hat = %.4f mT/A  (%s)' newline], lum, C.gI_hat, CALNAME);

s  = linspace(-SPAN, SPAN, NQ).';
GR = zeros(NQ,3);   GS = zeros(NQ,3);   G18 = zeros(NQ,3);
nscR = zeros(1,3);  nscS = zeros(1,3);
for a = 1:3
    Pq = zeros(NQ,3);   Pq(:,a) = s;
    I  = zeros(6,1);    I(POLE(a)) = 1;
    G  = gR(Pq, I);         GR(:,a) = DSGN(a) * G(:,a);
    G  = gS(Pq, I, 'du');   GS(:,a) = DSGN(a) * G(:,a);
    Lt = build_L(Pq, lum, C.e, cfg.Pc_base);
    Gc = base_model(Lt, C.KI_bar, C.Fmap(:,POLE(a)), C.gI_hat^2 / lum);
    G18(:,a) = DSGN(a) * squeeze(Gc(a,:,1));
    d3R = gR(Pq, I, 'd3');   nscR(a) = nsc_(d3R(:,a));
    d3S = gS(Pq, I, 'd3');   nscS(a) = nsc_(d3S(:,a));
    if a == 1
        bbR = sum(bR(Pq, I).^2, 2);            % b.b on x_a, P1 excited   [mT^2]
        bbS = sum(bS(Pq, I).^2, 2);
    end
    fprintf(['%s: RBF %7.4f..%7.4f  sph %7.4f..%7.4f  18par %7.4f..%7.4f mT^2/um' ...
             '   ripples R/S = %d/%d' newline], ANM{a}, min(GR(:,a)), max(GR(:,a)), ...
            min(GS(:,a)), max(GS(:,a)), min(G18(:,a)), max(G18(:,a)), nscR(a), nscS(a));
    fprintf(['      18-par vs RBF at the far end: %+.2f %%' newline], ...
            (G18(end,a)-GR(end,a))/GR(end,a)*100);
end
fprintf(['b.b on x_a: RBF %.4f..%.4f   sph %.4f..%.4f mT^2' newline], ...
        min(bbR), max(bbR), min(bbS), max(bbS));
assert(all(nscR == 0) && all(nscS == 0), 'x150:gate', ...
       'a chosen setting does not pass the gate at NQ = %d', NQ);

save(fullfile(DAT,'x150.mat'), 's','bbR','bbS','GR','GS','G18','nscR','nscS', ...
     'RHO','LAM','LDEG','R0','SPAN','NQ','POLE','DSGN','ANM','CALNAME','lum');

%% ---- 4. figures --------------------------------------------------------------
FS = 60;  FSLAB = 36;  FSLEG = 45;  CANV = 14.5;  LWBOX = 5.0;  LW = 7;
CR = [0.85 0.10 0.10];   CS = [0.05 0.10 0.95];   C18 = [0 0 0];
NM18 = 'Eighteen parameters';
% The legend names the METHOD only.  rho / lam / L belong in the printed output and
% in the .mat (user's spec, 2026-09-17) -- spelled out here they made the legend box
% wide enough to cover the curves it was labelling.
NMR = 'RBF';
NMS = 'Spherical harmonics';
XR  = [-SPAN SPAN];   XT = [-SPAN/2 0 SPAN/2];

% ---- 4a. b.b on x_a ----------------------------------------------------------
[YR, YT] = yaxis3_(min([bbR;bbS]), max([bbR;bbS]));
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.1895 0.1550 0.7155 0.7850]);  hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
h1 = plot(ax, s, bbR, '-',  'Color', CR, 'LineWidth', LW);
h2 = plot(ax, s, bbS, '--', 'Color', CS, 'LineWidth', LW-2);
style_(ax, FS, LWBOX);
xlim(ax, XR);  set(ax,'XTick',XT);   ylim(ax, YR);  set(ax,'YTick',YT);
xends_(ax, XR, YR, FS);
xlabel(ax, '$\mathbf{x_a\;(\mu m)}$',       'Interpreter','latex','FontSize',FSLAB);
ylabel(ax, '$\mathbf{b\cdot b\;(mT^{2})}$', 'Interpreter','latex','FontSize',FSLAB);
leg_(ax, [h1 h2], {NMR, NMS}, FSLEG, LWBOX, 'northwest');
hold(ax,'off');
print_(fig, fullfile(FIG,'x150_bb.png'), CANV);

% ---- 4b. the three gradients, stacked top to bottom --------------------------
% The legend sits ABOVE all three panels (user's call), so the panels give up the
% top 0.19 of the canvas for it.  ONE COLUMN, three rows: three entries side by side
% would need about 17 in at 45 pt on a 14.5 in canvas, and two columns leaves a
% ragged empty cell -- a single left-aligned column is the tidy arrangement that
% fits.  FSLAB drops to 32 because each panel is now 2.76 in tall and the two-line
% y label would otherwise run past its own frame.
% FS drops to 44 here (the single-panel b.b figure keeps the standard 60): the tick
% numbers are sized for a full-height frame, and three stacked panels are each only
% 2.76 in tall, so 60 pt reads as oversized against them.
% [MODIFIED 2026-09-17 user's call] EVERY panel carries the x tick numbers and the x
% axis title, not only the bottom one -- this deliberately overrides rule 8.  The gap
% has to grow to 0.077 to hold both, which costs panel height, so FSLAB drops to 26
% to keep the two-line y label inside its own 2.3 in frame.
% GAP is 0.12, not the 0.10 that looked sufficient on paper.  TickLength is measured
% against the LONGER axis, so on a panel 10 in wide and 2.4 in tall the outward ticks
% on the frame stick out about 0.16 in -- the ticks along the TOP edge of the panel
% below reach up into the gap and were being struck by the title sitting above them.
H = 0.1617;   GAP = 0.1200;   Y0 = 0.1550;   YLG = 0.9750;   FSLAB = 26;   FS = 44;
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
for a = 1:3
    yb = Y0 + (3-a)*(H+GAP);                        % a = 1 on top
    ax = axes('Parent',fig,'Position',[0.1895 yb 0.7155 H]);  hold(ax,'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    % [MODIFIED 2026-09-17 user's call] the eighteen-parameter curve is dropped from
    % this figure; it is still computed and stored, just not drawn.
    h1 = plot(ax, s, GR(:,a),  '-',  'Color', CR,  'LineWidth', LW);
    h2 = plot(ax, s, GS(:,a),  '--', 'Color', CS,  'LineWidth', LW-2);
    style_(ax, FS, LWBOX);
    vall = [GR(:,a); GS(:,a)];
    [YRa, YTa] = yaxis3_(min(vall), max(vall));
    xlim(ax, XR);  set(ax,'XTick',XT);   ylim(ax, YRa);  set(ax,'YTick',YTa);
    % Two lines rather than one: at 36 pt a one-line "d(b.b)/dx_a (mT^2/um)" is
    % 5.5 in long and each of the three panels is only 3.6 in tall.
    ylabel(ax, {['$\mathbf{d(b\cdot b)/d' ANM{a} '}$'], '$\mathbf{(mT^{2}/\mu m)}$'}, ...
           'Interpreter','latex','FontSize',FSLAB);
    if a == 1
        % Built at a normal location and then MOVED, never with 'northoutside':
        % the outside locations resize the axes they belong to, and the three
        % panels would stop matching (figure-style.md / 2026-09-02 trap).
        lg = leg_(ax, [h1 h2], {NMR, NMS}, FSLEG, LWBOX, 'northwest');
        drawnow;                                   % Position is stale until drawn
        lg.Units = 'normalized';
        lg.Position(1) = 0.1895 + (0.7155 - lg.Position(3))/2;   % centred on the frame
        lg.Position(2) = YLG - lg.Position(4);                   % flush to the top
        ax.Position = [0.1895 yb 0.7155 H];        % undo any resize the legend caused
    end
    xends_(ax, XR, YRa, FS);
    xttl_(ax, YRa, FSLAB, H*CANV, ANM{a});
    hold(ax,'off');
end
print_(fig, fullfile(FIG,'x150_g3.png'), CANV);

% ---- 4c. the difference, RBF - spherical harmonics ---------------------------
% One series per panel, so no legend; the y label says what it is.
% [MODIFIED 2026-09-18 user's call] the plotted quantity is now a PERCENTAGE, the
% pointwise relative difference against the harmonics' own gradient:
%
%     100 * ( grad_RBF(p) - grad_SPH(p) ) / grad_SPH(p)
%
% The sign is kept, so the shape reads the same as the old absolute version.  The
% harmonics are the denominator only because something has to be: neither model is
% truth, so this is a relative DIFFERENCE, not RBF's error.  grad_SPH runs 0.19 to
% 1.10 mT^2/um and never approaches zero, so the pointwise division is safe.
%   Note the per-axis summary printed below is sum|D| / sum(grad_SPH), a ratio of
%   totals -- it is NOT the mean of this curve, and the two differ because the
%   gradient varies by a factor 5.8 along each axis.
D = GR - GS;
H = 0.1817;   GAP = 0.1200;          % no legend here, so the panels get that space
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
for a = 1:3
    yb = Y0 + (3-a)*(H+GAP);
    ax = axes('Parent',fig,'Position',[0.1895 yb 0.7155 H]);  hold(ax,'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    Dp = 100 * D(:,a) ./ GS(:,a);              % per cent, pointwise
    plot(ax, s, Dp, '-', 'Color', [0 0 0], 'LineWidth', LW);
    style_(ax, FS, LWBOX);
    [YRa, YTa] = yaxis3_(min(Dp), max(Dp));
    xlim(ax, XR);  set(ax,'XTick',XT);   ylim(ax, YRa);  set(ax,'YTick',YTa);
    ax.YAxis.Exponent = 0;                     % a percentage needs no shared factor
    ylabel(ax, {['$\mathbf{(\nabla_{RBF}-\nabla_{SPH})/\nabla_{SPH}}$'], ...
                ['$\mathbf{on\;' ANM{a} '\;(\%)}$']}, ...
           'Interpreter','latex','FontSize',FSLAB);
    xends_(ax, XR, YRa, FS);
    xttl_(ax, YRa, FSLAB, H*CANV, ANM{a});
    hold(ax,'off');
end
print_(fig, fullfile(FIG,'x150_dif.png'), CANV);

%% ---- the summary, as a percentage ----------------------------------------------
%     PCT(a) = sum_i |grad_RBF - grad_SPH| / sum_i grad_SPH * 100
% a ratio of totals, the same form as the range-comparison table, with the harmonics
% as the reference.  PCTPT is the mean of the plotted pointwise curve, kept alongside
% because the two are different statistics and are easy to confuse.
fprintf([newline '%6s %14s %14s %12s %12s' newline], ...
        'axis','mean|RBF-sph|','max|RBF-sph|','sum/sum %','mean ptwise %');
PCT = zeros(1,3);   PCTPT = zeros(1,3);
for a = 1:3
    PCT(a)   = sum(abs(D(:,a))) / sum(GS(:,a)) * 100;
    PCTPT(a) = mean(abs(D(:,a) ./ GS(:,a))) * 100;
    fprintf('%6s %14.3e %14.3e %11.3f %% %11.3f %%\n', ANM{a}, mean(abs(D(:,a))), ...
            max(abs(D(:,a))), PCT(a), PCTPT(a));
end
mean_error = mean(abs(D(:)));
PCTA = sum(abs(D(:))) / sum(GS(:)) * 100;
fprintf([newline 'all three axes : mean|RBF-sph| = %.4e mT^2/um,  sum/sum = %.3f %%' newline], ...
        mean_error, PCTA);
save(fullfile(DAT,'x150_dif.mat'), 's','D','GR','GS','mean_error','PCT','PCTPT','PCTA', ...
     'ANM','RHO','LAM','LDEG');

%% ---- local helpers -----------------------------------------------------------
function style_(ax, FS, LWBOX)
% Ticks stay on the default (tex) interpreter: latex ignores FontWeight, so a latex
% tick label can never be bold (figure-style.md, trap 3).  Only axis titles are latex.
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
           'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
           'XTickLabelRotation',0,'YTickLabelRotation',0);
end

function lg = leg_(ax, h, names, FSLEG, LWBOX, loc)
    lg = legend(ax, h, names, 'Interpreter','tex', 'Location', loc, 'NumColumns',2);
    lg.FontSize = FSLEG;   lg.FontWeight = 'bold';
    lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.Color = 'w';
    % 90 pt, not the usual 55: MATLAB scales a dash's length with the line's width,
    % so a thick dashed sample fits only one gap-plus-dash into 55 pt and reads as a
    % short line indented from the solid one above it, rather than as a dashed line.
    lg.ItemTokenSize = [60 25];    % 90 leaves a visible blank gap at two entries       % the default 30 pt does not grow with FontSize
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
% The shared factor, written above the frame and RIGHT-aligned at the frame's left
% edge -- which is exactly where MATLAB right-aligns the y tick numbers, so the two
% line up.
    text(ax, XR(1), YR(2) + 0.06*diff(YR), '\times10^{-3}', ...
         'HorizontalAlignment','right','VerticalAlignment','bottom', ...
         'FontSize',FS,'FontWeight','bold','Clipping','off');
end

function xends_(ax, XR, YR, FS)
% rule 4: the x ends carry numbers but no tick mark, so they are placed by hand
    for xv = XR
        text(ax, xv, YR(1)-0.022*diff(YR), sprintf('%g', xv), ...
             'HorizontalAlignment','center','VerticalAlignment','top', ...
             'FontSize',FS,'FontWeight','bold','Clipping','off');
    end
end

function print_(fig, out, CANV)
% print + explicit PaperPosition: exportgraphics crops the margins and a square
% canvas then comes out not square (figure-style.md, trap 2)
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
    print(fig, out, '-dpng', '-r200');   close(fig);
    fprintf(['wrote %s' newline], out);
end

function n = nsc_(v)
%NSC_  number of sign changes in v, ignoring exact zeros
    g = sign(v);  nz = find(g~=0);  n = sum(diff(g(nz)) ~= 0);
end

function [lim, tk] = yaxis3_(a, bmax)
% Three equally spaced ticks with the SAME gap to each end of the frame, so the
% frame is exactly 4 steps tall and neither end lands on a tick (rules 4 and 5).
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
