%% rbf_x150.m -- R = 150 um, all nodes: which (rho,lam) pass the smoothness gate,
%%                then b.b on x_a, and the gradient on all three actuator axes.
%
%  The gate is the one used throughout this study: d3(b.b)/ds^3 keeps ONE SIGN on
%  all three actuator axes (x_a <- P1, y_a <- P3, z_a <- P6) over |s| <= 150 um.
%  It is a property of the MODEL, not of any single axis -- a pair that ripples on
%  y_a is rejected even if x_a looks clean.
%
%  Reads the pass/fail table written by temp_code/rbf_strat.m (RUN_R = 150,
%  RUN_NCAP = Inf, default 1188-pair grid) and refits the winning pair to draw the
%  curves.  Computing and drawing sit together here because this is temp_code
%  (plot-scripts-pure, "何時不適用").
%
%  SIGNS (user's spec, 2026-09-17).  z_a is the P5 direction, so the excited P6 sits
%  at NEGATIVE z_a and its gradient is negative there.  The gradient is a vector
%  component, so the z_a panel is drawn NEGATED -- that makes all three panels read
%  the same way, "rises towards the excited pole".  b.b is a scalar and is NOT
%  negated anywhere: flipping it would hide which side the pole is on.
%
%  Output: data/long2016_hexapole_halfcut/.mat/rbf_x150.mat
%          figures/long2016_hexapole_halfcut/rbf_x150_bb.png    b.b on x_a
%          figures/long2016_hexapole_halfcut/rbf_x150_g3.png    gradient, 3 axes
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

R0   = 150;
SPAN = 150;                 % the segment drawn, |s| <= SPAN
NS   = 1201;
POLE = [1 3 6];             % x_a <- P1, y_a <- P3, z_a <- P6
DSGN = [1 1 -1];            % display sign of the gradient (see SIGNS above)
ANM  = {'x_a','y_a','z_a'};

%% ---- 1. the pass map ---------------------------------------------------------
S   = load(fullfile(DAT,'rbf_v150.mat'));      % written by rbf_strat.m, RUN_R = 150
TAB = S.TAB;   RHOL = S.RHOL;   LAML = S.LAML;
np  = numel(TAB.rho);
fprintf(['grid: %d rho x %d lam = %d pairs, %d pass the gate (%.1f %%)' newline], ...
        numel(RHOL), numel(LAML), np, nnz(TAB.pass), nnz(TAB.pass)/np*100);
rmin = min(TAB.rho(TAB.pass));
lnev = LAML(~ismember(LAML, unique(TAB.lam(TAB.pass))));
fprintf(['smallest rho with any passing lam: %g um' newline], rmin);
fprintf(['lam that never pass at any rho: %s' newline], strjoin(compose('%g', lnev), ' '));

[~, iw] = min(TAB.nmae);                       % NaN where the gate failed
RHO = TAB.rho(iw);   LAM = TAB.lam(iw);
fprintf(['winner (lowest NMAE among the %d that pass): rho = %g um, lam = %.0e, ' ...
         'NMAE = %.4f %%' newline], nnz(TAB.pass), RHO, LAM, TAB.nmae(iw));

%% ---- 2. refit the winner and evaluate on the three axes ----------------------
[b, grad, info] = rbf_field(R0, struct('rho',RHO, 'lam',LAM));
assert(info.Np == 1771, 'rbf_x150:nodes', ...
       'expected the 1771 nodes inside R = 150 um, got %d', info.Np);

s  = linspace(-SPAN, SPAN, NS).';
GX = zeros(NS,3);   D3 = zeros(NS,3);   nsc = zeros(1,3);
for a = 1:3
    Pq = zeros(NS,3);   Pq(:,a) = s;
    I  = zeros(6,1);    I(POLE(a)) = 1;        % that axis's pole at 1 A
    G  = grad(Pq, I);        GX(:,a) = DSGN(a) * G(:,a);
    d3 = grad(Pq, I, 'd3');  D3(:,a) = d3(:,a);
    nsc(a) = nsc_(d3(:,a));
    if a == 1
        bb = sum(b(Pq, I).^2, 2);              % b.b on x_a, P1 excited   [mT^2]
    end
    fprintf(['%s: d(b.b)/d%s = %8.4f .. %8.4f mT^2/um   d3 sign changes %d' newline], ...
            ANM{a}, ANM{a}, min(GX(:,a)), max(GX(:,a)), nsc(a));
end
fprintf(['b.b on x_a: %.4f .. %.4f mT^2' newline], min(bb), max(bb));
assert(all(nsc == 0), 'rbf_x150:gate', 'the winning pair does not pass the gate');

save(fullfile(DAT,'rbf_x150.mat'), 's','bb','GX','D3','nsc','RHO','LAM','R0','SPAN', ...
     'POLE','DSGN','ANM','TAB','RHOL','LAML','rmin');

%% ---- 3. figures --------------------------------------------------------------
FS = 60;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;  LW = 7;
COL = [0.85 0.10 0.10];
XR  = [-SPAN SPAN];   XT = [-SPAN/2 0 SPAN/2];      % rule 5: 3 ticks, equal margins

% ---- 3a. b.b alone -----------------------------------------------------------
[YR, YT] = yaxis3_(min(bb), max(bb));
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.1895 0.1550 0.7155 0.7850]);  hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
plot(ax, s, bb, '-', 'Color', COL, 'LineWidth', LW);
style_(ax, FS, LWBOX);
xlim(ax, XR);  set(ax,'XTick',XT);   ylim(ax, YR);  set(ax,'YTick',YT);
xends_(ax, XR, YR, FS);
xlabel(ax, '$\mathbf{x_a\;(\mu m)}$',        'Interpreter','latex','FontSize',FSLAB);
ylabel(ax, '$\mathbf{b\cdot b\;(mT^{2})}$',  'Interpreter','latex','FontSize',FSLAB);
hold(ax,'off');
print_(fig, fullfile(FIG,'rbf_x150_bb.png'), CANV);

% ---- 3b. the three gradients, stacked top to bottom --------------------------
% Three panels in a square canvas leave each one 3.6 in tall, and a one-line
% "d(b.b)/dx_a (mT^2/um)" at 36 pt is 5.5 in long -- it would run off both ends of
% its own panel.  Split over two lines, which costs left margin (there is room)
% instead of panel height (there is not).
H = 0.2450;   GAP = 0.0450;   Y0 = 0.1550;
AX = gobjects(1,3);
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
for a = 1:3
    yb = Y0 + (3-a)*(H+GAP);                        % a = 1 on top
    AX(a) = axes('Parent',fig,'Position',[0.1895 yb 0.7155 H]);  hold(AX(a),'on');
    try, AX(a).Toolbar = []; end
    try, AX(a).Interactions = []; end
    plot(AX(a), s, GX(:,a), '-', 'Color', COL, 'LineWidth', LW);
    style_(AX(a), FS, LWBOX);
    [YRa, YTa] = yaxis3_(min(GX(:,a)), max(GX(:,a)));
    xlim(AX(a), XR);  set(AX(a),'XTick',XT);   ylim(AX(a), YRa);  set(AX(a),'YTick',YTa);
    ylabel(AX(a), {['$\mathbf{d(b\cdot b)/d' ANM{a} '}$'], '$\mathbf{(mT^{2}/\mu m)}$'}, ...
           'Interpreter','latex','FontSize',FSLAB);
    if a < 3
        set(AX(a),'XTickLabel',[]);   xlabel(AX(a),'');   % rule 8
    else
        xends_(AX(a), XR, YRa, FS);
        xlabel(AX(a), '$\mathbf{s\;(\mu m)}$', 'Interpreter','latex','FontSize',FSLAB);
    end
    hold(AX(a),'off');
end
print_(fig, fullfile(FIG,'rbf_x150_g3.png'), CANV);

%% ---- local helpers -----------------------------------------------------------
function style_(ax, FS, LWBOX)
% Ticks and legend stay on the default (tex) interpreter: latex ignores FontWeight,
% so a latex tick label can never be bold (figure-style.md, trap 3).  Only the axis
% titles are latex.
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
           'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
           'XTickLabelRotation',0,'YTickLabelRotation',0);
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
% print + an explicit PaperPosition, because exportgraphics crops the margins and a
% square canvas then comes out not square (figure-style.md, trap 2)
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
