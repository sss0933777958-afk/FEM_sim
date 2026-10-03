%% dif2.m -- what widening the data range 150 -> 200 um does to the gradient, asked
%%            of each basis separately: one difference figure for the RBF, one for the
%%            spherical harmonics.
%
%  Both bases are set by the SAME rule -- among the settings that pass the smoothness
%  gate, the one with the lowest working-region NMAE:
%      RBF   (rho, lam) from rbf_v150.mat / rbf_v200.mat   (the gate + min-NMAE pick
%                                                            those files already store)
%      sph   L          from s_minnmae.mat                  (BEST, same rule)
%  so the two panels answer the same question and can be read side by side.
%
%  Each model is fitted on every node inside its own radius and evaluated on the same
%  1001 equidistant points per actuator axis over |s| <= 150 um.  z_a is negated (the
%  excited P6 lies at negative z_a), as in every other figure here.
%
%  Output: data/long2016_hexapole_halfcut/.mat/dif2.mat
%          figures/long2016_hexapole_halfcut/dif_rbf.png
%          figures/long2016_hexapole_halfcut/dif_sph.png      [ADDED 2026-09-18]

clearvars;  clc;
HERE = fileparts(mfilename('fullpath'));  FMX = fileparts(fileparts(HERE));
MAIN = fileparts(fileparts(fileparts(FMX)));
FLX  = fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT = fullfile(FMX,'utils','data');
FIG = fullfile(FMX,'figures','long2016_hexapole_halfcut','current');

% [ADDED 2026-09-18 user's call] FIXPAR pins the SAME setting at both ranges, so the
% only thing that changes between the two curves is the data range.  With FIXPAR the
% gate's own per-range pick is ignored (rho went 490 -> 340 and L went 9 -> 8, so the
% earlier comparison mixed two variables); outputs get a _fix suffix and do not
% overwrite the per-range-pick versions.
FIXPAR   = true;
FIX_RHO  = 490;     FIX_LAM = 1e-8;     FIX_L = 9;
SUF      = '';   if FIXPAR, SUF = '_fix';  end

RL = [150 200];   SPAN = 150;   NQ = 1001;
POLE = [1 3 6];   DSGN = [1 1 -1];   ANM = {'x_a','y_a','z_a'};
s = linspace(-SPAN,SPAN,NQ).';

%% ---- the settings each basis was given --------------------------------------
RHO = zeros(1,2);  LAM = zeros(1,2);  LDEG = zeros(1,2);
if FIXPAR
    RHO(:) = FIX_RHO;   LAM(:) = FIX_LAM;   LDEG(:) = FIX_L;
    fprintf(['FIXED setting at both ranges: RBF rho = %g um, lam = %.0e;  sph L = %d' ...
             newline], FIX_RHO, FIX_LAM, FIX_L);
else
    for k = 1:2
        S = load(fullfile(DAT, sprintf('rbf_v%d.mat', RL(k))));
        [~,iw] = min(S.TAB.nmae);   RHO(k) = S.TAB.rho(iw);   LAM(k) = S.TAB.lam(iw);
        fprintf(['RBF @ R = %3d : rho = %g um, lam = %.0e, NMAE = %.4f %%' newline], ...
                RL(k), RHO(k), LAM(k), S.TAB.nmae(iw));
    end
    SM = load(fullfile(DAT,'s_minnmae.mat'));
    for k = 1:2
        LDEG(k) = SM.BEST(k);
        T = SM.TAB{k};   row = T(T(:,1)==LDEG(k), :);
        fprintf(['sph @ R = %3d : L = %d (K = %d), NMAE = %.4f %%' newline], ...
                RL(k), row(1), row(2), row(4));
    end
end

%% ---- the gradients ----------------------------------------------------------
GR = cell(1,2);   GS = cell(1,2);
for k = 1:2
    R0 = RL(k);
    [~, gR, iN] = rbf_field(R0, struct('rho',RHO(k), 'lam',LAM(k)));
    [~, gS]     = sph_field(R0, struct('P',iN.P,'B',iN.B6), ...
                            struct('L',LDEG(k), 'Rn',R0, 'quiet',true, 'condmax',1500));
    A = zeros(NQ,3);   B = zeros(NQ,3);
    for a = 1:3
        Pq = zeros(NQ,3);   Pq(:,a) = s;
        I  = zeros(6,1);    I(POLE(a)) = 1;
        Q = gR(Pq, I);        A(:,a) = DSGN(a)*Q(:,a);
        Q = gS(Pq, I, 'du');  B(:,a) = DSGN(a)*Q(:,a);
    end
    GR{k} = A;   GS{k} = B;
end
DR = GR{2} - GR{1};        % RBF : 200 minus 150
DS = GS{2} - GS{1};        % sph : 200 minus 150

%% ---- mean_error --------------------------------------------------------------
% [MODIFIED 2026-09-18 user's call] the normalisation the user wrote out:
%
%     sum_i |grad_A(p_i) - grad_B(p_i)| / sum_i grad_B(p_i) * 100 %
%
% i.e. the SUBTRACTED model's own gradient is the reference, summed over the same
% 1001 points.  (Point counts match, so this equals mean|D| / mean(grad_B).)  The
% gradients here are the displayed ones, positive on all three axes after the z_a
% negation, so no absolute value is needed in the denominator.
CMP = { DR,          GR{1}, 'RBF,  R = 200 minus R = 150   / grad RBF150'
        DS,          GS{1}, 'sph,  R = 200 minus R = 150   / grad sph150'
        GR{1}-GS{1}, GS{1}, 'RBF minus sph, both at R = 150 / grad sph150'
        GR{2}-GS{2}, GS{2}, 'RBF minus sph, both at R = 200 / grad sph200' };
ME = zeros(1,4);   PCT = zeros(4,3);   PCTA = zeros(1,4);
for q = 1:size(CMP,1)
    D = CMP{q,1};   G0 = CMP{q,2};
    fprintf([newline '%s' newline], CMP{q,3});
    fprintf(['%6s %16s %16s %16s %16s' newline], ...
            'axis','mean(diff)','mean|diff|','max|diff|','sum|D|/sum(ref)');
    for a = 1:3
        PCT(q,a) = sum(abs(D(:,a))) / sum(G0(:,a)) * 100;
        fprintf('%6s %16.3e %16.3e %16.3e %15.3f %%\n', ANM{a}, mean(D(:,a)), ...
                mean(abs(D(:,a))), max(abs(D(:,a))), PCT(q,a));
    end
    ME(q)   = mean(abs(D(:)));
    PCTA(q) = sum(abs(D(:))) / sum(G0(:)) * 100;
    fprintf(['mean_error = %.4e mT^2/um   (%.3f %%)' newline], ME(q), PCTA(q));
end
save(fullfile(DAT,['dif2' SUF '.mat']), 's','GR','GS','DR','DS','CMP','ME','PCT','PCTA', ...
     'RHO','LAM','LDEG','RL','SPAN','NQ','POLE','DSGN','ANM');

%% ---- figures -------------------------------------------------------------------
FS = 44;  FSLAB = 26;  FSLEG = 45;  CANV = 14.5;  LWBOX = 5.0;  LW = 7;
XR = [-SPAN SPAN];   XT = [-SPAN/2 0 SPAN/2];
C150 = [0.85 0.10 0.10];   C200 = [0.05 0.10 0.95];
NM = {'R = 150 \mum', 'R = 200 \mum'};

% ---- (a) the gradients themselves, the two ranges overlaid --------------------
% One figure per basis.  These carry a legend, so the panels are a little shorter
% than the difference figures below (which need no legend).
HG = 0.1617;   GAPG = 0.1200;   Y0G = 0.1550;   YLG = 0.9750;
for q = 1:2
    if q == 1, GG = GR;  out = ['grad_rbf' SUF '.png'];  else, GG = GS;  out = ['grad_sph' SUF '.png'];  end
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    for a = 1:3
        yb = Y0G + (3-a)*(HG+GAPG);
        ax = axes('Parent',fig,'Position',[0.1895 yb 0.7155 HG]);  hold(ax,'on');
        try, ax.Toolbar = []; end
        try, ax.Interactions = []; end
        h1 = plot(ax, s, GG{1}(:,a), '-',  'Color', C150, 'LineWidth', LW);
        h2 = plot(ax, s, GG{2}(:,a), '--', 'Color', C200, 'LineWidth', LW-3);
        style_(ax, FS, LWBOX);
        v = [GG{1}(:,a); GG{2}(:,a)];
        [YRa, YTa] = yaxis3_(min(v), max(v));
        xlim(ax, XR);  set(ax,'XTick',XT);   ylim(ax, YRa);  set(ax,'YTick',YTa);
        ylabel(ax, {['$\mathbf{d(b\cdot b)/d' ANM{a} '}$'], '$\mathbf{(mT^{2}/\mu m)}$'}, ...
               'Interpreter','latex','FontSize',FSLAB);
        if a == 1
            lg = legend(ax, [h1 h2], NM, 'Interpreter','tex', ...
                        'Location','northwest', 'NumColumns',2);
            lg.FontSize = FSLEG;   lg.FontWeight = 'bold';
            lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.Color = 'w';
            lg.ItemTokenSize = [60 25];
            drawnow;   lg.Units = 'normalized';
            lg.Position(1) = 0.1895 + (0.7155 - lg.Position(3))/2;
            lg.Position(2) = YLG - lg.Position(4);
            ax.Position = [0.1895 yb 0.7155 HG];   % undo any resize the legend caused
        end
        xends_(ax, XR, YRa, FS);
        xttl_(ax, FSLAB, ANM{a});
        hold(ax,'off');
    end
    print_(fig, fullfile(FIG,out), CANV);
end

% ---- (b) the differences ------------------------------------------------------
H = 0.1817;   GAP = 0.1200;   Y0 = 0.1550;
for q = 1:2
    if q == 1, D = DR;  out = ['dif_rbf' SUF '.png'];  else, D = DS;  out = ['dif_sph' SUF '.png'];  end
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    for a = 1:3
        yb = Y0 + (3-a)*(H+GAP);
        ax = axes('Parent',fig,'Position',[0.1895 yb 0.7155 H]);  hold(ax,'on');
        try, ax.Toolbar = []; end
        try, ax.Interactions = []; end
        Dp = D(:,a) * 1e3;
        plot(ax, s, Dp, '-', 'Color', [0 0 0], 'LineWidth', LW);
        style_(ax, FS, LWBOX);
        [YRa, YTa] = yaxis3_(min(Dp), max(Dp));
        xlim(ax, XR);  set(ax,'XTick',XT);   ylim(ax, YRa);  set(ax,'YTick',YTa);
        ax.YAxis.Exponent = 0;
        if a == 1, expo_(ax, XR, YRa, FS); end
        ylabel(ax, {['$\mathbf{\Delta\,d(b\cdot b)/d' ANM{a} '}$'], '$\mathbf{(mT^{2}/\mu m)}$'}, ...
               'Interpreter','latex','FontSize',FSLAB);
        xends_(ax, XR, YRa, FS);
        xttl_(ax, FSLAB, ANM{a});
        hold(ax,'off');
    end
    print_(fig, fullfile(FIG,out), CANV);
end

%% ---- local helpers -----------------------------------------------------------
function style_(ax, FS, LWBOX)
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
           'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
           'XTickLabelRotation',0,'YTickLabelRotation',0);
end

function xttl_(ax, FSLAB, nm)
% Placed from the axes' own TightInset, not with xlabel(): in a stacked layout MATLAB
% pushes an xlabel into the next panel's rectangle, where its opaque background hides
% it.  TightInset(2) is the room the tick numbers need, so this clears them at any
% font size.
    fig = ancestor(ax,'figure');   drawnow;
    ti = ax.TightInset;   px = ax.Position;   hgt = 0.026;
    annotation(fig,'textbox',[px(1), px(2)-ti(2)-hgt, px(3), hgt], ...
               'String',['$\mathbf{' nm '\;(\mu m)}$'],'Interpreter','latex', ...
               'FontSize',FSLAB,'EdgeColor','none','FitBoxToText','off', ...
               'HorizontalAlignment','center','VerticalAlignment','middle');
end

function expo_(ax, XR, YR, FS)
% One shared factor, above the top frame and right-aligned at its left edge, which is
% where MATLAB right-aligns the y tick numbers.
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
            t2 = round(ctr/st)*st;   tk = t2 + [-st 0 st];
            lim = [tk(1)-st, tk(end)+st];
            if lim(1) <= a + 1e-12 && lim(2) >= bmax - 1e-12, return, end
        end
    end
    error('yaxis3_:none','no step covers [%g %g]', a, bmax);
end
