% bdotb  b.b and its gradient along each actuator axis, RBF vs solid harmonics.
%
%   Fitted on all 1771 FEM nodes with |p| <= 150 um, six excitations at once, with
%   the parameters that won the R = 150 row of the range scan:
%       RBF   rho = 460 um, lam = 1e-7   (red)
%       sph   L = 7, K = 63              (blue)
%   Both models were handed the SAME nodes -- rbf_field loads them and its info.P /
%   info.B6 go straight into sph_field.
%
%       x_a  <- P1 excited      y_a  <- P3      z_a  <- P6
%
%   Each axis is drawn as it is, nothing sign-flipped: P6 lies on -z_a, so on that
%   axis b.b rises towards NEGATIVE s and its derivative is consequently negative
%   there.  That sign is the real direction the force points (towards P6).
%
%   Three sets of figures:
%       bdotb_<ax>          both models overlaid          (3, with legend)
%       bdotb_dif_<ax>      RBF minus sph, same units     (3, black: it belongs to
%                                                          neither model, so neither
%                                                          model's colour is used)
%       dbb_<ax>            d(b.b)/ds, both overlaid      (3, with legend)
%       dbb_dif_<ax>        RBF minus sph, gradient       (3, black)
%   Panels of the same set share one y frame so they can be laid side by side.
%
%   Reads only (plot-scripts-pure); bdotb.mat is written by temp_code/bdotb_fit.m.
%
% Output: figures/long2016_hexapole_halfcut/bdotb_{,dif_}{xa,ya,za}.png
%         figures/long2016_hexapole_halfcut/dbb_{,dif_}{xa,ya,za}.png

HERE = fileparts(mfilename('fullpath'));
FMX  = fileparts(fileparts(HERE));
DAT  = fullfile(FMX,'utils','data');
FIG  = fullfile(FMX,'figures','long2016_hexapole_halfcut','current');

S  = load(fullfile(DAT,'bdotb.mat'));
AX = S.AX;
fprintf(['RBF rho = %g um, lam = %.0e   |   spherical harmonics L = %d' newline], S.RHO, S.LAM, S.LSPH);

ST = struct('FS',60, 'FSLEG',45, 'FSLAB',36, 'CANV',14.5, 'LWBOX',5.0, 'LW',7);
cR = [0.85 0.10 0.10];   cS = [0.05 0.10 0.95];   cD = [0 0 0];
[XR, XT] = axis_sym_(S.R0);

% ---- 1. b.b, both models overlaid ---------------------------------------------
[YR, YT] = ylim_odd_(max(arrayfun(@(a) max([a.bbR; a.bbS]), AX)));
for a = 1:numel(AX)
    A = AX(a);
    out = fullfile(FIG, sprintf('bdotb_%s.png', erase(A.name,'_')));
    % the two agree to 0.09 %, so equal solid lines would hide one completely:
    % RBF goes down thick and solid, sph on top thin and dashed, and where they
    % coincide the dashes read as a broken line on a thick band.
    draw_(A.s, {A.bbR, A.bbS}, {cR, cS}, {'RBF','Spherical harmonics'}, ...
          XR, XT, YR, YT, '$\mathbf{b\cdot b\;(mT^{2})}$', false, out, ST, ...
          {'-','--'}, {ST.LW+4, ST.LW-2}, A.name);
    fprintf(['bdotb     %s (P%d): RBF %.1f..%.1f, sph %.1f..%.1f mT^2 -> %s' newline], ...
            A.name, A.pole, min(A.bbR), max(A.bbR), min(A.bbS), max(A.bbS), out);
end

% ---- 2. the difference, RBF - sph ---------------------------------------------
dmax = max(arrayfun(@(a) max(abs(a.bbR - a.bbS)), AX));
[YRd, YTd] = ylim_sym_(dmax);
for a = 1:numel(AX)
    A = AX(a);   d = A.bbR - A.bbS;
    out = fullfile(FIG, sprintf('bdotb_dif_%s.png', erase(A.name,'_')));
    draw_(A.s, {d}, {cD}, {}, XR, XT, YRd, YTd, ...
          '$\mathbf{\Delta(b\cdot b)\;(mT^{2})}$', true, out, ST, [], [], A.name);
    fprintf(['bdotb_dif %s (P%d): %+.4f .. %+.4f mT^2 (%.3f %% of peak) -> %s' newline], ...
            A.name, A.pole, min(d), max(d), max(abs(d))/max(A.bbR)*100, out);
end

% ---- 3. gradient, both models overlaid (z_a sign-flipped) ---------------------
%   z_a is negated here, unlike the b.b panels.  The gradient is a vector COMPONENT,
%   so flipping its sign on the axis whose pole sits at -z_a makes all three panels
%   read the same way: positive = the force pulls towards the excited pole.  (b.b is
%   a scalar and was left alone, because there the sign carries no direction and
%   flipping would only hide which end the pole is on.)  Same convention as grad3.m.
%   After the flip all three are positive throughout, so one zero-based frame serves
%   every panel and none of them wastes half the box.
SGN  = [1 1 -1];
gmax = max(arrayfun(@(a) max(abs([a.dbbR; a.dbbS])), AX));
gs   = nice_step_(gmax/4);
YRg  = [0 4*gs];   YTg = (1:3)*gs;
for a = 1:numel(AX)
    A  = AX(a);   yR = SGN(a)*A.dbbR;   yS = SGN(a)*A.dbbS;
    out = fullfile(FIG, sprintf('dbb_%s.png', erase(A.name,'_')));
    draw_(A.s, {yR, yS}, {cR, cS}, {'RBF','Spherical harmonics'}, ...
          XR, XT, YRg, YTg, ['$\mathbf{d(b\cdot b)/d' A.name '\;(mT^{2}/\mu m)}$'], false, out, ST, ...
          {'-','--'}, {ST.LW+4, ST.LW-2}, A.name);
    if SGN(a) < 0, tag = ' (negated)'; else, tag = ''; end
    fprintf(['dbb       %s (P%d)%s: RBF %.3f..%.3f, sph %.3f..%.3f mT^2/um -> %s' newline], ...
            A.name, A.pole, tag, min(yR), max(yR), min(yS), max(yS), out);
end

% ---- 4. gradient difference, RBF - sph ----------------------------------------
dgmax = max(arrayfun(@(a) max(abs(a.dbbR - a.dbbS)), AX));
[YRdg, YTdg] = ylim_sym_(dgmax);
for a = 1:numel(AX)
    A = AX(a);   d = SGN(a)*(A.dbbR - A.dbbS);
    out = fullfile(FIG, sprintf('dbb_dif_%s.png', erase(A.name,'_')));
    draw_(A.s, {d}, {cD}, {}, XR, XT, YRdg, YTdg, ...
          ['$\mathbf{\Delta[d(b\cdot b)/d' A.name ']\;(mT^{2}/\mu m)}$'], true, out, ST, [], [], A.name);
    fprintf(['dbb_dif   %s (P%d): %+.4f .. %+.4f mT^2/um (%.2f %% of peak) -> %s' newline], ...
            A.name, A.pole, min(d), max(d), max(abs(d))/max(abs(A.dbbR))*100, out);
end

% ---- local helpers ------------------------------------------------------------
function draw_(s, Y, COL, LAB, XR, XT, YR, YT, ylab, zeroline, out, ST, STY, LWS, axnm)
%   STY / LWS are optional per-series line styles and widths.
%DRAW_  one panel, figure-style: 60 pt bold ticks, 5.0 pt box, square canvas printed
%   with an explicit PaperPosition (exportgraphics would crop the margins and break
%   the square), x ends labelled by hand, y ends unlabelled, legend only if LAB given.
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 ST.CANV ST.CANV]);
    ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]);  hold(ax,'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    if zeroline, yline(ax, 0, 'k:', 'LineWidth', ST.LWBOX*0.6);  end

    if nargin < 13 || isempty(STY), STY = repmat({'-'},   1, numel(Y)); end
    if nargin < 14 || isempty(LWS), LWS = repmat({ST.LW}, 1, numel(Y)); end
    h = gobjects(1,numel(Y));
    for q = 1:numel(Y)
        h(q) = plot(ax, s, Y{q}, STY{q}, 'Color', COL{q}, 'LineWidth', LWS{q});
    end

    set(ax,'FontSize',ST.FS,'FontWeight','bold','LineWidth',ST.LWBOX, ...
           'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
           'XTickLabelRotation',0,'YTickLabelRotation',0);
    xlim(ax, XR);   set(ax,'XTick', XT);
    ylim(ax, YR);   set(ax,'YTick', YT);
    for xv = XR                                   % rule 4: label both x ends by hand
        text(ax, xv, YR(1)-0.022*diff(YR), sprintf('%g', xv), ...
             'HorizontalAlignment','center','VerticalAlignment','top', ...
             'FontSize',ST.FS,'FontWeight','bold','Clipping','off');
    end
    if nargin < 15 || isempty(axnm), axnm = 's'; end
    % NOTE concatenation, not sprintf: sprintf would eat \m, \c and \; as escapes
    xlabel(ax, ['$\mathbf{' axnm '\;(\mu m)}$'], 'Interpreter','latex','FontSize',ST.FSLAB);
    ylabel(ax, ylab,                    'Interpreter','latex','FontSize',ST.FSLAB);
    if ~isempty(LAB)
        % put the legend on whichever top corner the curves leave free: z_a peaks on
        % the left while x_a / y_a peak on the right, so a fixed corner collides on
        % one of them.  Compare the outer thirds and take the lower side.
        yl = max(cellfun(@(v) max(v(1:round(numel(v)/3))),      Y));
        yr = max(cellfun(@(v) max(v(round(2*numel(v)/3):end)),  Y));
        if yl > yr, loc = 'northeast'; else, loc = 'northwest'; end
        lg = legend(ax, h, LAB, 'Interpreter','tex', 'Location',loc, 'NumColumns',1);
        lg.FontSize = ST.FSLEG;  lg.FontWeight = 'bold';
        lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = ST.LWBOX;  lg.Color = 'w';
        lg.ItemTokenSize = [55 25];
    end
    hold(ax,'off');
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 ST.CANV ST.CANV],'PaperSize',[ST.CANV ST.CANV]);
    print(fig, out, '-dpng', '-r200');   close(fig);
end

function [lim, tk] = axis_sym_(R)
%AXIS_SYM_  symmetric x axis [-R, R], three even ticks, ends padded by the step
    s = R/2;   lim = [-R R];   tk = [-s 0 s];
end

function [lim, tk] = ylim_odd_(maxv)
%YLIM_ODD_  figure-style rules 4+5 from zero: top = 4*s, ticks (1:3)*s, so neither
%   0 nor the top edge carries a tick, and tick spacing equals the end padding.
    s = nice_step_(maxv/4);
    lim = [0 4*s];   tk = (1:3)*s;
end

function [lim, tk] = ylim_sym_(maxabs)
%YLIM_SYM_  same rules, symmetric about zero: frame [-2s, 2s], ticks [-s 0 s].
    s = nice_step_(maxabs/2);
    lim = [-2*s 2*s];   tk = [-s 0 s];
end

function s = nice_step_(need)
%NICE_STEP_  smallest two-significant-digit step >= need, preferring a round one when
%   it wastes no more than 20 % (figure-style: avoids leaving the panel half empty).
    cand = [1 1.1 1.2 1.25 1.5 1.6 1.75 2 2.25 2.5 3 4 5 6 7.5 8];
    nice = [1 1.5 2 2.5 3 4 5];
    for kk = floor(log10(need)) : floor(log10(need))+1
        ok = cand(cand*10^kk >= need*(1-1e-12)) * 10^kk;
        if isempty(ok), continue, end
        s = min(ok);
        okn = nice(nice*10^kk >= need*(1-1e-12)) * 10^kk;
        if ~isempty(okn) && min(okn) <= 1.2*s, s = min(okn); end
        return
    end
    error('nice_step_:none','no step for %g', need);
end
