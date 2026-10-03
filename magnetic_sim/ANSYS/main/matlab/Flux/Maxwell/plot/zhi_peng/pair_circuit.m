%% pair_circuit.m -- magnetic-circuit side view of the Pair_pole fields with the calibrated charges
%  Loads data/zhi_peng/.mat/pair_slice.mat (arrows, plate outline; written by utils/zhi_peng/pair_slice.m)
%  and the two pair calibrations, then draws two square panels + one colorbar per figure by the nine figure
%  rules (tick numbers 60, box 5.0, horizontal endpoints numbered by hand and vertical ones not, three equally
%  spaced ticks with end gaps equal to the spacing, canvas in inches written with print).
%  Left panel = upper pole P2, right panel = lower pole P1; dashed arrow = actuator axis x_a.
%  Two figures (the colour scale is shared so they can be compared):
%    circuit_pair_P1.png    ONE excitation  : P1 excited, both panels show the P1 field;
%                           charges from calib_pair_P1_R150_xa
%    circuit_pair_P1P2.png  TWO excitations : right panel = P1 field, left panel = P2 field (each pole excited
%                           alone); charges from calib_pair_P1P2_R150_xa
%  Only the charge with bias e_hat is marked (pink): charge 1 (+u_hat side) in the P1 panel, charge 2
%  (-u_hat side) in the P2 panel, with its 3-D distance to the origin. Positions are projected on the
%  y = 0 plane (rho = x, z).
%  Output: figures/zhi_peng/
clear;
here = fileparts(mfilename('fullpath'));                   % .../plot/zhi_peng
CAL  = fileparts(fileparts(here));                         % .../Flux/Maxwell
MATD = fullfile(CAL, 'data', 'zhi_peng', '.mat');
OUT  = fullfile(CAL, 'figures', 'zhi_peng');
if ~exist(OUT, 'dir'), mkdir(OUT); end

LAB_L = [-0.70  0.02];                                     % distance-label position [mm], left panel: free air under the plate
LAB_R = [ 0.70 -0.02];                                     % right panel: free air over the plate (mirror image)

sl = load(fullfile(MATD, 'pair_slice.mat'));               % S1 right/P1-field, S2 left/P2-field, S3 left/P1-field
rA = load(fullfile(MATD, 'calib_pair_P1_R150_xa.mat'));    % one excitation
rB = load(fullfile(MATD, 'calib_pair_P1P2_R150_xa.mat'));  % two excitations
[CL, CTK] = cbar_clean([0 max([sl.S1.CLIM, sl.S2.CLIM, sl.S3.CLIM])], 6);   % one colour scale for both figures

cases = {'P1',   sl.S3, sl.S1, rA;
         'P1P2', sl.S2, sl.S1, rB};
for q = 1:size(cases, 1)
    SL = cases{q,2};   SR = cases{q,3};   r = cases{q,4};
    [SR.qbx, SR.qbz, SR.qbr] = charge_xz(r, 1);
    [SL.qbx, SL.qbz, SL.qbr] = charge_xz(r, 2);
    fprintf('%-5s P1 panel (%.3f, %.3f) mm |%.2f| um | P2 panel (%.3f, %.3f) mm |%.2f| um\n', cases{q,1}, ...
            SR.qbx, SR.qbz, SR.qbr, SL.qbx, SL.qbz, SL.qbr);
    render(SL, SR, CL, CTK, {LAB_L, LAB_R}, fullfile(OUT, ['circuit_pair_' cases{q,1} '.png']));
end

% ============================================================================
function render(SL, SR, CL, CTK, LAB, outp)
    FS = 60;  FSD = 45;  FSLAB = 50;  LWBOX = 5.0;
    PS = 12;  LM = 2.7;  MG = 2.9;  BM = 1.9;  TM = 0.6;       % panel side, margins [in]
    CG = 0.5; CW = 0.40; CR = 3.6;
    W  = LM + PS + MG + PS + CG + CW + CR;   H = BM + PS + TM;
    fig = figure('Color','w','Units','inches','Position',[0.3 0.3 W H],'Visible','off');
    SS = {SL, SR};   X0 = [LM, LM+PS+MG];
    for k = 1:2
        S  = SS{k};
        ax = axes(fig,'Units','inches','Position',[X0(k) BM PS PS]);
        try, ax.Toolbar = []; end                                %#ok<TRYNC>
        try, ax.Interactions = []; end                           %#ok<TRYNC>
        hold(ax,'on');
        patch(ax, S.pox, S.poz, [0.82 0.84 0.88], 'FaceAlpha',0.30, ...
              'EdgeColor',[0.28 0.30 0.36], 'LineWidth',4.0);
        nb = 28;  edges = linspace(0,CL(2),nb+1);  cmap = turbo(nb);  lw = linspace(1.0,3.6,nb);
        for q = 1:nb
            if q < nb, m = S.Bm_mT>=edges(q) & S.Bm_mT<edges(q+1); else, m = S.Bm_mT>=edges(q); end
            if any(m)
                quiver(ax, S.Xs(m), S.Zs(m), S.Uq(m), S.Wq(m), 0, 'Color',cmap(q,:), ...
                       'LineWidth',lw(q), 'MaxHeadSize',0.35);
            end
        end
        draw_axis_line(ax, S);
        set(findobj(ax,'Type','line'), 'LineWidth', 4.5);        % the axis line and its head
        plot(ax, 0, 0, '+', 'Color','k', 'MarkerSize',44, 'LineWidth',5.5);
        CP = [1 0.30 0.65];
        plot(ax, S.qbx, S.qbz, 'o', 'MarkerFaceColor',CP, 'MarkerEdgeColor','k', 'MarkerSize',36, 'LineWidth',3);
        % distance to the origin (3-D, um); no background box
        text(ax, LAB{k}(1), LAB{k}(2), sprintf('%.2f \\mum', S.qbr), ...
             'Color',CP*0.85, 'FontSize',FSD, 'FontWeight','bold', 'Interpreter','tex', ...
             'HorizontalAlignment','center', 'VerticalAlignment','middle');
        axis(ax,'equal');  xlim(ax,S.XL);  ylim(ax,S.ZL);
        box(ax,'on');  grid(ax,'off');
        set(ax, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, ...
                'TickDir','in', 'TickLength',[.018 .018], 'XTick',S.XT, 'YTick',S.ZT, ...
                'Units','inches', 'Position',[X0(k) BM PS PS]);
        colormap(ax, turbo);  clim(ax, CL);
        yoff = S.ZL(1) - 0.022*diff(S.ZL);
        for xv = S.XL
            text(ax, xv, yoff, sprintf('%g',xv), 'HorizontalAlignment','center', ...
                 'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
        end
        hold(ax,'off');
    end
    cb = colorbar(ax, 'Units','inches');
    cb.Position = [X0(2)+PS+CG, BM, CW, PS];   ax.Position = [X0(2) BM PS PS];
    cb.FontSize = FS;  cb.FontWeight = 'bold';  cb.LineWidth = LWBOX;
    cb.Limits = CL;    cb.Ticks = CTK;
    cb.Label.Interpreter = 'latex';  cb.Label.String = '$\mathbf{\|b\|\;(mT)}$';
    cb.Label.FontSize = FSLAB;
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 W H], 'PaperSize',[W H]);
    print(fig, outp, '-dpng', '-r100');
    close(fig);
    fprintf('wrote %s\n', outp);
end

% ============================================================================
function [x, z, dist] = charge_xz(r, k)
% Charge k of a pair calibration -> (rho, z) of the y = 0 plane [mm] and its 3-D distance to the origin [um].
%   actuator frame: p_act = l_hat*(+/-u_hat + e_hat_k)  (+u_hat for charge 1, -u_hat for charge 2);
%   measure/WP frame: p = R_act' * p_act.
    s  = 3 - 2*k;
    pa = r.l_hat * (s * r.u_hat(:) + r.e_hat(:,k));
    pm = r.R_act.' * pa;
    x = pm(1)*1e3;   z = pm(3)*1e3;   dist = norm(pa)*1e6;
end

% ============================================================================
function draw_axis_line(ax, S)
% actuator x_a axis, dashed, through the origin along the in-plane projection of the P1 axis, clipped to the
% frame; the arrowhead sits at the +x_a end, fully inside the frame.
    d = S.axd(:);
    ts = [];
    if abs(d(1))>1e-9, ts=[ts S.XL(1)/d(1) S.XL(2)/d(1)]; end
    if abs(d(2))>1e-9, ts=[ts S.ZL(1)/d(2) S.ZL(2)/d(2)]; end
    tin = [];
    for t = ts
        p = t*d;
        if p(1)>=S.XL(1)-1e-6 && p(1)<=S.XL(2)+1e-6 && p(2)>=S.ZL(1)-1e-6 && p(2)<=S.ZL(2)+1e-6
            tin(end+1) = t; %#ok<AGROW>
        end
    end
    if numel(tin) < 2, return; end
    tlo = min(tin);  thi = max(tin);
    plo = tlo*d;  pe = thi*d;
    col = [0.20 0.30 0.45];
    hl = 0.06*diff(S.XL);  a = 26*pi/180;  b = -d;
    Rm = @(th) [cos(th) -sin(th); sin(th) cos(th)];
    Wg = [Rm(a)*b, Rm(-a)*b];
    mg = 0.012*diff(S.XL);
    for it = 1:400
        tp = [pe, pe + Wg*hl];
        if all(tp(1,:) >= S.XL(1)+mg) && all(tp(1,:) <= S.XL(2)-mg) && ...
           all(tp(2,:) >= S.ZL(1)+mg) && all(tp(2,:) <= S.ZL(2)-mg), break; end
        pe = pe - 0.01*hl*d;
    end
    plot(ax,[plo(1) pe(1)],[plo(2) pe(2)],'--','Color',col,'LineWidth',2.5);
    for k = 1:2
        w = Wg(:,k);
        plot(ax,[pe(1) pe(1)+w(1)*hl],[pe(2) pe(2)+w(2)*hl],'-','Color',col,'LineWidth',2.5);
    end
end

% ============================================================================
function [cl2, tk] = cbar_clean(cl, ntgt)
% Clean colour-axis bounds with equally spaced ticks (both ends numbered): widen [lo hi] to an integer
% multiple of a nice step, number of intervals closest to ntgt.
    cand = [0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000];
    bs = [];   bd = inf;
    for s = cand
        lo = floor(cl(1)/s + 1e-9)*s;   hi = ceil(cl(2)/s - 1e-9)*s;
        n  = round((hi-lo)/s);
        if n >= 3 && n <= 8 && abs(n-ntgt) < bd,  bd = abs(n-ntgt);  bs = [lo hi s];  end
    end
    if isempty(bs), bs = [cl(1) cl(2) (cl(2)-cl(1))/ntgt]; end
    tk  = round((bs(1):bs(3):bs(2))/0.1)*0.1;
    cl2 = [tk(1) tk(end)];
end
