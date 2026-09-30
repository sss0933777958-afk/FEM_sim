%% err_map.m -- force-error magnitude on the three planes through the origin.
%  Reads emap_R<R>_P<tag>_<model tag>.mat written by utils/force_map.m and draws
%  two figures that share ONE colour scale:
%      err_map_<..>_xy_xz.png   x_a-y_a (left) and x_a-z_a (right)
%      err_map_<..>_yz.png      y_a-z_a
%  No fitting, no field reading: load -> pick -> plot (plot-scripts-pure).
%
%  STYLE (user's call, 2026-09-29): the polar layout of the gain / isotropy
%  figures (paper_fig_plot/plot_svd_polar.m) -- no rectangular frame, no tick
%  numbers and no axis titles; the disc carries rings and 30-degree spokes, the
%  poles are named on the rim, two arrows from the origin show the in-plane
%  axes, colormap jet. The disc radius is not written on the figure; it goes in
%  the caption.
%  Figure rules that still apply: numbers 60 and bold (colorbar and pole names),
%  colorbar numbered at both ends with even steps, every panel square, canvas
%  in inches and at most 14.5 in tall, written with print.
clc;
if exist('PLOT_OVR','var'), OVR_ = PLOT_OVR; else, OVR_ = struct(); end
clearvars -except OVR_;
RSEL = 150;   PTAG = 'all';   TAG = 'eighteen';
if isfield(OVR_,'RSEL'), RSEL = OVR_.RSEL; end
if isfield(OVR_,'PTAG'), PTAG = OVR_.PTAG; end       % 'all' | '1' .. '6'
if isfield(OVR_,'TAG'),  TAG  = OVR_.TAG;  end
% PLOT_OVR.CCAP = upper colour limit [pN]. Values above it are drawn in the top
% colour. Needed once the disc reaches the pole tips: a few points next to a tip
% carry errors hundreds of times the rest and would leave the disc one colour.
CCAP = [];   if isfield(OVR_,'CCAP'), CCAP = OVR_.CCAP; end
MODEL = 'long2016_hexapole_halfcut';
here  = fileparts(mfilename('fullpath'));
FMX   = fileparts(fileparts(here));
S     = load(fullfile(FMX,'data',MODEL,'.mat',sprintf('emap_R%d_P%s_%s.mat',RSEL,PTAG,TAG)));
OUT   = fullfile(FMX,'figures',MODEL);

FS = 60;  FSLAB = 50;  LWBOX = 5.0;
RIN = 4.1;                        % disc radius on paper [in]
PS  = 12;                         % square panel side [in]; holds the pole names
BM  = 0.7;  TM = 0.7;  LM = 0.2;  CG = 0.2;  CW = 0.40;  CR = 3.3;
GRY = [.35 .35 .35];   PCOL = [0.5 0 0.12];

% ---- one colour scale for all three planes ----------------------------------
if isempty(CCAP), [CL, CTK] = cbar_ticks_([S.emin S.emax]);
else,             [CL, CTK] = cbar_ticks_([0 CCAP]);
end
R   = S.R;
LIM = R * (PS/2) / RIN;           % data half-width of the panel
U   = R / RIN;                    % data units per inch
fprintf('colour scale %g .. %g pN (data %.4f .. %.4f pN)%s', CL, S.emin, S.emax, newline);

% in-plane axes -> the poles at +/- ends (the ideal charge grid is +/- the axes)
POLE = struct('x_a',[1 2], 'y_a',[3 4], 'z_a',[5 6]);

JOBS = {{'xy','xz'}, {'yz'}};
for q = 1:numel(JOBS)
    J  = JOBS{q};   np = numel(J);
    W  = LM + np*PS + CG + CW + CR;   H = BM + PS + TM;
    fig = figure('Color','w','Units','inches','Position',[0.3 0.3 W H],'Visible','off');
    drawnow;
    assert(max(abs(fig.Position(3:4) - [W H])) < 0.01, 'err_map:canvas', ...
           'the canvas came out %.2f x %.2f in instead of %.2f x %.2f', fig.Position(3:4), W, H);
    for k = 1:np
        D  = S.(J{k});
        x0 = LM + (k-1)*PS;
        ax = axes(fig,'Units','inches','Position',[x0 BM PS PS]);
        try, ax.Toolbar = []; end                                    %#ok<TRYNC>
        try, ax.Interactions = []; end                               %#ok<TRYNC>
        hold(ax,'on');
        % E is indexed (first axis, second axis); imagesc wants rows = vertical
        im = imagesc(ax, S.g, S.g, D.E.');
        im.AlphaData = ~isnan(D.E.');
        colormap(ax, jet(256));   clim(ax, CL);

        th = linspace(0, 2*pi, 721);
        for rr = (1:4)*R/4
            plot(ax, rr*cos(th), rr*sin(th), '-', 'Color',GRY, 'LineWidth',2.4);
        end
        for aa = 0:30:330
            plot(ax, [0 R*cosd(aa)], [0 R*sind(aa)], '-', 'Color',GRY, 'LineWidth',2.0);
        end
        arrow_(ax, [1 0], 0.45*R);
        arrow_(ax, [0 1], 0.45*R);

        % pole names on the rim: marker on the circle, boxed name outside it.
        % The name is pushed out by its own half-extent, which differs between
        % the horizontal and the vertical direction at this type size.
        GAP = 0.30;   HWX = 0.80;   HWY = 0.50;                       % [in]
        pk  = [POLE.(D.ax{1}); POLE.(D.ax{2})];                      % [+h -h; +v -v]
        dir = [1 0; -1 0; 0 1; 0 -1];   nm = [pk(1,1) pk(1,2) pk(2,1) pk(2,2)];
        for m = 1:4
            d = dir(m,:);
            plot(ax, R*d(1), R*d(2), 'o', 'MarkerSize',26, 'MarkerFaceColor',[.90 .90 .90], ...
                 'MarkerEdgeColor',PCOL, 'LineWidth',4.0);
            ro = R + U*(GAP + abs(d(1))*HWX + abs(d(2))*HWY);
            text(ax, ro*d(1), ro*d(2), sprintf('P%d', nm(m)), 'FontSize',FS, ...
                 'FontWeight','bold', 'Color',PCOL, 'HorizontalAlignment','center', ...
                 'VerticalAlignment','middle', 'BackgroundColor','w', 'Margin',4, ...
                 'EdgeColor',PCOL, 'LineWidth',2.5);
        end

        axis(ax,'xy');   axis(ax,'off');
        xlim(ax, [-LIM LIM]);   ylim(ax, [-LIM LIM]);
        set(ax, 'DataAspectRatio',[1 1 1], 'Units','inches', 'Position',[x0 BM PS PS]);
        hold(ax,'off');
    end
    cb = colorbar(ax, 'Location','manual', 'Units','inches');
    cb.FontSize = FS;   cb.FontWeight = 'bold';   cb.LineWidth = LWBOX;
    cb.Limits = CL;     cb.Ticks = CTK;
    cb.Label.Interpreter = 'latex';   cb.Label.FontSize = FSLAB;
    cb.Label.String = '$\mathbf{Error\;(pN)}$';
    % positions last, the colorbar after the axes: styling re-lays it out
    CBP = [x0+PS+CG, BM + PS/2 - RIN, CW, 2*RIN];     % as tall as the disc
    ax.Position = [x0 BM PS PS];   drawnow;
    cb.Position = CBP;             drawnow;
    assert(max(abs(cb.Position - CBP)) < 1e-6 && ...
           max(abs(ax.Position - [x0 BM PS PS])) < 1e-6, 'err_map:layout', ...
           'the colorbar or the panel moved off its assigned position');

    ct = '';   if ~isempty(CCAP), ct = sprintf('_c%g', CCAP); end
    name = fullfile(OUT, sprintf('err_map_R%d_P%s_%s%s.png', RSEL, PTAG, strjoin(J,'_'), ct));
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 W H], 'PaperSize',[W H]);
    print(fig, name, '-dpng', '-r100');
    close(fig);
    fprintf('saved %s%s', name, newline);
end

% ---- arrow from the origin along d, drawn by hand so its head keeps its shape
function arrow_(ax, d, Lx)
    d = d(:).'/norm(d);   n = [-d(2) d(1)];   tip = Lx*d;
    hl = 0.15*Lx;   hw = 0.06*Lx;   b = tip - hl*d;
    plot(ax, [0 tip(1)], [0 tip(2)], '-', 'Color','k', 'LineWidth',6);
    plot(ax, [tip(1) b(1)+hw*n(1)], [tip(2) b(2)+hw*n(2)], '-', 'Color','k', 'LineWidth',6);
    plot(ax, [tip(1) b(1)-hw*n(1)], [tip(2) b(2)-hw*n(2)], '-', 'Color','k', 'LineWidth',6);
end

% ---- colour limits: both ends numbered, even steps, smallest span that covers
%      the data (the rule of plot_svd_polar.m's cbar_ticks, with finer steps
%      added because the error is below 1 pN)
function [cl2, tk] = cbar_ticks_(cl)
    cand = [0.01 0.02 0.03 0.05 0.1 0.2 0.3 0.5 1 2 3 5 10 20 30 50 100 200 300 500 1000];
    bs = [];   bsp = inf;
    for s = cand
        lo = floor(cl(1)/s + 1e-9)*s;   hi = ceil(cl(2)/s - 1e-9)*s;
        n  = round((hi-lo)/s);
        if n < 3 || n > 6, continue; end
        if (hi - lo) < bsp - 1e-9,  bsp = hi - lo;  bs = [lo hi s];  end
    end
    assert(~isempty(bs), 'err_map:cbar', 'no tick step fits [%g %g]', cl(1), cl(2));
    tk  = round((bs(1):bs(3):bs(2))/0.01)*0.01;
    cl2 = [tk(1) tk(end)];
end
