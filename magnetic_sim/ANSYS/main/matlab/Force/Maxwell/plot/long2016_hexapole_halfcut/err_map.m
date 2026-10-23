% err_map.m -- relative force error of the fitted LS inverse model on one mid-plane, inside the fitting nodes (+-30 um),
% one direction (PH, TH) [deg] (continuous map, interpolated shading).  The origin is the black dot.
% Reads utils/data/long2016_hexapole_halfcut/err_map_p<PH>_t<TH>_<PL>.mat (utils/scripts/.../err_map.m).  Load -> plot only.
% Output: figures/long2016_hexapole_halfcut/hexapole/err_map_p<PH>_t<TH>_<PL>.png
% MODE = 'dir': 3D surface of e(theta, phi) at one position (err_map_dir.mat) -> err_map_dir.png
% MODE = 'loss': NN training loss per epoch (file NNF, default nn_phys.mat) up to epoch NE ([] = all, 'conv' = the
%               saved convergence epoch ep_conv: the horizontal axis then ends there) -> nn_loss.png
% MODE = 'R'  : at each of the 27 fitting nodes the largest e over the 2522 directions vs R = |p|, LS (ls_fit_chain.mat,
%               chain-label fit) and NN (file NNF, physics loss) overlaid -> err_R.png
if ~exist('MODE','var'), MODE = 'plane'; end
if ~exist('VIEW','var'), VIEW = [-37.5 30]; end
if ~exist('PH','var'), PH = 0; end
if ~exist('TH','var'), TH = 0; end
if ~exist('PL','var'), PL = 'yz'; end
if ~exist('NE','var'), NE = []; end
if ~exist('NNF','var'), NNF = 'nn_phys.mat'; end
clearvars -except PH TH PL MODE VIEW NE NNF;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(here));   % .../Force/Maxwell
FIGD = fullfile(CALROOT, 'figures', 'long2016_hexapole_halfcut', 'hexapole');

% ---- style (figure-style.md) ---------------------------------------------------------
FS = 60;  FSLAB = FS;  LWBOX = 5.0;  CANV = 14.5;  MSO = 22;

if strcmp(MODE, 'loss')                        % ---- NN (physics loss) training loss per epoch, log scale ----
    NN  = load(fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut', NNF));
    Le  = NN.loss(:,1);  ne = numel(Le);
    CB  = [0.05 0.10 0.95];
    nice = sort(reshape([1 1.5 2 2.5 3 4 5].' .* 10.^(0:4), [], 1));
    sx  = min(nice(nice >= ne/4));  XL = [0 4*sx];  XT = sx*(1:3);            % three ticks, end gap = step
    if strcmp(NE, 'conv')                                                      % axis ends at the convergence epoch
        ne = NN.ep_conv;  Le = Le(1:ne);  XL = [0 ne];  XT = round(ne*(1:3)/4);
    elseif ~isempty(NE)                                                        % NE: last epoch drawn
        ne = min(NE, ne);  Le = Le(1:ne);  sx = min(nice(nice >= ne/4));  XL = [0 4*sx];  XT = sx*(1:3);
    end
    YL  = [1e-4 1e0];  YT = [1e-3 1e-2 1e-1];                                  % decades, end gap = one decade
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes(fig, 'Units','normalized', 'Position',[0.215 0.175 0.700 0.725]);  hold(ax, 'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    plot(ax, 1:ne, Le, '-', 'LineWidth',5, 'Color',CB);
    set(ax, 'YScale','log', 'XLim',XL, 'YLim',YL, 'XTick',XT, 'YTick',YT, 'FontSize',FS, 'FontWeight','bold', ...
            'LineWidth',LWBOX, 'Box','on', 'TickDir','in', 'Layer','top', 'YMinorTick','off');
    xlabel(ax, '$\mathbf{Epoch}$', 'Interpreter','latex', 'FontSize',FSLAB);
    ylabel(ax, '$\mathbf{Loss}$', 'Interpreter','latex', 'FontSize',FSLAB);
    for x = XL                                                                 % horizontal-axis end numbers
        text(ax, x, YL(1) / 10^(0.022*log10(YL(2)/YL(1))), num2str(x), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold');
    end
    OUT = fullfile(FIGD, 'nn_loss.png');
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
    print(fig, OUT, '-dpng', '-r200');
    close(fig);
    fprintf('wrote %s | %d epochs, loss %.4g -> %.4g\n', OUT, ne, Le(1), Le(end));
    return
end

if strcmp(MODE, 'R')                           % ---- worst direction at every fitting node vs distance R ----
    DD  = fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut');
    LS  = load(fullfile(DD, 'ls_fit_chain.mat'), 'e_ij', 'pbar', 'l_hat');      % LS, chain labels
    NN  = load(fullfile(DD, NNF), 'e_ij', 'P');                                % NN, physics loss
    RL  = vecnorm(LS.pbar, 2, 2) * LS.l_hat;   EL = max(LS.e_ij, [], 1).';      % R [um], e of the worst direction [%]
    RN  = vecnorm(NN.P, 2, 2);                 EN = max(NN.e_ij, [], 1).';
    CB  = [0.05 0.10 0.95];  CR = [0.85 0.10 0.10];  FSLEG = 45;  MSZ = 30;
    XT  = 15:15:45;  XL = [0 60];                                              % three ticks, end gap = step
    nice = sort(reshape([1 1.5 2 2.5 3 4 5].' .* 10.^(-2:3), [], 1));
    sy  = min(nice(nice >= 1.02*max([EL; EN])/4));  YT = sy*(1:3);  YL = [0 4*sy];
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes(fig, 'Units','normalized', 'Position',[0.215 0.175 0.700 0.725]);  hold(ax, 'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    h(1) = plot(ax, RL, EL, 'o', 'MarkerSize',MSZ, 'LineWidth',4, 'Color',CB);
    h(2) = plot(ax, RN, EN, 's', 'MarkerSize',MSZ, 'LineWidth',4, 'Color',CR);
    set(ax, 'XLim',XL, 'YLim',YL, 'XTick',XT, 'YTick',YT, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, ...
            'Box','on', 'TickDir','in', 'Layer','top');
    xlabel(ax, '$\mathbf{R\;(\mu m)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    ylabel(ax, '$\mathbf{e\;(\%)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    for x = XL                                                                 % horizontal-axis end numbers
        text(ax, x, YL(1) - 0.022*diff(YL), num2str(x), 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
             'FontSize',FS, 'FontWeight','bold');
    end
    lg = legend(ax, h, {'LS', 'NN'}, 'FontSize',FSLEG, 'FontWeight','bold', 'Location','northeast');
    lg.LineWidth = LWBOX;  lg.ItemTokenSize = [55 25];
    OUT = fullfile(FIGD, 'err_R.png');
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
    print(fig, OUT, '-dpng', '-r200');
    close(fig);
    fprintf('wrote %s | LS: max %.3g %%, mean %.3g %% | NN: max %.3g %%, mean %.3g %% (27 nodes, worst direction)\n', ...
            OUT, max(EL), mean(EL), max(EN), mean(EN));
    return
end

if strcmp(MODE, 'dir')                         % ---- e over every direction at one position (3D) ----
    D   = load(fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut', 'err_map_dir.mat'));
    [TT, PP] = meshgrid(D.th, D.ph);
    SZ  = 1;  while 4*SZ < max(D.ED(:)), SZ = SZ * 2; end     % e ticks SZ, 2SZ, 3SZ; frame [0 4SZ]
    LIM = [0 360; -90 90; 0 4*SZ];
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes(fig, 'Position', [0.17 0.16 0.72 0.74]);  hold(ax, 'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    surf(ax, TT, PP, D.ED, D.ED, 'EdgeColor','none', 'FaceColor','interp');
    colormap(ax, jet(256));  clim(ax, LIM(3,:));
    grid(ax, 'off');  box(ax, 'off');  pbaspect(ax, [1 1 0.8]);
    xlim(ax, LIM(1,:));  ylim(ax, LIM(2,:));  zlim(ax, LIM(3,:));
    view(ax, VIEW(1), VIEW(2));
    set(ax, 'XTick',90:90:270, 'YTick',-45:45:45, 'ZTick',SZ*(1:3), 'FontSize',FS, 'FontWeight','bold', ...
            'LineWidth',LWBOX, 'TickLength',[0.025 0.025], 'XTickLabelRotation',0, 'YTickLabelRotation',0, ...
            'ZTickLabelRotation',0);
    draw_box_(ax, LIM, LWBOX);
    ax.Clipping = 'off';
    xlabel(ax, '$\mathbf{\theta\;(deg)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    ylabel(ax, '$\mathbf{\varphi\;(deg)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    zlabel(ax, '$\mathbf{e\;(\%)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    OUT = fullfile(FIGD, 'err_map_dir.png');
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
    print(fig, OUT, '-dpng', '-r200');
    close(fig);
    fprintf('wrote %s | p = [%g %g %g] um | e %.3g..%.3g %%\n', OUT, D.P0, min(D.ED(:)), max(D.ED(:)));
    return
end

D   = load(fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut', sprintf('err_map_p%g_t%g_%s.mat', PH, TH, PL)));
AXP = [0.20 0.28 0.56 0.56];                  % square plot box; room for the colour bar on the right
LIM = [-1 1] * D.A;  TK = (-1:1) * D.A/2;      % frame = outermost fitting nodes, three ticks, end gap = tick step

NAME = {'xy', 'xz', 'yz'};  AX = [1 2; 1 3; 2 3];  AXN = {'x_a', 'y_a', 'z_a'};
[lo, hi, sc] = cbar_ticks_([min(D.E(:)) max(D.E(:))]);
q = find(strcmp(NAME, D.PL));
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes(fig, 'Units','normalized', 'Position', AXP);  hold(ax, 'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
[UU, VV] = meshgrid(D.u, D.u);                  % interpolated shading: continuous map, no cell blocks
surface(ax, UU, VV, zeros(size(UU)) - 1, D.E, 'EdgeColor','none', 'FaceColor','interp');
plot(ax, 0, 0, 'o', 'MarkerSize',MSO, 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
colormap(ax, jet(256));  clim(ax, [lo hi]);
cb = colorbar(ax);  cb.Limits = [lo hi];  cb.Ticks = lo:sc:hi;  cb.LineWidth = LWBOX;
cb.FontSize = FS;  cb.FontWeight = 'bold';
ylabel(cb, '$\mathbf{e\;(\%)}$', 'Interpreter','latex', 'FontSize',FSLAB);
xlim(ax, LIM);  ylim(ax, LIM);
set(ax, 'XTick',TK, 'YTick',TK, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, 'Box','on', 'Layer','top');
ax.Position = AXP;
xlabel(ax, ['$\mathbf{' AXN{AX(q,1)} '\;(\mu m)}$'], 'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, ['$\mathbf{' AXN{AX(q,2)} '\;(\mu m)}$'], 'Interpreter','latex', 'FontSize',FSLAB);
yl = ylim(ax);
for x = LIM                                     % horizontal-axis end numbers (not in XTick)
    text(ax, x, yl(1) - 0.012*diff(yl), num2str(x), 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
         'FontSize',FS, 'FontWeight','bold');
end
OUT = fullfile(FIGD, sprintf('err_map_p%g_t%g_%s.png', PH, TH, PL));
set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
print(fig, OUT, '-dpng', '-r200');
close(fig);
fprintf('wrote %s | plane %s | e %.3g..%.3g %%\n', OUT, D.PL, min(D.E(:)), max(D.E(:)));

% ---------------------------------------------------------------------------------------
function [lo, hi, s] = cbar_ticks_(r)
% colour-bar range and step: 3..4 segments, narrowest range; on a tie the larger step (fewer segments)
    best = Inf;
    for s0 = reshape([1 2 3 4 5].' * 10.^(-3:2), 1, [])
        s0 = round(s0, 10);
        l = round(floor(r(1)/s0)*s0, 10);  h = round(ceil(r(2)/s0)*s0, 10);  m = round((h - l)/s0);
        if m >= 3 && m <= 4 && h - l <= best + 1e-12, best = h - l;  lo = l;  hi = h;  s = s0; end
    end
end

function draw_box_(ax, lim, lw)
% box off + the 9 cube edges drawn by hand (the 3 edges meeting the corner nearest the camera are dropped).
    [Xc,Yc,Zc] = ndgrid(lim(1,:), lim(2,:), lim(3,:));  C = [Xc(:) Yc(:) Zc(:)];
    E = [];
    for i = 1:8
        for j = i+1:8
            if nnz(abs(C(i,:)-C(j,:)) > 0) == 1, E = [E; i j];  end %#ok<AGROW>
        end
    end
    cp = campos(ax);  [~, near] = min(sum((C - cp).^2, 2));
    E = E(~(E(:,1)==near | E(:,2)==near), :);
    for k = 1:size(E,1)
        p1 = C(E(k,1),:);  p2 = C(E(k,2),:);
        plot3(ax, [p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'k-', 'LineWidth', lw);
    end
end
