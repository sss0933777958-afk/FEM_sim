% nn_grid.m -- training positions of the NN inverse model in the actuator frame (3D scatter).
% Reads utils/data/long2016_hexapole_halfcut/nn_grid.mat (written by utils/scripts/long2016_hexapole_halfcut/nn_grid.m).
% Load -> plot only (plot-scripts-pure).  Single series -> no legend; the origin (sphere centre) is a black dot.
% ZL: draw only the layers z_a in ZL [um] (a few adjacent layers, joined by the vertical edges, so the
%     cubes show without points of many depths overlapping); ZL = [] draws every layer.
% Output: figures/long2016_hexapole_halfcut/hexapole/nn_grid[_all][SFX].png   (_all: every layer)
if ~exist('VIEW','var'), VIEW = [-37.5 30]; end   % off the lattice rows, so the points do not stack
if ~exist('ZL','var'), ZL = [-30 0 30]; end
if ~exist('SFX','var'), SFX = ''; end            % file-name suffix, e.g. '_side' for VIEW = [0 0]
clearvars -except VIEW ZL SFX;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(here));   % .../Force/Maxwell
S   = load(fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut', 'nn_grid.mat'));
FN  = 'nn_grid';  if isempty(ZL), FN = [FN '_all']; end
FN  = [FN SFX '.png'];
OUT = fullfile(CALROOT, 'figures', 'long2016_hexapole_halfcut', 'hexapole', FN);

% ---- style (figure-style.md) ---------------------------------------------------------
FS = 60;  FSLAB = 46;  LWBOX = 5.0;  CANV = 14.5;  MS = 9;  MSO = 22;
CB = [0.05 0.10 0.95];  GR = [0.55 0.55 0.55];  LWE = 1.5;    % lattice-edge colour and width
AXP = [0.14 0.12 0.78 0.80];
STEP = 90;                                     % tick step [um]: ticks -90 0 90, frame +-180 (end gap = step)
LIM = repmat([-2 2]*STEP, 3, 1);  TK = (-1:1)*STEP;

P = S.P_act;
if ~isempty(ZL), P = P(ismembertol(P(:,3), ZL, 1e-9, 'DataScale',1), :); end
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes(fig, 'Position', AXP);  hold(ax, 'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
% lattice edges: join each point to its neighbour one step H along x_a, y_a, z_a (presentation only)
K = round(P / S.H);                            % integer lattice index
E = eye(3);
for k = 1:3
    [hit, j] = ismember(K + E(k,:), K, 'rows');
    i = find(hit);  j = j(hit);
    xs = [P(i,1) P(j,1) nan(numel(i),1)].';  ys = [P(i,2) P(j,2) nan(numel(i),1)].';  zs = [P(i,3) P(j,3) nan(numel(i),1)].';
    plot3(ax, xs(:), ys(:), zs(:), '-', 'Color',GR, 'LineWidth',LWE);
end
plot3(ax, P(:,1), P(:,2), P(:,3), 'o', 'MarkerSize',MS, 'MarkerFaceColor',CB, 'MarkerEdgeColor',CB);
plot3(ax, 0, 0, 0, 'o', 'MarkerSize',MSO, 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
grid(ax, 'off');  box(ax, 'off');  daspect(ax, [1 1 1]);
xlim(ax, LIM(1,:));  ylim(ax, LIM(2,:));  zlim(ax, LIM(3,:));
view(ax, VIEW(1), VIEW(2));
set(ax, 'XTick',TK, 'YTick',TK, 'ZTick',TK, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, ...
        'TickLength',[0.025 0.025], 'XTickLabelRotation',0, 'YTickLabelRotation',0, 'ZTickLabelRotation',0);
draw_box_(ax, LIM, LWBOX);                     % 9 cube edges; the 3 at the camera corner dropped
ax.Clipping = 'off';
xlabel(ax, '$\mathbf{x_a\;(\mu m)}$', 'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{y_a\;(\mu m)}$', 'Interpreter','latex', 'FontSize',FSLAB);
zlabel(ax, '$\mathbf{z_a\;(\mu m)}$', 'Interpreter','latex', 'FontSize',FSLAB);

set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
print(fig, OUT, '-dpng', '-r200');
close(fig);
fprintf('wrote %s | %d of %d points (z_a = %s um)\n', OUT, size(P,1), size(S.P_act,1), mat2str(ZL));

% ---------------------------------------------------------------------------------------
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
