% f_shell.m -- 3-D force envelope at the origin (actuator frame): for every direction r_hat the largest force
% component that a unit-norm current (||I||_2 = 1) can produce along r_hat, lambda_max(Q_r), Q_r = sum_k r_hat_k Q_k,
% drawn at lambda_max * r_hat.  lambda_max is continuous in r_hat, so the surface is closed (no tears; plotting the
% full force of the top eigenvector instead tears where the two largest eigenvalues cross).
% Loads utils/data/long2016_hexapole_halfcut/f_shell.mat only (written by utils/scripts/long2016_hexapole_halfcut/f_shell.m).
% Mesh lines every MSTEP deg in phi and theta; surface coloured by the force magnitude |f| = lambda_max [pN] (jet),
% colour bar styled as err_map.m.
% FRAME = 'global' (default): rotated to the global (Maxwell) frame, f = R_act' * f_a, so the envelope stands tilted as
%         the pole axes do (each 54.74 deg from z) -> f_shell.png;  'actuator': axes f_{x_a}, f_{y_a}, f_{z_a} -> f_shell_act.png
% Output: figures/long2016_hexapole_halfcut/hexapole/f_shell[_act].png
% TAG = position tag of the input (utils/data/long2016_hexapole_halfcut/f_shell_<TAG>.mat -> f_shell_<TAG>[_act].png); '' = origin.
% SRC = 'lambda' (default, above) | 'rho': maximum force gain rho(p, r_hat) = 1 / ||I_opt_hat||^2 [pN/A^2] of the
%       minimum-norm unit-current problem (utils/data/long2016_hexapole_halfcut/rho_map[_<TAG>].mat, all directions global minima), drawn as the
%       vector rho * r_hat in the actuator frame (axes rho_{x_a}, rho_{y_a}, rho_{z_a}; FRAME ignored) -> rho[_<TAG>].png
if ~exist('FRAME','var'), FRAME = 'global'; end
if ~exist('TAG','var'),   TAG   = ''; end
if ~exist('SRC','var'),   SRC   = 'lambda'; end
% VIEW = [az el] camera; VTAG = output suffix for a non-default view (e.g. 'bot' -> rho_bot.png)
if ~exist('VIEW','var'),  VIEW  = [65 25]; end
if ~exist('VTAG','var'),  VTAG  = ''; end
clearvars -except FRAME TAG SRC VIEW VTAG;  close all;
SFX = '';  if ~isempty(TAG), SFX = ['_' TAG]; end
OSFX = SFX;  if ~isempty(VTAG), OSFX = [SFX '_' VTAG]; end   % output name only
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(here));   % .../Force/Maxwell
if strcmp(SRC, 'rho')
    R   = load(fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut', ['rho_map' SFX '.mat']), 'RHO', 'ph', 'th', 'GLOB', 'R_act');
    assert(all(R.GLOB(:)), 'f_shell:glob', 'rho_map.mat has directions not certified as global minima');
    [S.TH, S.PHI] = meshgrid(R.th, R.ph);  S.LAM = R.RHO;  S.R_act = R.R_act;   % rho [pN/A^2]
    FRAME = 'actuator';  UNIT = 'pN/A^2';  CBL = '\rho';  AXP = [0.18 0.13 0.55 0.76];
    OUT = fullfile(CALROOT, 'figures', 'long2016_hexapole_halfcut', 'hexapole', ['rho' OSFX '.png']);
else
    S   = load(fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut', ['f_shell' SFX '.mat']));
    FN  = struct('global',['f_shell' OSFX '.png'], 'actuator',['f_shell' OSFX '_act.png']);
    OUT = fullfile(CALROOT, 'figures', 'long2016_hexapole_halfcut', 'hexapole', FN.(FRAME));
    UNIT = 'pN';  CBL = '|f|';  AXP = [0.18 0.13 0.55 0.76];
end

% ---- style (figure-style.md) ---------------------------------------------------------
FS = 60;  FSLAB = 46;  LWBOX = 5.0;  CANV = 14.5;  LWM = 1.2;  MSO = 40;
MSTEP = 10;                                   % mesh spacing [deg]

LAM = S.LAM;                                  % lambda_max [pN]
X = LAM .* cosd(S.PHI) .* cosd(S.TH);  Y = LAM .* cosd(S.PHI) .* sind(S.TH);  Z = LAM .* sind(S.PHI);   % actuator frame
if strcmp(FRAME, 'global')                    % presentation rotation only: f = R_act' * f_a
    G = S.R_act.' * [X(:) Y(:) Z(:)].';
    X = reshape(G(1,:), size(X));  Y = reshape(G(2,:), size(X));  Z = reshape(G(3,:), size(X));
    AX = {'f_{x}', 'f_{y}', 'f_{z}'};
elseif strcmp(SRC, 'rho')
    AX = {'\rho_{x_a}', '\rho_{y_a}', '\rho_{z_a}'};
else
    AX = {'f_{x_a}', 'f_{y_a}', 'f_{z_a}'};
end
[LIM, TK] = axis_ticks_({X, Y, Z});           % equal span and step on all three axes, centred on the data
[CL, CTK] = cbar_ticks_([min(LAM(:)) max(LAM(:))]);   % colour range from the data (same helper as err_map.m)

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes(fig, 'Position', AXP);  hold(ax, 'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
surf(ax, X, Y, Z, LAM, 'FaceColor','interp', 'EdgeColor','none', 'FaceLighting','none');   % no shading: colour = |f| only
% mesh lines: constant phi (rows) and constant theta (columns)
ip = find(mod(S.PHI(:,1), MSTEP) == 0).';  it = find(mod(S.TH(1,:), MSTEP) == 0);
for i = ip, plot3(ax, X(i,:), Y(i,:), Z(i,:), 'k-', 'LineWidth',LWM); end
for j = it, plot3(ax, X(:,j), Y(:,j), Z(:,j), 'k-', 'LineWidth',LWM); end
grid(ax, 'off');  box(ax, 'off');  daspect(ax, [1 1 1]);
xlim(ax, LIM(1,:));  ylim(ax, LIM(2,:));  zlim(ax, LIM(3,:));
view(ax, VIEW(1), VIEW(2));                           % default [65 25] = same view as sph_shell.m
set(ax, 'XTick',TK(1,:), 'YTick',TK(2,:), 'ZTick',TK(3,:), 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, ...
        'TickLength',[0.025 0.025], 'XTickLabelRotation',0, 'YTickLabelRotation',0, 'ZTickLabelRotation',0);
draw_box_(ax, LIM, LWBOX);                         % 9 cube edges; the 3 at the camera corner dropped
ax.Clipping = 'off';
xlabel(ax, ['$\mathbf{' AX{1} '\;(' UNIT ')}$'], 'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, ['$\mathbf{' AX{2} '\;(' UNIT ')}$'], 'Interpreter','latex', 'FontSize',FSLAB);
zlabel(ax, ['$\mathbf{' AX{3} '\;(' UNIT ')}$'], 'Interpreter','latex', 'FontSize',FSLAB);
drawnow;  lb = ax.XLabel;  lb.Units = 'normalized';        % a long x title must not run off the left edge
xl = AXP(1) + lb.Extent(1) * AXP(3);                       % its left end, figure-normalized
if xl < 0.01, lb.Position(1) = lb.Position(1) + (0.01 - xl) / AXP(3); end
colormap(ax, jet(256));  clim(ax, CL);
cb = colorbar(ax, 'Position',[0.78 0.22 0.032 0.58]);
cb.FontSize = FS;  cb.FontWeight = 'bold';  cb.LineWidth = LWBOX;
cb.Limits = CL;  cb.Ticks = CTK;
cb.Label.Interpreter = 'latex';  cb.Label.String = ['$\mathbf{' CBL '\;(' UNIT ')}$'];  cb.Label.FontSize = FSLAB;

set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
print(fig, OUT, '-dpng', '-r200');
close(fig);
fprintf('wrote %s | lambda_max %.3f .. %.3f pN | colour %s step %g\n', OUT, min(LAM(:)), max(LAM(:)), mat2str(CL), CTK(2)-CTK(1));

% ---------------------------------------------------------------------------------------
function draw_box_(ax, lim, lw)
% box off + the 9 cube edges drawn by hand (the 3 edges meeting the corner nearest the camera are dropped),
% same as sph_shell.m; lim = 3x2 axis limits.
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

function [cl2, tk] = cbar_ticks_(cl)
% Colour-bar range and ticks from the data range (copied from err_map.m): the step with 3..6 equal segments whose
% [floor(lo/s)*s, ceil(hi/s)*s] is tightest; both ends are ticks.
    cand = [0.01 0.02 0.03 0.05 0.1 0.2 0.3 0.5 1 2 3 5 10 20 30 50 100 200 300 500 1000];
    bs = [];   bsp = inf;
    for s = cand
        lo = floor(cl(1)/s + 1e-9)*s;   hi = ceil(cl(2)/s - 1e-9)*s;
        n  = round((hi-lo)/s);
        if n < 3 || n > 6, continue; end
        if (hi - lo) < bsp - 1e-9,  bsp = hi - lo;  bs = [lo hi s];  end
    end
    assert(~isempty(bs), 'f_shell:cbar', 'no tick step fits [%g %g]', cl(1), cl(2));
    tk  = round((bs(1):bs(3):bs(2))/0.01)*0.01;
    cl2 = [tk(1) tk(end)];
end

function [lim, tk] = axis_ticks_(V)
% One step s for all three axes (equal aspect): each axis spans 4s, ticks at c-s, c, c+s, c = data centre on the s grid.
    rg = cellfun(@(v) max(v(:)) - min(v(:)), V);   md = cellfun(@(v) (max(v(:)) + min(v(:)))/2, V);
    cand = [1 2 3 4 5 6 8 10 12 15 20 25 30 40 50 60 80 100 120 150 200];
    for s = cand(cand >= max(rg)/4)
        c = round(md/s)*s;  lo = c - 2*s;  hi = c + 2*s;
        ok = all(cellfun(@(v) min(v(:)), V) >= lo + 0.02*s) && all(cellfun(@(v) max(v(:)), V) <= hi - 0.02*s);
        if ok, break; end
    end
    lim = [lo(:) hi(:)];  tk = c(:) + [-s 0 s];
end
