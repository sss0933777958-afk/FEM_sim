%% plot_disc_sampling.m — visualise the B_surface disc sampling pattern
%
%  Shows the 37-point concentric-ring grid used by compute_B_bar_matrix.m
%  to area-average B.n over the Ø0.3 mm Hall sensor disc.
%
%  Output: magnetic_sim/ANSYS/main/figures/long2016_hexapole_halfcut/Bsurf_disc_sampling.png

clear; clc; close all;

%% --- disc sample points (matches compute_B_bar_matrix.m) ---
%  centre + 15 rings at 0.01..0.15 mm (radial spacing 0.01 mm);
%  ring k has 6k points (proportional to radius). Total 1 + 6*120 = 721.
sensor_radius = 0.15;            % mm  (Ø0.3 mm)
n_rings = 15;
ring_dr = sensor_radius / n_rings;   % 0.01 mm
disc_local = [0, 0];
ring_id = 0;                      % 0 = centre
for k = 1:n_rings
    r_k = ring_dr * k;
    n_k = 6*k;
    phi = (0:n_k-1) * 2*pi / n_k;
    disc_local = [disc_local; r_k*cos(phi(:)), r_k*sin(phi(:))];
    ring_id = [ring_id; k*ones(n_k,1)];
end
N_disc_points = size(disc_local,1);
assert(N_disc_points == 721, 'point count mismatch');

%% --- figure ---
fig = figure('Color','w','Position',[120 120 720 720]);
ax = axes(fig); hold(ax,'on');

th = linspace(0,2*pi,200);
% boundary disc (r = 0.15 mm)
fill(sensor_radius*cos(th), sensor_radius*sin(th), [0.93 0.95 1.0], ...
     'EdgeColor',[0.20 0.35 0.65],'LineWidth',2.0);
% ring guide circles (dashed)
for k = 1:n_rings
    r_k = ring_dr*k;
    plot(ax, r_k*cos(th), r_k*sin(th), '--', 'Color',[0.70 0.70 0.75],'LineWidth',0.5);
end

% sample points, coloured by ring index (turbo colormap)
cmap = turbo(n_rings+1);
for g = 0:n_rings
    m = ring_id == g;
    plot(ax, disc_local(m,1), disc_local(m,2), 'o', ...
        'MarkerSize',4.5,'MarkerFaceColor',cmap(g+1,:), ...
        'MarkerEdgeColor','none');
end
colormap(ax, cmap); caxis(ax,[0 sensor_radius]);
cb = colorbar(ax); cb.Label.String = 'ring radius  [mm]';

axis(ax,'equal'); box(ax,'on'); grid(ax,'on');
lim = sensor_radius*1.25;
xlim(ax,[-lim lim]); ylim(ax,[-lim lim]);
set(ax,'FontSize',11);
xlabel(ax,'sensor-local u  [mm]','FontSize',12);
ylabel(ax,'sensor-local v  [mm]','FontSize',12);
title(ax,{sprintf('B_{surface} disc sampling — %d-point concentric-ring grid', N_disc_points), ...
          ['\rm\fontsize{10}centre + 15 rings (0.01..0.15 mm, \Deltar = 0.01 mm); ' ...
           'ring k = 6k pts; plain mean of B\cdotn \approx area average']}, ...
      'FontSize',13);

out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
if ~exist(out_dir,'dir'); mkdir(out_dir); end
out_path = fullfile(out_dir,'Bsurf_disc_sampling.png');
exportgraphics(fig, out_path, 'Resolution',300);
fprintf('points: %d (centre 1 + 15 rings, ring k = 6k pts)\n', N_disc_points);
fprintf('Saved: %s\n', out_path);
