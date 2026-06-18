%% validate_KI.m — validate the cube100 charge fit on a 100-500 um shell
%  (out-of-fit extrapolation). Picks 100 air nodes pooled over the 6 coils in
%  100um < max(|x|,|y|,|zwp|) <= 500um, computes per-point model-vs-FEM error %,
%  and plots error % vs point (sorted by distance from WP, colored by distance).
clear; clc; close all;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');

S = load('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit\fit_KI_full.mat');
cnst = S.cnst; Khat = S.Khat; ell = S.ell; gB = S.gB; I = S.I_actual;
apdl_to_paper_idx = S.apdl_to_paper_idx;
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
fig_dir      = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);                  % 3x6 unit pole-tip directions
r_in = 100e-6;  r_out = 300e-6;              % validation shell

% ---- pool shell air nodes across 6 coils (FEM negated -> source convention) ----
P = []; B = []; pj_all = [];
for k = 1:6
    cn  = sprintf('coil%d',k);
    d   = import_ansys_data(fullfile(results_root, cn, 'standard'),'wp',cn);
    air = filter_iron_nodes(d.x,d.y,d.z,cnst,struct('visualize',false));
    zwp = d.z - cnst.SPH_OFST;
    rr  = sqrt(d.x.^2 + d.y.^2 + zwp.^2);    % radial distance from WP
    m   = air & rr>r_in & rr<=r_out;
    P   = [P; d.x(m), d.y(m), zwp(m)];                 %#ok<AGROW>
    B   = [B; -[d.bx(m), d.by(m), d.bz(m)]];           %#ok<AGROW> negate = source conv.
    pj_all = [pj_all; repmat(apdl_to_paper_idx(k), sum(m), 1)]; %#ok<AGROW>
end
fprintf('shell air nodes pooled (6 coils, %d-%d um): %d\n', r_in*1e6, r_out*1e6, size(P,1));

% ---- pick 100 random validation samples ----
rng(0);
sel  = randperm(size(P,1), min(100,size(P,1)));
Psel = P(sel,:); Bsel = B(sel,:); pjsel = pj_all(sel);
dist = vecnorm(Psel,2,2);                    % distance from WP [m]

% ---- model field + error per sample ----
err = zeros(numel(sel),1);
for s = 1:numel(sel)
    p  = Psel(s,:)/ell;                      % normalized field point
    w  = gB * Khat(:,pjsel(s)) * I;          % 6x1 charge weights
    bm = [0 0 0];
    for i = 1:6
        dv = p - dhat(:,i)';
        bm = bm + w(i) * dv / (dv*dv')^1.5;
    end
    err(s) = 100 * norm(bm - Bsel(s,:)) / norm(Bsel(s,:));
end

[dist_s, ord] = sort(dist);  err_s = err(ord);
fprintf('error over %d pts: mean=%.2f%%  median=%.2f%%  max=%.2f%%\n', ...
        numel(sel), mean(err), median(err), max(err));

% ---- figure: error% vs point (sorted by distance, colored by distance) ----
fig = figure('Color','w','Position',[100 100 920 520]);
scatter(1:numel(sel), err_s, 44, dist_s*1e6, 'filled', 'MarkerEdgeColor','k','LineWidth',0.3);
hold on; yline(mean(err),'--r','LineWidth',1.4);
text(2, mean(err), sprintf(' mean = %.1f%%', mean(err)), 'Color','r', ...
     'VerticalAlignment','bottom','FontSize',11);
colormap(turbo); cb = colorbar; cb.Label.String = 'distance from WP  [\mum]';
xlabel('validation point  (sorted by distance from WP)','FontSize',12);
ylabel('error  [%]','FontSize',12);
title(sprintf('Charge-model validation: %d points in %d-%d \\mum shell  (fit cube = \\pm100 \\mum)', ...
      numel(sel), r_in*1e6, r_out*1e6), 'FontSize',12,'FontWeight','bold');
grid on; box on; xlim([0 numel(sel)+1]); ylim([0 max(err_s)*1.08]);
out = fullfile(fig_dir, sprintf('validate_KI_cube100_shell%d_%d.png', r_in*1e6, r_out*1e6));
exportgraphics(fig, out, 'Resolution',150);
fprintf('saved %s\n', out);
