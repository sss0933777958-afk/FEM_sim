%% error_per_point_R100.m
%  Per-point fitting error of the K_I R=100um ball fit, using the SAME error
%  definition as fit_J_cube40.tex:
%     dB_i = B_model(p_i) - B_FEM(p_i),   err_i = |dB_i|/|B_FEM_i| * 100 [%]
%     overall Error = mean_i(err_i)
%  Also outputs a "Percent fitting error vs index of points" star plot.
clear; clc; close all;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut';
matf = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit\fit_KI_ball\fit_KI_R100.mat';
fig_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';

S = load(matf);                      % Khat, ell, gB, I_actual, apdl_to_paper_idx, Nmax
cnst = mt_constants();
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);
R = 100e-6;  Nmax = S.Nmax;

%% rebuild the R=100 ball sample set + per-point error (pooled over 6 coils)
err = [];
for k = 1:6
    cname = sprintf('coil%d', k);
    d   = import_ansys_data(fullfile(results_root, cname), 'wp', cname);
    air = filter_iron_nodes(d.x, d.y, d.z, cnst, struct('visualize', false));
    zwp = d.z - cnst.SPH_OFST;
    P   = [d.x(air), d.y(air), zwp(air)];
    Bn  = -[d.bx(air), d.by(air), d.bz(air)];        % source convention (= fit)
    idx = find(sum(P.^2,2) < R^2);
    rng(0); if numel(idx) > Nmax, idx = idx(randperm(numel(idx), Nmax)); end
    p = P(idx,:);  Bfem = Bn(idx,:);
    pn = p / S.ell;  w = S.gB * S.Khat(:, S.apdl_to_paper_idx(k)) * S.I_actual;
    Bmod = zeros(size(p));
    for i = 1:6
        dd = pn - dhat(:,i).';   r3 = vecnorm(dd,2,2).^3;
        Bmod = Bmod + w(i) * dd ./ r3;
    end
    err = [err; vecnorm(Bmod-Bfem,2,2) ./ vecnorm(Bfem,2,2) * 100]; %#ok<AGROW>
end

overall = mean(err);
fprintf('\nR=100um  K_I fit  |  %d points\n', numel(err));
fprintf('Overall Error (mean |dB|/|B_FEM|) = %.2f%%   (median %.2f%%, max %.2f%%)\n', ...
        overall, median(err), max(err));

%% star plot (evenly-spaced subsample so it reads like the reference panel)
Nsub = min(240, numel(err));
sidx = round(linspace(1, numel(err), Nsub));
figure('Color','w','Position',[100 100 560 460]);
plot(1:Nsub, err(sidx), '-*', 'Color','b', 'MarkerSize',5, 'LineWidth',0.6); hold on;
yline(overall, 'r--', sprintf('mean = %.2f%%', overall), 'LineWidth',1.2, ...
      'LabelHorizontalAlignment','left');
xlabel('index of points'); ylabel('error [%]');
title('Percent fitting error'); grid on; box on; set(gca,'FontSize',12);
xlim([0 Nsub]); ylim([0 max(err(sidx))*1.1]);

out = fullfile(fig_dir,'fit_error_per_point_R100.png');
exportgraphics(gcf, out, 'Resolution',200);
fprintf('saved %s\n', out);
