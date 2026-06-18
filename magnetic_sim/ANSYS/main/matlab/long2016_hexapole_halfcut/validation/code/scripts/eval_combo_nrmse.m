%% eval_combo_nrmse.m — NRMSE(max) of charge-model vs FEM, per current combination,
%  as a function of the fit's sampling radius R_param.
%  Test set: ALL air nodes in the R<=50 um ball (348 pts, shared positions across coils).
%  FEM combo field = sum_k c_k * Bn_k (superposition; coils at 1A FEM unit, negated source).
%  Model: params {ell,gB,Khat} from dense sweep (I_actual=1) -> currents plug in directly.
%  NRMSE_max = sqrt(mean ||B_model-B_FEM||^2) / max||B_FEM|| * 100   (nrmse.pdf, max-field).
%  5 figures: combo_nrmse_C1..C5.png (NRMSE vs R_param).
clear; clc; close all;

addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
data_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit\fitting_trend';
fig_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';

cnst = mt_constants();
apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);

%% load 6 coils (air, WP frame, negated = source convention)
fprintf('loading 6 coils ...\n');
P = cell(1,6); Bn = cell(1,6);
for k = 1:6
    cname = sprintf('coil%d', k);
    d   = import_ansys_data(fullfile(results_root, cname), 'wp', cname);
    air = filter_iron_nodes(d.x, d.y, d.z, cnst, struct('visualize', false));
    zwp = d.z - cnst.SPH_OFST;
    P{k}  = [d.x(air), d.y(air), zwp(air)];
    Bn{k} = -[d.bx(air), d.by(air), d.bz(air)];
end

%% test set = all air nodes in R<=50 um ball (positions identical across coils)
R_test = 50e-6;
idx = find(sum(P{1}.^2,2) <= R_test^2);
Ptest = P{1}(idx,:);
N = numel(idx);
Bk = zeros(N,3,6);                 % per-coil field at test nodes (aligned by index)
for k = 1:6, Bk(:,:,k) = Bn{k}(idx,:); end
fprintf('test set = %d nodes in R<=50 um ball\n', N);

%% load dense fit params (I_actual=1)
S = load(fullfile(data_dir,'KI_trend_sweep_dense.mat'));   % R_um, ell_s[m], gB_s, Khat_s
nR = numel(S.R_um);

%% 5 current combinations (APDL coil order, Amps)
combos = [1 1 1 1 1 1;
          2 2 2 2 2 2;
          2 2 2 1 1 1;
          1 1 1 2 2 2;
          2 1 2 1 2 1];
labels = {'C1','C2','C3','C4','C5'};
blue = [0.2 0.45 0.85];

NRMSE = nan(size(combos,1), nR);
for c = 1:size(combos,1)
    Cc = combos(c,:);
    % FEM superposed field for this combo
    Bfem = zeros(N,3);
    for k = 1:6, Bfem = Bfem + Cc(k)*Bk(:,:,k); end
    maxFEM = max(vecnorm(Bfem,2,2));
    % current vector in paper-pole order
    Ivec = zeros(6,1);
    for k = 1:6, Ivec(apdl_to_paper_idx(k)) = Cc(k); end

    for ri = 1:nR
        Khat = S.Khat_s(:,:,ri);
        if any(isnan(Khat(:))), continue; end
        ell = S.ell_s(ri); gB = S.gB_s(ri);
        w = gB * Khat * Ivec;                 % I_actual=1 -> currents direct
        PN = Ptest / ell;
        Bmod = zeros(N,3);
        for i = 1:6
            D  = PN - dhat(:,i).';
            r3 = (sum(D.^2,2)).^1.5;
            Bmod = Bmod + w(i) * (D ./ r3);
        end
        rmse = sqrt(mean(sum((Bmod - Bfem).^2,2)));
        NRMSE(c,ri) = rmse / maxFEM * 100;
    end
    fprintf('%s  curr=%s  maxFEM=%.3e T  min NRMSE=%.2f%% @ R=%d um\n', ...
        labels{c}, mat2str(Cc), maxFEM, min(NRMSE(c,:)), S.R_um(find(NRMSE(c,:)==min(NRMSE(c,:)),1)));
end

%% one figure per combination
R = S.R_um;
for c = 1:size(combos,1)
    f = figure('Color','w','Position',[100 100 760 480]);
    plot(R, NRMSE(c,:), '-o','LineWidth',1.4,'MarkerSize',3.5, ...
         'MarkerFaceColor',blue,'Color',blue);
    grid on; box on; set(gca,'FontSize',11); xlim([0 500]);
    xlabel('Sampling radius  $R$  [$\mu$m]  (params from)','Interpreter','latex','FontSize',13);
    ylabel('$\mathrm{NRMSE}_{\max}$  [\%]','Interpreter','latex','FontSize',13);
    title(sprintf('NRMSE on $R\\le50\\,\\mu$m ball, %s: $I=[%s]$ A  ($N=%d$)', ...
          labels{c}, strjoin(string(combos(c,:)),','), N), 'Interpreter','latex','FontSize',12);
    exportgraphics(f, fullfile(fig_dir, sprintf('combo_nrmse_%s.png',labels{c})), 'Resolution',150);
end
save(fullfile(data_dir,'combo_nrmse.mat'), 'R','NRMSE','combos','labels','N','R_test');
fprintf('\nwrote combo_nrmse_C1..C5.png + combo_nrmse.mat\n');
