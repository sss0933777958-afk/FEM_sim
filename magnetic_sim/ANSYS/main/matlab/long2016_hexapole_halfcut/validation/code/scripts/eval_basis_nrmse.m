%% eval_basis_nrmse.m — per-coil (basis) NRMSE in the R<=50 um ball vs R_param.
%  Linearity: any current combo's error = sum_k c_k * (coil-k error), and single-coil
%  NRMSE is current-magnitude invariant. So the 6 single-coil NRMSE curves (+ worst
%  envelope) fully characterize model fidelity for ANY current -> decides calibration R.
%  NRMSE_max = sqrt(mean||B_model-B_FEM||^2)/max||B_FEM|| * 100  (nrmse.pdf, max-field).
clear; clc; close all;

addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
data_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit\fitting_trend';
fig_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';

cnst = mt_constants();
apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];          % apdl coil k -> paper pole
paper_name = {'P1','P2','P3','P4','P5','P6'};
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);

%% load 6 coils (air, WP frame, negated = source convention)
fprintf('loading 6 coils ...\n');
P1 = []; Bn = cell(1,6);
for k = 1:6
    cname = sprintf('coil%d', k);
    d   = import_ansys_data(fullfile(results_root, cname), 'wp', cname);
    air = filter_iron_nodes(d.x, d.y, d.z, cnst, struct('visualize', false));
    zwp = d.z - cnst.SPH_OFST;
    Pk  = [d.x(air), d.y(air), zwp(air)];
    if k==1, P1 = Pk; end
    Bn{k} = -[d.bx(air), d.by(air), d.bz(air)];
end

%% test set = all air nodes in R<=50 um ball (positions identical across coils)
R_test = 50e-6;
idx = find(sum(P1.^2,2) <= R_test^2);
Ptest = P1(idx,:);  N = numel(idx);
fprintf('test set = %d nodes in R<=50 um ball\n', N);

%% load dense fit params (I_actual=1)
S = load(fullfile(data_dir,'KI_trend_sweep_dense.mat'));
nR = numel(S.R_um);

%% per-coil NRMSE vs R_param
NRMSE = nan(6, nR);
for kc = 1:6                                   % apdl coil index
    pj = apdl_to_paper_idx(kc);                % paper column
    Bfem = Bn{kc}(idx,:);                       % single-coil FEM field (1A)
    maxFEM = max(vecnorm(Bfem,2,2));
    Ivec = zeros(6,1); Ivec(pj) = 1;            % only this coil on (1A)
    for ri = 1:nR
        Khat = S.Khat_s(:,:,ri);
        if any(isnan(Khat(:))), continue; end
        ell = S.ell_s(ri); gB = S.gB_s(ri);
        w = gB * Khat * Ivec;                   % = gB*Khat(:,pj)
        PN = Ptest / ell;
        Bmod = zeros(N,3);
        for i = 1:6
            D  = PN - dhat(:,i).';
            r3 = (sum(D.^2,2)).^1.5;
            Bmod = Bmod + w(i) * (D ./ r3);
        end
        NRMSE(kc,ri) = sqrt(mean(sum((Bmod - Bfem).^2,2))) / maxFEM * 100;
    end
end

R = S.R_um;
Eenv = max(NRMSE,[],1);          % worst-case envelope
Amean = mean(NRMSE,1);           % mean across coils
[Emin, imin] = min(Eenv);
fprintf('worst-envelope min = %.3f%% at R=%d um\n', Emin, R(imin));

%% figure
cols = lines(6);
f = figure('Color','w','Position',[100 100 860 540]); hold on;
hc = gobjects(1,6);
for kc = 1:6
    hc(kc) = plot(R, NRMSE(kc,:), '-o', 'LineWidth',1.4, 'MarkerSize',3, ...
                  'Color',cols(kc,:), 'MarkerFaceColor',cols(kc,:));
end
grid on; box on; set(gca,'FontSize',11); xlim([0 500]);
yline(2,'--','Color',[0.6 0.2 0.6],'LineWidth',1.4,'Label','threshold 2%','FontSize',10, ...
      'LabelHorizontalAlignment','left');
xline(350,':','Color',[0.6 0.2 0.6],'LineWidth',1.2,'Label','R_{hi}=350','FontSize',9);
xlabel('Sampling radius  $R$  [$\mu$m]  (params from)','Interpreter','latex','FontSize',13);
ylabel('per-coil $\mathrm{NRMSE}_{\max}$  [\%]','Interpreter','latex','FontSize',13);
title(sprintf('Per-coil basis NRMSE on $R\\le50\\,\\mu$m ball ($N=%d$)', N), ...
      'Interpreter','latex','FontSize',12);
legnames = arrayfun(@(kc) sprintf('coil%d (%s)', kc, paper_name{apdl_to_paper_idx(kc)}), 1:6, 'uni',0);
legend(hc, legnames, 'Location','northwest','FontSize',10, 'Interpreter','none');
exportgraphics(f, fullfile(fig_dir,'basis_nrmse_6lines.png'), 'Resolution',150);

save(fullfile(data_dir,'combo_basis_nrmse.mat'), 'R','NRMSE','Eenv','Amean','N','R_test','apdl_to_paper_idx');
fprintf('wrote basis_nrmse_vs_R.png + combo_basis_nrmse.mat\n');
