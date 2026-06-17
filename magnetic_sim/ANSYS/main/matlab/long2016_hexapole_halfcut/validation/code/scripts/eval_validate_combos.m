%% eval_validate_combos.m — closed-loop check: 10 random current combos @ R*=170 params,
%  NRMSE on the R<=50 um ball. Varies BOTH ratio and overall magnitude to show NRMSE
%  does not amplify with current size. Uses AS-FITTED Khat (raw gauge) with raw currents.
clear; clc; close all;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut';
data_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit\fitting_trend';
out_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit';
fig_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';

cnst = mt_constants();
apdl_to_paper_idx = [1,3,6,5,2,4];
tip = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);

%% R*=170 calibration (AS-FITTED gauge, for prediction)
% frozen final calibration (10-fit mean @ R*=250); use as-fitted Khat_raw for prediction
F = load('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit\calibration_final.mat');
R_star = F.R_star; ell = F.ell_hat; gB = F.gB_hat; Khat = F.Khat_raw;
fprintf('R*=%d (frozen 10-fit mean): ell=%.4f mm, gB=%.4e, ||K||_F=%.3f\n', R_star, ell*1e3, gB, norm(Khat,'fro'));

%% 6 coils + 50 um ball test set
P1=[]; Bn=cell(1,6);
for k=1:6
    cn=sprintf('coil%d',k);
    d=import_ansys_data(fullfile(results_root, cn, 'standard'),'wp',cn);
    air=filter_iron_nodes(d.x,d.y,d.z,cnst,struct('visualize',false));
    zwp=d.z-cnst.SPH_OFST; Pk=[d.x(air),d.y(air),zwp(air)];
    if k==1, P1=Pk; end
    Bn{k}=-[d.bx(air),d.by(air),d.bz(air)];
end
idx=find(sum(P1.^2,2)<=(50e-6)^2); Ptest=P1(idx,:); N=numel(idx);
Bk=zeros(N,3,6); for k=1:6, Bk(:,:,k)=Bn{k}(idx,:); end
fprintf('test set N=%d (R<=50um ball)\n', N);

%% 10 random current combos: random ratio x random overall scale (magnitude spread)
rng(7);
G=100;
combos = zeros(G,6); NRMSE=zeros(1,G); normI=zeros(1,G);
PN = Ptest/ell;
for gI=1:G
    c = 0.2 + 1.8*rand(1,6);                   % U(0.2,2) A per coil (realistic operating range)
    combos(gI,:)=c; normI(gI)=norm(c);
    % FEM superposition
    Bfem=zeros(N,3); for k=1:6, Bfem=Bfem+c(k)*Bk(:,:,k); end
    maxFEM=max(vecnorm(Bfem,2,2));
    % model
    Iv=zeros(6,1); for k=1:6, Iv(apdl_to_paper_idx(k))=c(k); end
    w=gB*Khat*Iv; Bmod=zeros(N,3);
    for ii=1:6, D=PN-dhat(:,ii).'; r3=(sum(D.^2,2)).^1.5; Bmod=Bmod+w(ii)*(D./r3); end
    NRMSE(gI)=sqrt(mean(sum((Bmod-Bfem).^2,2)))/maxFEM*100;
end
fprintf('\n grp   ||I||[A]   NRMSE%%\n');
for gI=1:G, fprintf('  %2d   %7.2f   %.3f\n', gI, normI(gI), NRMSE(gI)); end
fprintf('  NRMSE range %.3f - %.3f%% (mean %.3f%%)\n', min(NRMSE),max(NRMSE),mean(NRMSE));

%% figure
f=figure('Color','w','Position',[100 100 820 500]);
plot(1:G, NRMSE, '-o','LineWidth',1.0,'MarkerSize',3, ...
     'MarkerFaceColor',[0.2 0.45 0.85],'Color',[0.2 0.45 0.85]); hold on;
yline(mean(NRMSE),'--','Color',[0.5 0.5 0.5],'LineWidth',1.0, ...
      'Label',sprintf('mean %.2f%%',mean(NRMSE)),'FontSize',9);
grid on; box on; set(gca,'FontSize',11); ylim([0 2]); xlim([0 G+1]);
xlabel('current combination \#','Interpreter','latex','FontSize',13);
ylabel('$\mathrm{NRMSE}_{\max}$ on $R\le50\,\mu$m ball  [\%]','Interpreter','latex','FontSize',13);
title(sprintf('Closed-loop check: %d random currents @ $R^*{=}%d\\,\\mu$m calibration (50\\,$\\mu$m ball)', G, R_star), ...
      'Interpreter','latex','FontSize',12);
exportgraphics(f, fullfile(fig_dir,sprintf('validate_combos_nrmse_R%d.png',R_star)),'Resolution',150);
save(fullfile(out_dir,sprintf('validate_combos_R%d.mat',R_star)),'combos','normI','NRMSE','N','R_star');
fprintf('\nwrote validate_combos_nrmse.png + validate_combos.mat\n');
