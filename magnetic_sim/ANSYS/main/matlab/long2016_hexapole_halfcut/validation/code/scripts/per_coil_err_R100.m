%% per_coil_err_R100.m -- mean per-point error of the K_I R=100um fit, per coil/pole
clear; clc;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\hexapole-long2016\analysis');
rr = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut';
S = load('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit\fit_KI_ball\fit_KI_R100.mat');
cnst = mt_constants();
tip = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);
isL = logical(cnst.pole_is_lower);
R = 100e-6;
fprintf('coil  pole  layer   Npts   mean_err[%%]\n');
for k = 1:6
    d=import_ansys_data(fullfile(rr, sprintf('coil%d',k), 'standard'),'wp',sprintf('coil%d',k));
    air=filter_iron_nodes(d.x,d.y,d.z,cnst,struct('visualize',false));
    zwp=d.z-cnst.SPH_OFST; P=[d.x(air),d.y(air),zwp(air)]; Bn=-[d.bx(air),d.by(air),d.bz(air)];
    idx=find(sum(P.^2,2)<R^2); rng(0); if numel(idx)>S.Nmax, idx=idx(randperm(numel(idx),S.Nmax)); end
    p=P(idx,:); Bf=Bn(idx,:); pj=S.apdl_to_paper_idx(k);
    pn=p/S.ell; w=S.gB*S.Khat(:,pj)*S.I_actual; Bm=zeros(size(p));
    for i=1:6, dd=pn-dhat(:,i).'; Bm=Bm+w(i)*dd./vecnorm(dd,2,2).^3; end
    e=vecnorm(Bm-Bf,2,2)./vecnorm(Bf,2,2)*100;
    lay='Upper'; if isL(pj), lay='Lower'; end
    fprintf('  %d    P%d   %-5s  %5d   %.2f\n', k, pj, lay, numel(idx), mean(e));
end
