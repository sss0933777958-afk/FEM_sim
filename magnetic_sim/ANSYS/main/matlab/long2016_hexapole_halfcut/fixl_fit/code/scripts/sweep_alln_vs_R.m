%% sweep_alln_vs_R.m — deterministic calibration sweep: at each radius R use ALL
%  in-ball FEM nodes per coil (no subsampling, no seed, no redraw), one lsqnonlin fit.
%  Result at each R is unique -> reproducible by construction (re-running gives the
%  identical numbers). Removes the "how many redraws" question entirely.
%  Outputs param trends {ell,gB,||K||}(R) and accuracy NRMSE(R) on the 50um ball.
clear; clc; close all;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
data_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit\fitting_trend';
fig_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';

cnst = mt_constants();
apdl_to_paper_idx = [1,3,6,5,2,4];
I_actual = 1;
tip = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);

fprintf('loading 6 coils ...\n');
C = struct('P',{},'Bn',{},'pj',{});
for k=1:6
    cn=sprintf('coil%d',k);
    d=import_ansys_data(fullfile(results_root, cn, 'standard'),'wp',cn);
    air=filter_iron_nodes(d.x,d.y,d.z,cnst,struct('visualize',false));
    zwp=d.z-cnst.SPH_OFST;
    C(k).P=[d.x(air),d.y(air),zwp(air)]; C(k).Bn=-[d.bx(air),d.by(air),d.bz(air)];
    C(k).pj=apdl_to_paper_idx(k);
end

R_um = 50:5:500;

%% R<=50um operating-ball test set (for accuracy of each fit)
idx50 = find(sum(C(1).P.^2,2) <= (50e-6)^2);
Ptest50 = C(1).P(idx50,:); N50 = numel(idx50);
Bk50 = zeros(N50,3,6); for k=1:6, Bk50(:,:,k)=C(k).Bn(idx50,:); end
fprintf('operating ball N=%d nodes\n', N50);

Khat0 = eye(6)-ones(6)/6; ell0=0.5e-3; gB0=1;
freemask=true(6); freemask(1,1)=false;
opts=optimoptions('lsqnonlin','Display','off','MaxFunctionEvaluations',1e5, ...
    'MaxIterations',4e3,'FunctionTolerance',1e-20,'StepTolerance',1e-12);

nR=numel(R_um);
ell_R=nan(1,nR); gB_R=nan(1,nR); Kf_R=nan(1,nR); npts=nan(1,nR); NRMSE_R=nan(1,nR);
Ksave=zeros(6,6,nR);
fprintf('\n  R[um]  npts/coil(min)  ell[mm]  gB[1e-3]  ||K||_F  NRMSE50%%\n');
for ri=1:nR
    R=R_um(ri)*1e-6;
    coil=struct('p',{},'bfem',{},'pj',{});
    nmin=inf;
    for k=1:6
        ii=find(sum(C(k).P.^2,2)<R^2); nmin=min(nmin,numel(ii));
        Bn=C(k).Bn(ii,:);
        coil(k).p=C(k).P(ii,:); coil(k).bfem=[Bn(:,1);Bn(:,2);Bn(:,3)]; coil(k).pj=C(k).pj;
    end
    npts(ri)=nmin;
    xf=lsqnonlin(@(x) resid_all(x,coil,dhat,I_actual,freemask),[ell0*1e3;gB0;Khat0(freemask)],[],[],opts);
    [ell,gB,Khat]=unpack(xf,freemask);
    ell_R(ri)=ell; gB_R(ri)=gB; Kf_R(ri)=norm(Khat,'fro'); Ksave(:,:,ri)=Khat;
    % accuracy on 50um ball (worst coil)
    PN50=Ptest50/ell; ev=zeros(1,6);
    for kc=1:6
        Iv=zeros(6,1); Iv(apdl_to_paper_idx(kc))=1; w=gB*Khat*Iv; Bm=zeros(N50,3);
        for i=1:6, Dm=PN50-dhat(:,i).'; r3=(sum(Dm.^2,2)).^1.5; Bm=Bm+w(i)*(Dm./r3); end
        Bf=Bk50(:,:,kc); ev(kc)=sqrt(mean(sum((Bm-Bf).^2,2)))/max(vecnorm(Bf,2,2))*100;
    end
    NRMSE_R(ri)=max(ev);
    fprintf('  %4d   %8d        %.4f   %6.3f    %.4f   %6.3f\n', ...
        R_um(ri), nmin, ell*1e3, gB*1e3, Kf_R(ri), NRMSE_R(ri));
end

[bestN,bi]=min(NRMSE_R); R_best=R_um(bi);
fprintf('\n best accuracy: R*=%d um, NRMSE=%.3f%% (ell=%.4f mm, gB=%.4e, ||K||=%.4f)\n', ...
    R_best, bestN, ell_R(bi)*1e3, gB_R(bi), Kf_R(bi));

%% figure 1: parameter trends vs R (deterministic)
f=figure('Color','w','Position',[100 100 880 560]);
tl=tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
nexttile; plot(R_um,ell_R*1e3,'-o','LineWidth',1.4,'MarkerSize',3,'Color',[0.2 0.45 0.85]);
grid on; box on; ylabel('$\hat\ell$ [mm]','Interpreter','latex','FontSize',12); set(gca,'FontSize',10);
title('Deterministic fit (all in-ball nodes) vs sampling radius $R$','Interpreter','latex','FontSize',12);
xline(R_best,'-.','Color',[0.1 0.55 0.2],'LineWidth',1.3);
nexttile; plot(R_um,gB_R*1e3,'-s','LineWidth',1.4,'MarkerSize',3,'Color',[0.85 0.4 0.1]);
grid on; box on; ylabel('$\hat g_B$ [$\times10^{-3}$]','Interpreter','latex','FontSize',12); set(gca,'FontSize',10);
xline(R_best,'-.','Color',[0.1 0.55 0.2],'LineWidth',1.3);
nexttile; plot(R_um,Kf_R,'-^','LineWidth',1.4,'MarkerSize',3,'Color',[0.45 0.3 0.7]);
grid on; box on; ylabel('$\|\hat K\|_F$','Interpreter','latex','FontSize',12); set(gca,'FontSize',10);
xlabel('sampling radius $R$ [$\mu$m]','Interpreter','latex','FontSize',12);
xline(R_best,'-.','Color',[0.1 0.55 0.2],'LineWidth',1.3,'Label',sprintf('R^*=%d',R_best),'FontSize',9);
exportgraphics(f, fullfile(fig_dir,'sweep_alln_params_vs_R.png'),'Resolution',150);

%% figure 2: accuracy NRMSE(R) on 50um ball
f2=figure('Color','w','Position',[100 100 860 500]); hold on;
plot(R_um,NRMSE_R,'-o','LineWidth',1.6,'MarkerSize',3.5,'Color',[0.2 0.45 0.85],'MarkerFaceColor',[0.2 0.45 0.85]);
plot(R_best,bestN,'p','MarkerSize',14,'MarkerFaceColor',[0.1 0.55 0.2],'MarkerEdgeColor','k');
yline(2,'--','Color',[0.6 0.2 0.6],'LineWidth',1.1,'Label','2\%','Interpreter','latex','FontSize',9);
set(gca,'YScale','log'); grid on; box on; set(gca,'FontSize',11); xlim([40 510]);
xlabel('sampling radius $R$ [$\mu$m]','Interpreter','latex','FontSize',13);
ylabel('worst-coil $\mathrm{NRMSE}_{\max}$ on 50\,$\mu$m ball [\%]','Interpreter','latex','FontSize',12);
title(sprintf('Accuracy of deterministic all-node fit vs $R$  (best: $R^*{=}%d\\,\\mu$m, %.2f\\%%)',R_best,bestN), ...
      'Interpreter','latex','FontSize',12);
exportgraphics(f2, fullfile(fig_dir,'sweep_alln_nrmse_vs_R.png'),'Resolution',150);

save(fullfile(data_dir,'sweep_alln_vs_R.mat'),'R_um','ell_R','gB_R','Kf_R','Ksave','npts','NRMSE_R','R_best');
fprintf('\nwrote sweep_alln_params_vs_R.png + sweep_alln_nrmse_vs_R.png + sweep_alln_vs_R.mat\n');

%% local
function r = resid_all(x, coil, dhat, I, freemask)
    [ell, gB, Khat] = unpack(x, freemask); r=[];
    for k=1:numel(coil)
        pn=coil(k).p/ell; N=size(pn,1); B=zeros(3*N,1); w=gB*Khat(:,coil(k).pj)*I;
        for i=1:6
            dx=pn(:,1)-dhat(1,i); dy=pn(:,2)-dhat(2,i); dz=pn(:,3)-dhat(3,i);
            r3=(dx.^2+dy.^2+dz.^2).^1.5; B=B+w(i)*[dx./r3; dy./r3; dz./r3];
        end
        r=[r; B-coil(k).bfem]; %#ok<AGROW>
    end
end
function [ell, gB, Khat] = unpack(x, freemask)
    ell=x(1)*1e-3; gB=x(2); Khat=zeros(6); Khat(1,1)=5/6; Khat(freemask)=x(3:end);
end
