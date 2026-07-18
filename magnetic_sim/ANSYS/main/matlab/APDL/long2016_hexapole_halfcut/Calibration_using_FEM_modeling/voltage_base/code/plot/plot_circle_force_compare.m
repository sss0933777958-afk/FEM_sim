function plot_circle_force_compare()
% plot_circle_force_compare — 沿 R=50µm 圓的兩種力疊圖（coil1@1A）：
%   (綠) raw FEM  F=0.0225·∇(|B|²)（gap_calibrate coil1 raw 節點內插+差分；鋸齒＝measured-like）
%   (黑) voltage 18-param 解析模型 F=1.875e-7·(VᵀD̄ᵀ[∇(SᵀS)]D̄V)（平滑＝modeling）
%   同 pN 軸、3 subplot Fx/Fy/Fz vs 圓上 index。measure/WP frame、單位 pN/µm/mT。
%   輸出：voltage_base/figures/eighteen_param/circle_force_compare.png。

    here  = fileparts(mfilename('fullpath'));
    vbase = fileparts(fileparts(here));
    figdir = fullfile(vbase,'figures','eighteen_param'); if ~exist(figdir,'dir'); mkdir(figdir); end
    outpng = fullfile(figdir,'circle_force_compare.png');

    model = 'long2016_hexapole_halfcut';
    ML    = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut';
    addpath(fullfile(ML,'common'));
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    cnst = mt_constants();

    Nn = 720; th = linspace(0,2*pi,Nn).'; Rc = 50;
    Pm = [Rc*cos(th), Rc*sin(th), zeros(Nn,1)];        % µm, measure/WP frame
    idx = (1:Nn).';

    % ===== (綠) raw FEM 力 =====
    d   = import_ansys_data(ansys_path(model,'data','gap_calibrate','coil1'),'all','coil1');
    air = filter_iron_nodes(d.x,d.y,d.z,cnst,struct('visualize',false));
    zwp = d.z - cnst.SPH_OFST;
    P0  = [d.x(air), d.y(air), zwp(air)]*1e6;  B0 = 1e3*[d.bx(air), d.by(air), d.bz(air)];
    inb = sum(P0.^2,2) < 150^2;  P0 = P0(inb,:); B0 = B0(inb,:);
    IBx = scatteredInterpolant(P0(:,1),P0(:,2),P0(:,3),B0(:,1),'linear','none');
    IBy = scatteredInterpolant(P0(:,1),P0(:,2),P0(:,3),B0(:,2),'linear','none');
    IBz = scatteredInterpolant(P0(:,1),P0(:,2),P0(:,3),B0(:,3),'linear','none');
    B2  = @(P) IBx(P(:,1),P(:,2),P(:,3)).^2 + IBy(P(:,1),P(:,2),P(:,3)).^2 + IBz(P(:,1),P(:,2),P(:,3)).^2;
    dl1 = 12;  ex=[dl1 0 0]; ey=[0 dl1 0]; ez=[0 0 dl1];
    Rx = 0.0225*(B2(Pm+ex)-B2(Pm-ex))/(2*dl1);
    Ry = 0.0225*(B2(Pm+ey)-B2(Pm-ey))/(2*dl1);
    Rz = 0.0225*(B2(Pm+ez)-B2(Pm-ez))/(2*dl1);

    % ===== (黑) voltage 18-param 解析力 =====
    S = load(fullfile(vbase,'data','calib_D_gap_calibrate_bias.mat'));
    q = S.D_bar * S.Vmat_p(:,1);  Pc = S.Pc;  ell = S.ell_hat;  R_act = S.R_act;
    dl2 = 0.5;  gx=[dl2 0 0]; gy=[0 dl2 0]; gz=[0 0 dl2];
    Vx = 1.875e-7*(phi(Pm+gx,R_act,ell,Pc,q)-phi(Pm-gx,R_act,ell,Pc,q))/(2*dl2);
    Vy = 1.875e-7*(phi(Pm+gy,R_act,ell,Pc,q)-phi(Pm-gy,R_act,ell,Pc,q))/(2*dl2);
    Vz = 1.875e-7*(phi(Pm+gz,R_act,ell,Pc,q)-phi(Pm-gz,R_act,ell,Pc,q))/(2*dl2);

    % ===== 疊圖 =====
    f=figure('Color','w','Position',[100 100 1180 940]);
    lbl={'F_x','F_y','F_z'}; RAW={Rx,Ry,Rz}; VLT={Vx,Vy,Vz};
    for k=1:3
        subplot(3,1,k); hold on;
        hR=plot(idx,RAW{k},'-','Color',[0 0.7 0],'LineWidth',1.3);
        hV=plot(idx,VLT{k},'k-','LineWidth',1.6);
        grid on; set(gca,'FontSize',13,'FontWeight','bold','LineWidth',1.5);
        ylabel([lbl{k} ' (pN)'],'FontWeight','bold'); xlim([1 Nn]);
        if k==1
            title('coil1 @ 1A    force on R=50 \mum circle:  raw FEM  vs  voltage 18-param model','FontWeight','bold');
            legend([hR hV],{'raw FEM (0.0225\cdot\nabla|B|^2)','voltage 18-param model'},'Location','best','FontSize',11);
        end
        if k==3; xlabel('point index (along R=50 \mum circle)','FontWeight','bold'); end
    end
    exportgraphics(f,outpng,'Resolution',150); close(f);
    fprintf('wrote %s\n', outpng);
end

function v = phi(Pmat, R_act, ell, Pc, q)
    pb = (R_act*Pmat.')/ell;  N=size(Pmat,1); Btil=zeros(3,N);
    for i=1:6, Di=pb-Pc(:,i); r3=sum(Di.^2,1).^1.5; Btil=Btil+q(i)*(Di./r3); end
    v = sum(Btil.^2,1).';
end
