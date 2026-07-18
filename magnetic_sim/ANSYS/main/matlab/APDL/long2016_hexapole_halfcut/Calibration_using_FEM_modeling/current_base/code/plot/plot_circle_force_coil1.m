function plot_circle_force_coil1()
% plot_circle_force_coil1 — 沿 R=50µm 圓(measure/WP frame, z=0)算磁珠力 F=0.0225·∇(|B|²)
%   B = long2016 gap_calibrate coil1 @ 1A（raw FEM 節點內插）。單位：場 mT、長度 µm、力 pN（Unit Reference Sheet）。
%   |B|² 對 B 符號/frame 旋轉不變 → 用 measure frame raw B；∇(|B|²) 中央差分 → mT²/µm；×0.0225 → pN。
%   3 subplot：Fx/Fy/Fz (pN) vs 圓上取樣點 index。輸出 figures/shared/circle_force_coil1.png。

    model = 'long2016_hexapole_halfcut';
    ML    = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut';
    addpath(fullfile(ML,'common'));                                                   % ansys_path
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis'); % import_ansys_data/filter_iron_nodes/mt_constants
    here   = fileparts(mfilename('fullpath'));
    cbase  = fileparts(fileparts(here));
    figdir = fullfile(cbase,'figures','shared'); if ~exist(figdir,'dir'); mkdir(figdir); end
    outpng = fullfile(figdir,'circle_force_coil1.png');

    % ---- 載入 coil1 gap_calibrate（measure/WP frame, µm, mT）----
    cnst = mt_constants();
    d    = import_ansys_data(ansys_path(model,'data','gap_calibrate','coil1'),'all','coil1');
    air  = filter_iron_nodes(d.x,d.y,d.z,cnst,struct('visualize',false));
    zwp  = d.z - cnst.SPH_OFST;                                   % WP frame z
    Pm   = [d.x(air), d.y(air), zwp(air)]*1e6;                    % µm, WP frame（球殼中心=WP=原點）
    Bm   = 1e3*[d.bx(air), d.by(air), d.bz(air)];                % mT（Tesla→mT；|B|² 對符號不變）
    r2   = sum(Pm.^2,2); inb = r2 < 150^2;                       % R<=150µm 球殼
    Pm = Pm(inb,:); Bm = Bm(inb,:);
    fprintf('R<=150µm 球內節點 = %d\n', size(Pm,1));

    % ---- B 分量內插器（scatteredInterpolant, µm→mT）----
    IBx = scatteredInterpolant(Pm(:,1),Pm(:,2),Pm(:,3),Bm(:,1),'linear','none');
    IBy = scatteredInterpolant(Pm(:,1),Pm(:,2),Pm(:,3),Bm(:,2),'linear','none');
    IBz = scatteredInterpolant(Pm(:,1),Pm(:,2),Pm(:,3),Bm(:,3),'linear','none');
    B2  = @(x,y,z) IBx(x,y,z).^2 + IBy(x,y,z).^2 + IBz(x,y,z).^2; % |B|² [mT²]

    % ---- R=50µm 圓（z=0）取樣 ----
    Nn = 720; th = linspace(0,2*pi,Nn).';
    Rc = 50; cx = Rc*cos(th); cy = Rc*sin(th); cz = zeros(Nn,1);

    % ---- ∇(|B|²) 中央差分（δ µm）→ F=0.0225·∇(|B|²) [pN] ----
    dl = 12;  PRE = 0.0225;
    gx = (B2(cx+dl,cy,cz) - B2(cx-dl,cy,cz))/(2*dl);
    gy = (B2(cx,cy+dl,cz) - B2(cx,cy-dl,cz))/(2*dl);
    gz = (B2(cx,cy,cz+dl) - B2(cx,cy,cz-dl))/(2*dl);
    Fx = PRE*gx; Fy = PRE*gy; Fz = PRE*gz;                       % pN
    idx = (1:Nn).';
    fprintf('|F| range: Fx[%.3g,%.3g] Fy[%.3g,%.3g] Fz[%.3g,%.3g] pN\n', ...
        min(Fx),max(Fx),min(Fy),max(Fy),min(Fz),max(Fz));

    % ---- 3 subplot ----
    f = figure('Color','w','Position',[100 100 1150 900]);
    lbl = {'F_x','F_y','F_z'}; Fc = {Fx,Fy,Fz};
    for k = 1:3
        subplot(3,1,k); plot(idx,Fc{k},'k-','LineWidth',1.4); grid on;
        set(gca,'FontSize',14,'FontWeight','bold','LineWidth',1.5);
        ylabel([lbl{k} ' (pN)'],'FontWeight','bold'); xlim([1 Nn]);
        if k==1; title('coil1 @ 1A  force on R=50\mum xy circle   F = 0.0225\cdot\nabla(|B|^2)','FontWeight','bold'); end
        if k==3; xlabel('point index (along R=50\mum circle)','FontWeight','bold'); end
    end
    exportgraphics(f,outpng,'Resolution',150); close(f);
    fprintf('wrote %s\n', outpng);
end
