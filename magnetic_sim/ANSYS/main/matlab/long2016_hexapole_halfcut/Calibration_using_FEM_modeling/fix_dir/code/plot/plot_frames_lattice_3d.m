function plot_frames_lattice_3d()
% PLOT_FRAMES_LATTICE_3D  measure 與 actuator 座標系 + ℓ̂ 球殼 + 6 極點 幾何示意（fix-ℓ）。
%   WP 原點同時畫兩組三軸：
%     measure  [x_m y_m z_m] = 全域 x,y,z（灰）
%     actuator [x_a y_a z_a] = d̂_P1 / d̂_P3 / d̂_P5（藍，指向 P1/P3/P5）
%   半徑 = ℓ̂（≈867µm，fit_fixl 值）的半透明球殼；殼上 6 點 = 等效磁荷 ℓ̂·d̂_k，標 P1..P6。
%   兩層「同平面」三角面：下極 P1/P3/P6（紅）、上極 P2/P4/P5（藍）。
%   出兩張圖：
%     frames_lattice_3d.png       — 只 ℓ̂ 殼
%     frames_lattice_R150_3d.png  — 多包一顆 R=150µm 校正取樣球殼（綠）
%   單位 µm；框 = box on + daspect（cube 幾何、view-robust）。

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % mt_constants
    cnst = mt_constants();

    here   = fileparts(mfilename('fullpath'));
    fixdir = fileparts(fileparts(here));
    matf   = fullfile(fixdir, 'data', 'fit_fixl_R150um_gap200um_mueq.mat');
    figdir = fullfile(fixdir, 'figures');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    S = load(matf, 'ell');   ell = S.ell;                          % µm（≈867）
    tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
    dhat = tip ./ vecnorm(tip);                                    % 3×6 measure 單位方向
    Q    = ell * dhat;                                             % 6 顆等效磁荷位置 [µm]

    draw_scene(ell, dhat, Q, figdir, false, 0,   'frames_lattice_3d.png');
    draw_scene(ell, dhat, Q, figdir, true,  150, 'frames_lattice_R150_3d.png');
end

function draw_scene(ell, dhat, Q, figdir, showR150, R150, outname)
    fig = figure('Color','w','Position',[80 80 940 880]); hold on;

    %% 半透明球殼（半徑 ℓ̂）
    [sx,sy,sz] = sphere(40);
    surf(ell*sx, ell*sy, ell*sz, 'FaceColor',[0.6 0.7 0.85], 'FaceAlpha',0.10, ...
         'EdgeColor',[0.7 0.75 0.82], 'EdgeAlpha',0.25, 'LineWidth',0.3);

    %% R=150µm 校正取樣球殼（綠；線框樣式同大球，才看得出是 3D 球）
    if showR150
        surf(R150*sx, R150*sy, R150*sz, 'FaceColor',[0.45 0.78 0.55], 'FaceAlpha',0.12, ...
             'EdgeColor',[0.35 0.68 0.45], 'EdgeAlpha',0.35, 'LineWidth',0.3);
    end

    %% WP → charge 虛線（6 條，墊底）+ ℓ̂ 標記
    for k = 1:6
        plot3([0 Q(1,k)],[0 Q(2,k)],[0 Q(3,k)], '--', 'Color',[0.45 0.45 0.45], 'LineWidth',1.2);
    end
    mid = 0.50*Q(:,2);                                           % ℓ̂ 標在 P2 線（無 actuator 軸、右上開闊、遠離 WP）
    text(mid(1)+70, mid(2), mid(3), '$\hat{\ell}$', 'Interpreter','latex', ...
         'FontSize',22,'FontWeight','bold','Color',[0.15 0.15 0.15]);

    %% 兩層「同平面」三角面：下極 P1/P3/P6（z=−ℓ̂·c）、上極 P2/P4/P5（z=+ℓ̂·c）
    lo = [1 3 6];  up = [2 4 5];
    fill3(Q(1,lo),Q(2,lo),Q(3,lo), [0.95 0.55 0.45], 'FaceAlpha',0.20, 'EdgeColor',[0.80 0.35 0.25], 'LineWidth',1.3);
    fill3(Q(1,up),Q(2,up),Q(3,up), [0.45 0.65 0.95], 'FaceAlpha',0.20, 'EdgeColor',[0.25 0.45 0.85], 'LineWidth',1.3);

    %% measure 三軸（全域 x,y,z；灰）
    Lm = 0.50*ell;  cm = [0.25 0.25 0.25];  mlab = {'x_m','y_m','z_m'};
    Im = eye(3);
    for k = 1:3
        v = Lm*Im(:,k);
        quiver3(0,0,0, v(1),v(2),v(3), 0, 'Color',cm, 'LineWidth',2.5, 'MaxHeadSize',0.35);
        t = v*1.22;
        text(t(1),t(2),t(3), mlab{k}, 'FontSize',15,'FontWeight','bold','Color',cm, ...
             'HorizontalAlignment','center');
    end

    %% actuator 三軸（d̂_P1/d̂_P3/d̂_P5，指向 P1/P3/P5；藍）
    ca = [0.10 0.35 0.85];  alab = {'x_a','y_a','z_a'};  apole = [1 3 5];
    for k = 1:3
        v = 0.50*ell*dhat(:,apole(k));                            % 沿 d̂_Pk 但不到極點（避免 label 撞 Pk）
        quiver3(0,0,0, v(1),v(2),v(3), 0, 'Color',ca, 'LineWidth',2.5, 'MaxHeadSize',0.42);
        t = v*1.26;
        text(t(1),t(2),t(3), alab{k}, 'FontSize',15,'FontWeight','bold','Color',ca, ...
             'HorizontalAlignment','center');
    end

    %% 6 極點（等效磁荷 ℓ̂·d̂）+ P1..P6 標籤
    plot3(Q(1,:),Q(2,:),Q(3,:), 'o', 'MarkerSize',9, ...
          'MarkerFaceColor',[0.85 0.15 0.35], 'MarkerEdgeColor','k', 'LineWidth',1.0);
    for k = 1:6
        t = Q(:,k) + 155*dhat(:,k);
        text(t(1),t(2),t(3), sprintf('P%d',k), 'FontSize',14,'FontWeight','bold', ...
             'Color',[0.6 0 0.2], 'HorizontalAlignment','center');
    end

    %% WP 原點（只留 + 標記，不標 WP 字）
    plot3(0,0,0,'k+','MarkerSize',13,'LineWidth',2.2);

    %% 3D 框（box on + daspect：cube 幾何、view-robust）
    bh = 1000;
    grid off; box on; daspect([1 1 1]);
    xlim([-bh bh]); ylim([-bh bh]); zlim([-bh bh]);
    view(120, 25);
    set(gca,'FontSize',14,'FontWeight','bold','LineWidth',1.5);
    set(gca,'XTick',-1000:500:1000,'YTick',-1000:500:1000,'ZTick',-1000:500:1000);
    xlabel('x_m (\mum)','FontWeight','bold');
    ylabel('y_m (\mum)','FontWeight','bold');
    zlabel('z_m (\mum)','FontWeight','bold');
    ax = gca; ax.Toolbar.Visible = 'off';
    hold off;

    out = fullfile(figdir, outname);
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('saved %s\n', out);
    close(fig);
end
