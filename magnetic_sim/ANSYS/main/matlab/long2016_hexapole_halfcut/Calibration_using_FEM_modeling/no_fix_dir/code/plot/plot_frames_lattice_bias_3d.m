function plot_frames_lattice_bias_3d()
% PLOT_FRAMES_LATTICE_BIAS_3D  measure/actuator 座標系 + 固定 ℓ̂ 球殼 + 6 極 bias 磁荷（no_fix）。
%   風格同 fix_dir 的 frames_lattice_R150_3d，唯一差別：6 顆等效磁荷改用 18-param bias 模型
%   實際位置 pc = R_act'·(ℓ̂·(Pc_base+E(ê)))，**會離開固定的 ℓ̂ 球殼**（fix 版 6 點都在殼上）。
%   球殼固定：只畫一顆半徑 = ℓ̂（no_fix 自己的 ell）的參考殼，不隨點縮放；點相對它移動、離殼。
%   出兩張：frames_lattice_bias_3d.png（只 ℓ̂ 殼）、frames_lattice_bias_R150_3d.png（多綠 R=150µm 球）。
%   單位 µm；框 = box on + daspect（cube、view-robust）。

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % mt_constants
    here   = fileparts(mfilename('fullpath'));
    nofix  = fileparts(fileparts(here));                            % .../no_fix_dir
    addpath(fullfile(nofix,'code','function'));                     % make_Pc
    figdir = fullfile(nofix,'figures');   if ~exist(figdir,'dir'); mkdir(figdir); end

    %% ---- 載 bias fit（ell µm、ê 17×1）+ 重建 6 荷位置（不需 FEM）----
    S = load(fullfile(nofix,'data','fit_bias_R150um_gap200um_mueq.mat'), 'ell','e_hat');
    ell = S.ell;   e_hat = S.e_hat;                                 % ell µm（≈857）
    cnst = mt_constants();
    tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
    dhat = tip ./ vecnorm(tip);                                     % 3×6 measure 單位方向
    R_act   = [dhat(:,1), dhat(:,3), dhat(:,5)].';                  % actuator→measure（列=P1/P3/P5）
    Pc_base = [1 -1 0 0 0 0; 0 0 1 -1 0 0; 0 0 0 0 1 -1];
    Pc = make_Pc(e_hat, Pc_base);                                   % 3×6 actuator 正規化（含 bias）
    pc = R_act.' * (ell * Pc);                                      % measure 框 [µm]；E=0 時 = ell·dhat

    off_um = vecnorm(pc - ell*dhat, 2, 1);                          % 每極離殼(相對 fix 在軸位)偏移 [µm]
    rad_um = vecnorm(pc, 2, 1);                                     % 每極磁荷半徑 [µm]（≠ ell = 離殼）
    fprintf('ell(shell)=%.1f µm；6 荷半徑 |pc| = %s µm\n', ell, mat2str(round(rad_um),5));
    fprintf('每極偏移 Δ = %s µm（max %.1f）\n', mat2str(round(off_um),4), max(off_um));

    draw_scene(ell, dhat, pc, figdir, false, 0,   'frames_lattice_bias_3d.png');
    draw_scene(ell, dhat, pc, figdir, true,  150, 'frames_lattice_bias_R150_3d.png');
end

function draw_scene(ell, dhat, pc, figdir, showR150, R150, outname)
    fig = figure('Color','w','Position',[80 80 940 880]); hold on;

    %% 固定 ℓ̂ 半透明球殼（半徑 = ell，不隨點縮放）
    [sx,sy,sz] = sphere(40);
    surf(ell*sx, ell*sy, ell*sz, 'FaceColor',[0.6 0.7 0.85], 'FaceAlpha',0.10, ...
         'EdgeColor',[0.7 0.75 0.82], 'EdgeAlpha',0.25, 'LineWidth',0.3);

    %% R=150µm 校正取樣球殼（綠）
    if showR150
        surf(R150*sx, R150*sy, R150*sz, 'FaceColor',[0.45 0.78 0.55], 'FaceAlpha',0.12, ...
             'EdgeColor',[0.35 0.68 0.45], 'EdgeAlpha',0.35, 'LineWidth',0.3);
    end

    %% WP → bias 磁荷 虛線（6 條，墊底，指向離殼的 pc）
    for k = 1:6
        plot3([0 pc(1,k)],[0 pc(2,k)],[0 pc(3,k)], '--', 'Color',[0.45 0.45 0.45], 'LineWidth',1.2);
    end

    %% 兩層「同平面」三角面：下極 P1/P3/P6（紅）、上極 P2/P4/P5（藍）（用 bias 位置）
    lo = [1 3 6];  up = [2 4 5];
    fill3(pc(1,lo),pc(2,lo),pc(3,lo), [0.95 0.55 0.45], 'FaceAlpha',0.20, 'EdgeColor',[0.80 0.35 0.25], 'LineWidth',1.3);
    fill3(pc(1,up),pc(2,up),pc(3,up), [0.45 0.65 0.95], 'FaceAlpha',0.20, 'EdgeColor',[0.25 0.45 0.85], 'LineWidth',1.3);

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
        v = 0.50*ell*dhat(:,apole(k));
        quiver3(0,0,0, v(1),v(2),v(3), 0, 'Color',ca, 'LineWidth',2.5, 'MaxHeadSize',0.42);
        t = v*1.26;
        text(t(1),t(2),t(3), alab{k}, 'FontSize',15,'FontWeight','bold','Color',ca, ...
             'HorizontalAlignment','center');
    end

    %% 6 顆 bias 磁荷（離殼）+ P1..P6 標籤
    plot3(pc(1,:),pc(2,:),pc(3,:), 'o', 'MarkerSize',9, ...
          'MarkerFaceColor',[0.85 0.15 0.35], 'MarkerEdgeColor','k', 'LineWidth',1.0);
    for k = 1:6
        t = pc(:,k) + 155*dhat(:,k);
        text(t(1),t(2),t(3), sprintf('P%d',k), 'FontSize',14,'FontWeight','bold', ...
             'Color',[0.6 0 0.2], 'HorizontalAlignment','center');
    end

    %% WP 原點
    plot3(0,0,0,'k+','MarkerSize',13,'LineWidth',2.2);

    %% 3D 框（box on + daspect：cube、view-robust）
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
