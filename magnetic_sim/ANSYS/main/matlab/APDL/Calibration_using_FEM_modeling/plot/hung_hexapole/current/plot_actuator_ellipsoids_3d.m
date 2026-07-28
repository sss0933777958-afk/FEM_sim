function plot_actuator_ellipsoids_3d()
% PLOT_ACTUATOR_ELLIPSOIDS_3D  actuator 座標系（真實傾斜 d̂ 方向）+ WP 中心 + 近 P1 tip 的致動橢球。
%   在 measure frame 繪圖（同 frames_lattice / plot_ellipsoids_yaxis_3d）：actuator 三軸畫成真實
%   magic-angle 傾斜方向 x_a=d̂_P1 / y_a=d̂_P3 / z_a=d̂_P5（指向 P1/P3/P5），非正交框對齊。
%   致動橢球 = 單位電流球經 T(p)=(Dk/|Dk|³)·Ĥ_I 的像：半軸 σ_k（SVD）沿 U(:,k)（measure frame）。
%   兩點：WP 中心 p=0（近球形/等向）與 近 P1 tip p=350µm·d̂_P1（沿 x_a 傾斜拉長）。共用尺度。
%   只呈現 橢球 + actuator 三軸 + P1 tip（半透明+淡網格+光照＝立體實心感）。
%   ★ model-derived（hung）。box on + daspect（不變形）、view(120,25)。輸出 current_base/figures/actuator_frame_ellipsoids_3d.png。

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');   % hung mt_constants
    cnst = mt_constants();
    here   = fileparts(mfilename('fullpath'));
    fixdir = fileparts(fileparts(here));
    figdir = fullfile(fixdir,'figures','single_param');   if ~exist(figdir,'dir'); mkdir(figdir); end

    S = load(fullfile(fixdir,'data','fit_fixl_R150um_gap_200um.mat'),'ell','gB','Khat');
    ell_m = S.ell*1e-6;   Hhat = S.gB*S.Khat;                        % ell µm→m；Ĥ_I=gB·K̂（電流側、frame 無關）
    tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
    dhat = tip ./ vecnorm(tip);                                      % 3×6 measure 單位方向（電荷/軸/tip 都用它）
    Rn   = 500;                                                      % R_norm [µm]（極尖端半徑）

    %% ---- 兩取樣點（measure）：WP 中心、近 P1 tip（沿 d̂_P1 350µm）----
    C_um = [ [0;0;0], 350*dhat(:,1) ];                              % 3×2 橢球中心 [µm]
    npt  = size(C_um,2);
    Uc = cell(1,npt);  SVc = cell(1,npt);
    for i = 1:npt
        p  = C_um(:,i)*1e-6;
        Dk = p/ell_m - dhat;                                        % 3×6（measure）
        [U,Sig,~] = svd((Dk ./ vecnorm(Dk).^3) * Hhat);
        Uc{i} = U;  SVc{i} = diag(Sig);                            % U 在 measure frame
    end
    s1max = max(cellfun(@(s) s(1), SVc));
    sc    = 150 / s1max;                                           % 共用尺度：最大半軸 ≈ 150µm

    %% ---- 畫圖（精簡：只 橢球 + actuator 三軸 + P1 tip + R_norm 球殼）----
    [ex,ey,ez] = sphere(26);  sz0 = size(ex);
    fig = figure('Color','w','Position',[70 70 940 900]); hold on;

    % R_norm=500µm 工作球殼（淡線框、墊底；P1 tip 落在球面；FaceLighting none 不被 camlight 打亮）
    [gx,gy,gz] = sphere(40);
    surf(Rn*gx, Rn*gy, Rn*gz, 'FaceColor',[0.60 0.70 0.85], 'FaceAlpha',0.08, ...
         'EdgeColor',[0.70 0.75 0.82], 'EdgeAlpha',0.20, 'LineWidth',0.3, 'FaceLighting','none');

    % 淡虛線 WP→P1 tip（沿 d̂_P1）
    Ptip = Rn*dhat(:,1);
    plot3([0 Ptip(1)],[0 Ptip(2)],[0 Ptip(3)], '--', 'Color',[0.62 0.62 0.62], 'LineWidth',1.2);

    % actuator 三軸（藍；由原點沿真實傾斜方向 d̂_P1/d̂_P3/d̂_P5）— 先畫，讓半透明橢球疊其上
    ca = [0.10 0.35 0.85];  alab = {'x_a','y_a','z_a'};  apole = [1 3 5];  La = 400;
    aoff = [ 0 130 0 ; 0 60 0 ; 0 0 55 ].';                        % 標籤偏移（x_a 往 +y 離開 y=0 擁擠面）
    for k = 1:3
        v = La*dhat(:,apole(k));
        quiver3(0,0,0, v(1),v(2),v(3), 0, 'Color',ca, 'LineWidth',2.6, 'MaxHeadSize',0.4);
        t = v*1.12 + aoff(:,k);
        text(t(1),t(2),t(3), alab{k}, 'FontSize',17,'FontWeight','bold','Color',ca, ...
             'HorizontalAlignment','center');
    end

    % 兩顆半透明立體橢球（lighting 明暗 + 淡網格 → 立體、不擋軸；無內部主軸箭頭）
    for i = 1:npt
        c = C_um(:,i);  U = Uc{i};  sv = SVc{i};
        E = U * diag(sv*sc) * [ex(:) ey(:) ez(:)].';               % 3×N（半軸 σ_k·sc 沿 U，measure）
        surf(reshape(E(1,:),sz0)+c(1), reshape(E(2,:),sz0)+c(2), reshape(E(3,:),sz0)+c(3), ...
             'FaceColor',[0.40 0.58 0.86], 'FaceAlpha',0.55, ...
             'EdgeColor',[0.22 0.38 0.64], 'EdgeAlpha',0.20, 'LineWidth',0.3, ...
             'FaceLighting','gouraud', 'AmbientStrength',0.5, 'DiffuseStrength',0.8, ...
             'SpecularStrength',0.3);
    end

    % P1 tip 紅點 + 標籤
    plot3(Ptip(1),Ptip(2),Ptip(3), 'o', 'MarkerSize',12, 'MarkerFaceColor',[0.85 0.15 0.35], ...
          'MarkerEdgeColor','k', 'LineWidth',1.0);
    text(Ptip(1)+35, Ptip(2), Ptip(3)-40, 'P1 tip', 'FontSize',15,'FontWeight','bold', ...
         'Color',[0.6 0 0.2], 'HorizontalAlignment','center');

    % WP 原點 `+`（中心）
    plot3(0,0,0,'k+','MarkerSize',13,'LineWidth',2.2);

    % [ADDED] 兩顆橢球中心黑點 + 中心徑向距離標註（中心都在座標軸上：WP=原點、P1-tip=x_a 軸 350µm）
    %   標籤偏移：三 actuator 箭頭投影＝x_a→左下 / y_a→右下 / z_a→右上；左上為開放區。
    %   WP 標籤沿 −d̂_P3（投影到左上、避開 z_a 與橢球）；P1-tip 沿 +x 一點（左下開放）。
    coff = { 270*(-dhat(:,3)) + [0;0;20],  [150;-20;-25] };
    for i = 1:npt
        c  = C_um(:,i);
        rr = norm(c);                                              % 中心離 WP 原點徑向距離 [µm]（0 / 350）
        o  = coff{i};
        plot3(c(1),c(2),c(3), 'o', 'MarkerSize',9, 'MarkerFaceColor','k', ...
              'MarkerEdgeColor','k', 'LineWidth',0.8);
        text(c(1)+o(1), c(2)+o(2), c(3)+o(3), sprintf('r = %.0f \\mum', rr), ...
             'FontSize',13,'FontWeight','bold','Color','k','HorizontalAlignment','center');
    end

    %% ---- 框 / 風格（measure 座標 box；view 同 frames_lattice）----
    grid off; box on; daspect([1 1 1]);
    xlim([-560 560]); ylim([-560 560]); zlim([-560 560]);          % 含 R_norm=500µm 球殼 + 留邊
    view(120, 25);
    camlight('headlight');  camlight('left');  material dull;      % 光照讓橢球立體
    set(gca,'FontSize',15,'FontWeight','bold','LineWidth',2,'TickLength',[.012 .012]);
    set(gca,'XTick',-500:250:500,'YTick',-500:250:500,'ZTick',-500:250:500);   % 三軸標到 500µm
    xlabel('x_m (\mum)','FontWeight','bold');
    ylabel('y_m (\mum)','FontWeight','bold');
    zlabel('z_m (\mum)','FontWeight','bold');
    ax = gca; ax.Toolbar.Visible = 'off';
    hold off;

    out = fullfile(figdir,'actuator_frame_ellipsoids_3d.png');
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('sc=%.3f µm/(mT/A)；WP σ=[%.2f %.2f %.2f]、near-tip σ=[%.2f %.2f %.2f]；saved %s\n', ...
            sc, SVc{1}, SVc{2}, out);
    close(fig);
end
