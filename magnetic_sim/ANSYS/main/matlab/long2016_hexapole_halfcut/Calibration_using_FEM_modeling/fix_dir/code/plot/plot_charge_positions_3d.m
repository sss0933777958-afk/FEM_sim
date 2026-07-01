function plot_charge_positions_3d()
% PLOT_CHARGE_POSITIONS_3D  fix-l 點電荷模型：P1 & P2 磁荷位置示意 3D 圖。
%   畫 P1(下極,halfcut 半錐) 與 P2(上極,全錐) 的極輪廓 + 粉色等效磁荷點
%   (q = ell_hat * dhat, 落在極軸上, 距 WP = ell_hat) + 磁荷->WP 虛線(標 ell_hat)
%   + WP marker。無磁路箭頭、無場資料(純幾何)。
%   view(-60,-30)、①粗體框圖風格。ell_hat 取 gap200um_mueq(canonical) 擬合值。
%   輸出 fix_dir/figures/charge_positions_P1P2_3d.png。

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    cnst = mt_constants();

    %% 路徑(相對自身)：.../fix_dir/code/plot -> fix_dir
    here   = fileparts(mfilename('fullpath'));
    fixdir = fileparts(fileparts(here));
    matf   = fullfile(fixdir, 'data', 'fit_fixl_R150um_gap200um_mueq.mat');
    figdir = fullfile(fixdir, 'figures');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    S = load(matf, 'ell');  ell = S.ell;
    % ell 單位在單位改版後 = µm(~867)；舊檔可能是 m(~9e-4) 或 mm(~0.87) → 正規化到 mm
    if     ell > 100,  ellmm = ell*1e-3;   % µm
    elseif ell > 0.1,  ellmm = ell;        % mm
    else               ellmm = ell*1e3;    % m
    end
    fprintf('ell_hat = %.4f mm (%.1f um), variant=gap200um_mueq\n', ellmm, ellmm*1e3);

    tip = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];   % 3x6 [m], WP frame
    L   = cnst.POLE_CONE_LEN*1e3;  R = cnst.POLE_R*1e3;             % mm

    fig = figure('Color','w','Position',[80 80 900 820]); hold on;

    for p = [1 2]
        tp   = tip(:,p);  dhat = tp/norm(tp);
        ax   = cnst.pole_axis(:,p);  ax = ax/norm(ax);
        apex = tp*1e3;                 % 極尖端 (mm)
        chg  = ellmm*dhat;             % 等效磁荷點 (mm) = ell_hat * dhat

        %% 極輪廓 cone (下極 halfcut = 半錐 + 平切面；上極 = 全錐) —— 抄單極範本
        isLower = logical(cnst.pole_is_lower(p));
        tmax = min(1, 2.6/L);  tt = linspace(0,tmax,12);
        if isLower
            v1 = cross(ax,[0;0;1]); v1 = v1/norm(v1);  v2 = [0;0;1];
            ph = linspace(pi, 2*pi, 22);
        else
            vN = null(ax.'); v1 = vN(:,1); v2 = vN(:,2);  ph = linspace(0,2*pi,40);
        end
        [TT,PP] = meshgrid(tt,ph);
        CX = apex(1)+L*TT*ax(1)+(R*TT).*(cos(PP)*v1(1)+sin(PP)*v2(1));
        CY = apex(2)+L*TT*ax(2)+(R*TT).*(cos(PP)*v1(2)+sin(PP)*v2(2));
        CZ = apex(3)+L*TT*ax(3)+(R*TT).*(cos(PP)*v1(3)+sin(PP)*v2(3));
        surf(CX,CY,CZ,'FaceColor',[0.55 0.55 0.55],'FaceAlpha',0.30,'EdgeColor','none');
        if isLower
            bc = apex + L*tmax*ax;  e1 = bc + R*tmax*v1;  e2 = bc - R*tmax*v1;
            fill3([apex(1) e1(1) e2(1)],[apex(2) e1(2) e2(2)],[apex(3) e1(3) e2(3)], ...
                  [0.55 0.55 0.55],'FaceAlpha',0.30,'EdgeColor','none');
        end

        %% 磁荷 -> WP 虛線 + ell_hat 標註（P1/P2 都在 y=0 面，ℓ̂ 往 +y 抬離避免與 q/WP 重疊）
        plot3([chg(1) 0],[chg(2) 0],[chg(3) 0],'--','Color',[0.35 0.35 0.35],'LineWidth',1.6);
        mid = chg/2;
        text(mid(1), mid(2)+0.40, mid(3), '$\hat{\ell}$', 'Interpreter','latex', ...
             'FontSize',20,'FontWeight','bold','Color',[0.15 0.15 0.15], ...
             'HorizontalAlignment','center');

        %% 粉色磁荷點 + label（q 標籤沿磁極軸往外移，遠離中心 WP）
        plot3(chg(1),chg(2),chg(3),'o','MarkerSize',8, ...
              'MarkerFaceColor',[1 0.41 0.71],'MarkerEdgeColor','k','LineWidth',1.2);
        qpos = chg + 0.55*dhat;
        text(qpos(1), qpos(2), qpos(3), sprintf('q_{P%d}',p), ...
             'FontSize',16,'FontWeight','bold','Color',[0.75 0 0.45], ...
             'HorizontalAlignment','center');
    end

    %% WP marker (不畫磁路箭頭)
    plot3(0,0,0,'k+','MarkerSize',14,'LineWidth',2.2);
    text(0,0,0.18,'WP','FontSize',15,'FontWeight','bold');

    %% 3D 框體 — 手動畫框邊黑色，省略「最遠角」的 3 條（會橫穿內部的雜線），其餘 9 邊全黑等粗
    bh = 2.6;
    grid off; box off; daspect([1 1 1]);                     % daspect 等比且不動 explicit limits
    xlim([-bh bh]); ylim([-bh bh]); zlim([-bh bh]);
    view(-30,-20);   % 方位角 −30、仰角 −20（先設好 view，campos 才正確）
    set(gca,'FontSize',13,'FontWeight','bold','LineWidth',1.5);
    set(gca,'XTick',-2:1:2,'YTick',-2:1:2,'ZTick',-2:1:2);  % 三軸刻度一致（z 跟 x/y 同尺度）
    draw_box_edges(bh, 1.5);                                 % 省略最遠角 3 邊（頂部橫穿內部雜線），其餘 9 邊黑
    hold off;
    xlabel('x (mm)','FontWeight','bold');
    ylabel('y (mm)','FontWeight','bold');
    zlabel('z (mm)','FontWeight','bold');

    out = fullfile(figdir,'charge_positions_P1P2_3d.png');
    exportgraphics(fig, out, 'Resolution',600);
    fprintf('saved %s\n', out);
end

function draw_box_edges(bh, lw)
% 畫立方體框邊(±bh) 黑色 LineWidth lw，但**省略離相機最遠那個角的 3 條邊**
% （那 3 條會從頂部橫穿內部變雜線；其餘 9 邊 = 外框輪廓 + 近角邊，全黑等粗）
    s = [-bh bh];  [A,B,C] = ndgrid(s,s,s);  corners = [A(:) B(:) C(:)];   % 8 角
    cp = campos;  [~,ifar] = max(vecnorm(corners - cp, 2, 2));  far = corners(ifar,:);  % 最遠角
    for i = 1:8
        for j = i+1:8
            if sum(abs(corners(i,:)-corners(j,:)) > 1e-9) == 1              % 差一個座標 = 一條邊
                if all(abs(corners(i,:)-far)<1e-9) || all(abs(corners(j,:)-far)<1e-9)
                    continue;                                              % 跳過與最遠角相連的 3 邊
                end
                plot3([corners(i,1) corners(j,1)], [corners(i,2) corners(j,2)], ...
                      [corners(i,3) corners(j,3)], 'k-', 'LineWidth', lw);
            end
        end
    end
end
