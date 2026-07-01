function plot_charge_positions_3d()
% PLOT_CHARGE_POSITIONS_3D  no_fix_l（18-param bias）點電荷模型：P1 & P2 磁荷位置 3D 圖。
%   跟 fix_dir 同風格（P1 下極半錐 / P2 上極全錐 + 粉色磁荷點 + charge→WP 虛線標 ℓ̂ + WP，
%   手動框邊省略最遠角、box off + daspect、view −30/−20）。差異：磁荷取「18-param bias」擬合位置
%   pc = ℓ̂·(Pc_base + E(ê))（actuator frame → 轉 WP frame），並**標出相對於 fix_dir 的總偏移 Δ**。
%   gap200um_mueq、R=150µm；輸出 no_fix_dir/figures/charge_positions_P1P2_3d.png。

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');   % mt_constants...
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');  % ansys_path
    here   = fileparts(mfilename('fullpath'));
    nofix  = fileparts(fileparts(here));                    % .../no_fix_dir
    calroot= fileparts(nofix);                              % .../Calibration_using_FEM_modeling
    addpath(fullfile(nofix,'code','function'));             % load_coils_actuator/select_ball/fit_bias/make_Pc
    figdir = fullfile(nofix,'figures');  if ~exist(figdir,'dir'); mkdir(figdir); end

    cnst = mt_constants();
    tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];   % 3x6 [m] WP frame
    dhat = tip ./ vecnorm(tip);
    L = cnst.POLE_CONE_LEN*1e3;  R = cnst.POLE_R*1e3;                % mm

    %% ---- fix_dir 的 ℓ̂（canonical gap200um_mueq）→ fix 磁荷（WP 框, mm）----
    Sf = load(fullfile(calroot,'fix_dir','data','fit_fixl_R150um_gap200um_mueq.mat'),'ell');
    ell_fix = Sf.ell;  if ell_fix>100, ell_fix_mm = ell_fix*1e-3; elseif ell_fix>0.1, ell_fix_mm = ell_fix; else, ell_fix_mm = ell_fix*1e3; end
    pc_fix_mm = ell_fix_mm * dhat;                                   % 3x6 fix 磁荷 (WP, mm)

    %% ---- 重跑 no_fix 18-param bias 擬合 → no_fix 磁荷（WP 框, mm）----
    model = 'long2016_hexapole_halfcut';  apdl = [1,3,6,5,2,4];
    D = load_coils_actuator(model, cnst, apdl, 'all', 'gap200um_mueq');
    [P, Bstack] = select_ball(D, 150e-6);
    [ell_nofix, e_hat] = fit_bias(P, Bstack, D.Pc_base, 0.5e-3);     % ell_nofix [m]
    Pc      = make_Pc(e_hat, D.Pc_base);                            % 3x6 (actuator, 正規化)
    pc_meas = D.R_act.' * (ell_nofix * Pc);                         % actuator→WP [m]
    pc_mm   = pc_meas * 1e3;                                        % 3x6 no_fix 磁荷 (WP, mm)

    off_mm  = pc_mm - pc_fix_mm;  off_um = vecnorm(off_mm,2,1)*1e3;  % 每極總偏移 [µm]
    fprintf('ell_fix=%.1f µm, ell_nofix=%.1f µm\n', ell_fix_mm*1e3, ell_nofix*1e6);
    for p=[1 2], fprintf('  P%d 總偏移 Δ = %.1f µm\n', p, off_um(p)); end

    %% ---- 畫圖 ----
    fig = figure('Color','w','Position',[80 80 900 820]); hold on;
    for p = [1 2]
        tp = tip(:,p);  ax = cnst.pole_axis(:,p); ax = ax/norm(ax);
        apex = tp*1e3;  chg = pc_mm(:,p).';                        % no_fix 磁荷 (mm, 1x3)

        %% 極輪廓 cone（P1 下極 halfcut 半錐 + 平切面；P2 上極全錐）
        isLower = logical(cnst.pole_is_lower(p));
        tmax = min(1, 2.6/L);  tt = linspace(0,tmax,12);
        if isLower
            v1 = cross(ax,[0;0;1]); v1 = v1/norm(v1);  v2 = [0;0;1];  ph = linspace(pi,2*pi,22);
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

        %% 磁荷 → WP 虛線 + ℓ̂
        plot3([chg(1) 0],[chg(2) 0],[chg(3) 0],'--','Color',[0.35 0.35 0.35],'LineWidth',1.6);
        mid = chg/2;
        text(mid(1), mid(2)+0.40, mid(3), '$\hat{\ell}$', 'Interpreter','latex', ...
             'FontSize',20,'FontWeight','bold','Color',[0.15 0.15 0.15],'HorizontalAlignment','center');

        %% 粉色磁荷點 + q 標籤 + 總偏移 Δ 標註
        plot3(chg(1),chg(2),chg(3),'o','MarkerSize',8, ...
              'MarkerFaceColor',[1 0.41 0.71],'MarkerEdgeColor','k','LineWidth',1.2);
        cdir = chg/norm(chg);  qpos = chg + 0.55*cdir;
        text(qpos(1), qpos(2), qpos(3), sprintf('q_{P%d}',p), ...
             'FontSize',16,'FontWeight','bold','Color',[0.75 0 0.45],'HorizontalAlignment','center');
        dpos = qpos + [0, 0.6, -0.25];   % Δ 標註抬離 y=0 磁極軸平面，避免壓到磁荷點
        text(dpos(1), dpos(2), dpos(3), sprintf('\\Delta = %.0f \\mum', off_um(p)), ...
             'FontSize',13,'FontWeight','bold','Color',[0.1 0.1 0.1],'HorizontalAlignment','center');
    end

    %% WP marker
    plot3(0,0,0,'k+','MarkerSize',14,'LineWidth',2.2);
    text(0,0,0.18,'WP','FontSize',15,'FontWeight','bold');

    %% 3D 框體（同 fix：手動框邊、省略最遠角 3 邊、box off + daspect、三軸同刻度）
    bh = 2.6;
    grid off; box off; daspect([1 1 1]);
    xlim([-bh bh]); ylim([-bh bh]); zlim([-bh bh]);
    view(-30,-20);
    set(gca,'FontSize',13,'FontWeight','bold','LineWidth',1.5);
    set(gca,'XTick',-2:1:2,'YTick',-2:1:2,'ZTick',-2:1:2);
    draw_box_edges(bh, 1.5);
    hold off;
    xlabel('x (mm)','FontWeight','bold');
    ylabel('y (mm)','FontWeight','bold');
    zlabel('z (mm)','FontWeight','bold');

    out = fullfile(figdir,'charge_positions_P1P2_3d.png');
    exportgraphics(fig, out, 'Resolution',600);
    fprintf('saved %s\n', out);
end

function draw_box_edges(bh, lw)
% 畫立方體框邊(±bh) 黑色 LineWidth lw，省略離相機最遠那個角的 3 條邊（頂部橫穿內部雜線）
    s = [-bh bh];  [A,B,C] = ndgrid(s,s,s);  corners = [A(:) B(:) C(:)];
    cp = campos;  [~,ifar] = max(vecnorm(corners - cp, 2, 2));  far = corners(ifar,:);
    for i = 1:8
        for j = i+1:8
            if sum(abs(corners(i,:)-corners(j,:)) > 1e-9) == 1
                if all(abs(corners(i,:)-far)<1e-9) || all(abs(corners(j,:)-far)<1e-9), continue; end
                plot3([corners(i,1) corners(j,1)], [corners(i,2) corners(j,2)], ...
                      [corners(i,3) corners(j,3)], 'k-', 'LineWidth', lw);
            end
        end
    end
end
