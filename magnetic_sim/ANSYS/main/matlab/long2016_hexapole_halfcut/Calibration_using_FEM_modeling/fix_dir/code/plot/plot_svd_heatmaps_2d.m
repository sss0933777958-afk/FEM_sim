function plot_svd_heatmaps_2d()
% PLOT_SVD_HEATMAPS_2D  long2016 gain(C)/iso(κ) 的「連續極座標熱圖」（取代舊離散節點散點）。
%   C、κ 是位置相關的平滑擬合函數（charge model：C=∏σ_k、κ=σ₃/σ₁；非 raw FEM 場）→
%   在「三個 actuator 座標參考切面」上連續評估（極座標網格 r×θ、非 FEM 節點），
%   切面沿 actuator 軸 ea_x=P1、ea_y=P3、ea_z=P5 尖端方向（magic-angle 兩兩正交，與 plot_ref_planes_3d 一致）：
%     ref_{x_ay_a}：p=r(ea_x·cosθ + ea_y·sinθ)
%     ref_{y_az_a}：p=r(ea_y·cosθ + ea_z·sinθ)
%     ref_{x_az_a}：p=r(ea_x·cosθ + ea_z·sinθ)
%   sv=svd((p/ℓ̂−d̂)./|p/ℓ̂−d̂|³·Ĥ_I)（C、κ 為 rotation-invariant 純量場，Pc 維持 measure 框不旋轉）；
%   畫成極座標 pcolor（θ 0–360°、r 0–500µm）。
%   輸出 6 檔：fix_dir/figures/{gain,iso}_polar_{xaya,yaza,xaza}.png。

    here    = fileparts(mfilename('fullpath'));
    fixdir  = fileparts(fileparts(here));
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');   % mt_constants
    figdir  = fullfile(fixdir,'figures');  if ~exist(figdir,'dir'); mkdir(figdir); end

    cnst = mt_constants();                                 % long
    tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];  dhat = tip./vecnorm(tip);  % 3×6 measure 框
    ea_x = dhat(:,1);  ea_y = dhat(:,3);  ea_z = dhat(:,5);  % actuator 三軸（P1/P3/P5）
    S = load(fullfile(fixdir,'data','fit_fixl_R150um_gap_200um.mat'), 'ell','gB','Khat');
    ell_m = S.ell*1e-6;  Hhat = S.gB*S.Khat;

    make_polar(ell_m, Hhat, dhat, ea_x, ea_y, ea_z, figdir);
end

function make_polar(ell_m, Hhat, Pc, ea_x, ea_y, ea_z, figdir)
    R = 500;                                              % µm
    r_um = 0:3:R;   th = (0:2:360)*pi/180;
    [RR, TH] = meshgrid(r_um, th);                        % [nθ × nr]
    dhat = Pc;                                            % 3×6 極尖端單位向量（標記用）
    [Cxy,Kxy] = ck_grid(RR, TH, ea_x, ea_y, ell_m, Hhat, Pc);   % x_a-y_a 面
    [Cyz,Kyz] = ck_grid(RR, TH, ea_y, ea_z, ell_m, Hhat, Pc);   % y_a-z_a 面
    [Cxz,Kxz] = ck_grid(RR, TH, ea_x, ea_z, ell_m, Hhat, Pc);   % x_a-z_a 面

    Clab = '$\mathcal{C}\;[(\mathrm{mT/A})^{3}]$';   Klab = '$\kappa$';
    % 各面：(u,v)=面內基底、w=面法向（判 in-plane 極）。標記在 rim r=R 標對應磁極。
    render_polar(RR, TH, Cxy, Clab, R, false, dhat, ea_x, ea_y, ea_z, fullfile(figdir,'gain_polar_xaya.png'));
    render_polar(RR, TH, Kxy, Klab, R, false,  dhat, ea_x, ea_y, ea_z, fullfile(figdir,'iso_polar_xaya.png'));
    render_polar(RR, TH, Cyz, Clab, R, false, dhat, ea_y, ea_z, ea_x, fullfile(figdir,'gain_polar_yaza.png'));
    render_polar(RR, TH, Kyz, Klab, R, false,  dhat, ea_y, ea_z, ea_x, fullfile(figdir,'iso_polar_yaza.png'));
    render_polar(RR, TH, Cxz, Clab, R, false, dhat, ea_x, ea_z, ea_y, fullfile(figdir,'gain_polar_xaza.png'));
    render_polar(RR, TH, Kxz, Klab, R, false,  dhat, ea_x, ea_z, ea_y, fullfile(figdir,'iso_polar_xaza.png'));
end

function [C,K] = ck_grid(RR, TH, u, v, ell_m, Hhat, Pc)
% 逐格解析評估 C=∏σ、κ=σ₃/σ₁，切面由基底對 (u,v) 定義：p=r(u·cosθ + v·sinθ)。RR µm、TH rad。
    C = zeros(size(RR));  K = zeros(size(RR));
    for a = 1:numel(RR)
        r = RR(a)*1e-6;  t = TH(a);
        p = r*(u*cos(t) + v*sin(t));                      % 3×1（actuator 面上一點，measure 座標）
        Dk = p/ell_m - Pc;                                % 3×6
        sv = svd((Dk ./ (vecnorm(Dk).^3)) * Hhat);
        C(a) = prod(sv);  K(a) = sv(3)/sv(1);
    end
end

function render_polar(RR, TH, val, clab, R, flipcmap, dhat, u, v, w, out)
    X = RR.*cos(TH);  Y = RR.*sin(TH);
    fig = figure('Color','w','Position',[80 80 980 860]);
    ax  = axes(fig);  hold(ax,'on');
    surf(ax, X, Y, zeros(size(X)), val, 'EdgeColor','none');   % 連續填色（z=0 flat, color=val）
    view(ax,2);  shading(ax,'interp');
    cmap = jet;  if flipcmap, cmap = flipud(cmap); end
    colormap(ax, cmap);  caxis(ax, [min(val(:)) max(val(:))]);
    axis(ax,'equal');  axis(ax,'off');

    % ---- 極座標網格 overlay（畫在高 z、view(2) 蓋在填色之上）----
    zt  = max(val(:)) + 1;
    thg = linspace(0, 2*pi, 240);
    rla = 100*pi/180;                                          % 徑向標籤角度（避開 spoke）
    for rr = (0.25:0.25:1)*R
        plot3(ax, rr*cos(thg), rr*sin(thg), zt*ones(size(thg)), '-', 'Color',[.35 .35 .35], 'LineWidth',0.9);
        text(ax, rr*cos(rla), rr*sin(rla), zt, sprintf('%d\\mum',round(rr)), ...     % 環上標半徑
             'FontSize',12,'FontWeight','bold','Color','k','HorizontalAlignment','center', ...
             'BackgroundColor','w','Margin',1);
    end
    for a = 0:30:330
        ar = a*pi/180;
        plot3(ax, [0 R*cos(ar)], [0 R*sin(ar)], [zt zt], '-', 'Color',[.35 .35 .35], 'LineWidth',0.8);
        text(ax, 1.09*R*cos(ar), 1.09*R*sin(ar), sprintf('%d\\circ',a), ...
             'HorizontalAlignment','center','FontSize',15,'FontWeight','bold');
    end

    % ---- 磁極投影標記：面內極（|d̂·w|≈0）在 rim r=R 標紅點 + P k（尖端剛好落 R_norm 面）----
    for k = 1:size(dhat,2)
        pu = dhat(:,k).'*u;  pv = dhat(:,k).'*v;  pw = dhat(:,k).'*w;
        if abs(pw) < 0.1                                     % 面內極（門檻鬆到 0.1：面內殘差~1e-4 vs 面外~1.0；hung magic_angle=54.74 四捨五入殘差 1.3e-4）
            tk = atan2(pv, pu);                              % 面內角
            plot3(ax, R*cos(tk), R*sin(tk), zt+1, 'o', 'MarkerSize',11, ...
                  'MarkerFaceColor',[0.90 0.90 0.90], 'MarkerEdgeColor',[0.5 0 0.12], 'LineWidth',1.6);
            text(ax, 0.86*R*cos(tk), 0.86*R*sin(tk), zt+1, sprintf('P%d',k), ...
                 'FontSize',15,'FontWeight','bold','Color',[0.5 0 0.12], ...
                 'HorizontalAlignment','center','BackgroundColor','w','Margin',1, ...
                 'EdgeColor',[0.5 0 0.12]);
        end
    end
    xlim(ax,[-1.20*R 1.20*R]);  ylim(ax,[-1.20*R 1.20*R]);

    cint = 'tex'; if startsWith(strtrim(clab),'$'), cint = 'latex'; end
    cb = colorbar(ax); cb.FontSize = 16; cb.FontWeight = 'bold';
    cb.Label.Interpreter = cint; cb.Label.String = clab; cb.Label.FontWeight='bold'; cb.Label.FontSize=16;
    if numel(cb.Ticks) >= 4, cb.Ticks = cb.Ticks(1:2:end); end
    ax.Toolbar.Visible = 'off';
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('saved %s\n', out);
    close(fig);
end
