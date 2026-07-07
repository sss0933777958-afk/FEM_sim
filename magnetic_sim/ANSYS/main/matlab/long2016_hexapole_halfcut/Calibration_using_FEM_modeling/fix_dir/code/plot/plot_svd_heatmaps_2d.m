function plot_svd_heatmaps_2d()
% PLOT_SVD_HEATMAPS_2D  long2016 gain(C)/iso(κ) 的「連續極座標熱圖」（取代舊離散節點散點）。
%   C、κ 是位置相關的平滑擬合函數（charge model：C=∏σ_k、κ=σ₃/σ₁；非 raw FEM 場）→
%   在兩個正交參考切面上「連續評估」（極座標網格 r×θ、非 FEM 節點）：
%     ref_xy = z_m=0 平面：p=r[cosθ,sinθ,0]
%     ref_xz = y_m=0 平面：p=r[cosθ,0,sinθ]
%   sv=svd((p/ℓ̂−d̂)./|p/ℓ̂−d̂|³·Ĥ_I)（measure 框、框無關）；畫成極座標 pcolor（θ 0–360°、r 0–500µm）。
%   輸出 4 檔：fix_dir/figures/{gain,iso}_polar_{xy,xz}.png。

    here    = fileparts(mfilename('fullpath'));
    fixdir  = fileparts(fileparts(here));
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');   % mt_constants
    figdir  = fullfile(fixdir,'figures');  if ~exist(figdir,'dir'); mkdir(figdir); end

    cnst = mt_constants();                                 % long
    tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];  dhat = tip./vecnorm(tip);  % 3×6 measure 框
    S = load(fullfile(fixdir,'data','fit_fixl_R150um_gap200um_mueq.mat'), 'ell','gB','Khat');
    ell_m = S.ell*1e-6;  Hhat = S.gB*S.Khat;

    make_polar(ell_m, Hhat, dhat, figdir);
end

function make_polar(ell_m, Hhat, Pc, figdir)
    R = 500;                                              % µm
    r_um = 0:3:R;   th = (0:2:360)*pi/180;
    [RR, TH] = meshgrid(r_um, th);                        % [nθ × nr]
    [Cxy,Kxy] = ck_grid(RR, TH, 'xy', ell_m, Hhat, Pc);
    [Cxz,Kxz] = ck_grid(RR, TH, 'xz', ell_m, Hhat, Pc);

    Clab = '$\mathcal{C}\;[(\mathrm{mT/A})^{3}]$';   Klab = '$\kappa$';
    render_polar(RR, TH, Cxy, Clab, R, false, fullfile(figdir,'gain_polar_xy.png'));
    render_polar(RR, TH, Kxy, Klab, R, false,  fullfile(figdir,'iso_polar_xy.png'));
    render_polar(RR, TH, Cxz, Clab, R, false, fullfile(figdir,'gain_polar_xz.png'));
    render_polar(RR, TH, Kxz, Klab, R, false,  fullfile(figdir,'iso_polar_xz.png'));
end

function [C,K] = ck_grid(RR, TH, plane, ell_m, Hhat, Pc)
% 逐格解析評估 C=∏σ、κ=σ₃/σ₁（plane='xy'→z=0；'xz'→y=0）。RR µm、TH rad。
    C = zeros(size(RR));  K = zeros(size(RR));
    for a = 1:numel(RR)
        r = RR(a)*1e-6;  t = TH(a);
        if strcmp(plane,'xy'), p = r*[cos(t);sin(t);0];
        else,                  p = r*[cos(t);0;sin(t)]; end
        Dk = p/ell_m - Pc;                                % 3×6
        sv = svd((Dk ./ (vecnorm(Dk).^3)) * Hhat);
        C(a) = prod(sv);  K(a) = sv(3)/sv(1);
    end
end

function render_polar(RR, TH, val, clab, R, flipcmap, out)
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
