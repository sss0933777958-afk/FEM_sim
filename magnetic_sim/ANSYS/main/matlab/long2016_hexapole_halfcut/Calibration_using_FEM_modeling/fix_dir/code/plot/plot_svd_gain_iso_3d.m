function plot_svd_gain_iso_3d()
% PLOT_SVD_GAIN_ISO_3D  逐點致動 SVD 的「單一山丘」surf（fix-ℓ, z=0 平面, l_bar≤150µm）。
%   每點 p 解析算 T_i(p)=S_i(p)·Ĥ_I（Ĥ_I=gB·Khat，已校正、框無關），SVD 取:
%     gain = ‖T‖_F = √(σ1²+σ2²+σ3²)（電流→場整體增益，mT/A）
%     iso  = σ3/σ1（均勻度，無因次；→1 越接近正球=各向一致）
%   座標（使用者定案）：z=0 水平面（無 φ）；2D 底面 = 極座標 (l_bar 半徑, θ 方位)，l_bar 0→150µm；
%     高度 z = gain / iso（山丘）；顏色 = 高度值本身（gain/iso 大小），colorbar 標 gain/iso。
%   渲染 Cartesian 圓盤：X=r cosθ、Y=r sinθ、Z=純量 → 單一 bumpy 山丘。physical p=r·[cosθ;sinθ;0]。
%   ★ model-derived（非 raw FEM）：T_i 解析式。gain/iso 是奇異值純量、框無關（用 measure 取樣）。
%   風格：surf 平滑 + jet；框體 = **box on + pbaspect([1 1 1])**（同 plot_upperP2P5_circuit_3d 的框法；
%         軸異質 µm vs gain/iso → 用 pbaspect 不用 daspect）、粗體黑框 LineWidth 2；無 legend。
%   輸出：fix_dir/figures/svd_gain_3d.png、svd_iso_3d.png。

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % mt_constants
    cnst = mt_constants();

    here   = fileparts(mfilename('fullpath'));
    fixdir = fileparts(fileparts(here));
    matf   = fullfile(fixdir, 'data', 'fit_fixl_R150um_gap200um_mueq.mat');
    figdir = fullfile(fixdir, 'figures');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    %% 校正參數 → Ĥ_I（框無關），ℓ̂
    S = load(matf, 'ell','gB','Khat');
    ell_m  = S.ell * 1e-6;
    Hhat_I = S.gB * S.Khat;                                          % 6×6，mT/A
    tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
    dhat = tip ./ vecnorm(tip);                                     % 3×6 measure 單位方向

    %% 極座標網格（z=0 平面）：半徑 l_bar、方位 θ
    rr_um = 0:3:150;   thd = 0:5:360;                               % l_bar[µm]、θ[deg]
    [RR, TH] = meshgrid(rr_um, thd*pi/180);                        % size = [nth, nr]
    X = RR .* cos(TH);   Y = RR .* sin(TH);                        % µm（Cartesian 渲染）
    gain = zeros(size(RR));   iso = zeros(size(RR));
    for a = 1:numel(RR)
        r = RR(a)*1e-6;
        p  = r * [cos(TH(a)); sin(TH(a)); 0];                      % z=0 平面物理點 [m]
        Dk = p/ell_m - dhat;                                       % 3×6
        sv = svd((Dk ./ (vecnorm(Dk).^3)) * Hhat_I);              % 3 個奇異值 [mT/A]
        gain(a) = norm(sv);   iso(a) = sv(1)/sv(3);   % P=√Σσ²；iso=σ_max/σ_min=σ1/σ3
    end
    fprintf('中心 r=0：P=%.3f、iso=%.4f（對照 18.23 / 1.061）\n', gain(1,1), iso(1,1));

    %% ---- 兩張單一山丘 surf（顏色=高度值；gain 符號 P、iso=σ1/σ3）------------
    render_hill(X, Y, gain, 'Gain (mT/A)', fullfile(figdir,'svd_gain_3d.png'), false);  % z 標籤 Gain（資料 √Σσ²）
    render_hill(X, Y, iso,  'iso',       fullfile(figdir,'svd_iso_3d.png'),  true);     % iso 翻色階（低=紅、高=藍）
end

function render_hill(X, Y, Z, zlab, outpng, flipcmap)
    if nargin < 6, flipcmap = false; end
    zmn = min(Z(:));  zmx = max(Z(:));
    fig = figure('Color','w','Position',[80 80 1000 880]);
    surf(X, Y, Z, 'EdgeColor','none');                            % CData = Z（顏色=gain/iso 大小）
    cmap = jet;  if flipcmap, cmap = flipud(cmap); end            % iso 翻色階
    shading interp; colormap(cmap); caxis([zmn zmx]);

    grid off; box on;                                             % 框體 = box on（同 circuit_3d 框法）
    xlim([-150 150]); ylim([-150 150]); zlim([zmn zmx]);
    view(-40, 30);
    pbaspect([1 1 1]);                                            % 立方 plot box（軸異質 → pbaspect 不用 daspect）
    set(gca,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);   % 粗體黑框
    set(gca,'XTick',-150:75:150,'YTick',-150:75:150);
    zt = get(gca,'ZTick');                                       % MATLAB 自動 z 刻度
    tol = 0.055*(zmx - zmn);                                     % 丟掉太靠近端點的自動刻度，避免與 min/max 標籤重疊
    zt = unique([zmn, zt(zt>zmn+tol & zt<zmx-tol), zmx]);       % 保留內部自動刻度 + 多加 min/max
    set(gca,'ZTick', zt, 'ZTickLabel', arrayfun(@(v) sprintf('%.4g',v), zt, 'UniformOutput',false));
    xlabel('x (\mum)','FontWeight','bold');
    ylabel('y (\mum)','FontWeight','bold');
    zlabel(zlab,'FontWeight','bold');

    cb = colorbar; cb.FontSize = 16; cb.FontWeight = 'bold';
    cb.Ticks = cb.Ticks(1:2:end);
    cb.Label.String = zlab; cb.Label.FontWeight = 'bold'; cb.Label.FontSize = 16;

    ax = gca; ax.Toolbar.Visible = 'off';                         % 匯出不帶 axes 工具列
    exportgraphics(fig, outpng, 'Resolution', 150);
    fprintf('saved %s\n', outpng);
    close(fig);
end
