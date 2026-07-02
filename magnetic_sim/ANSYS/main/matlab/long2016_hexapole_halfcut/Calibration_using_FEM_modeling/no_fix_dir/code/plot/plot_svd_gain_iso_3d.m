function plot_svd_gain_iso_3d()
% PLOT_SVD_GAIN_ISO_3D  逐點致動 SVD 的「單一山丘」surf（no_fix 18-param bias, z=0 平面, l_bar≤150µm）。
%   每點 p 解析算 T_i(p)=S_i(p)·Ĥ_I（Ĥ_I=gB·Khat，已校正、框無關），SVD 取:
%     gain = ‖T‖_F = √(σ1²+σ2²+σ3²)（電流→場整體增益，mT/A）
%     iso  = σ_max/σ_min = σ1/σ3（均勻度，無因次；→1 越接近正球=各向一致）
%   與 fix_dir 版唯一差別：電荷用 **bias 位置** dhat_bias = R_act'·(Pc_base+E(ê))（把 bias 電荷轉 measure 框），
%   其餘（z=0 極座標底面、jet、iso 翻色階、gain z 標籤、z 軸自動刻度+min/max）全同 fix_dir。
%   ★ model-derived。gain/iso 是奇異值純量、框無關。輸出 no_fix_dir/figures/svd_gain_3d.png、svd_iso_3d.png。

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % mt_constants
    cnst = mt_constants();

    here   = fileparts(mfilename('fullpath'));
    nofix  = fileparts(fileparts(here));
    addpath(fullfile(nofix,'code','function'));                     % make_Pc
    matf   = fullfile(nofix, 'data', 'fit_bias_R150um_gap200um_mueq.mat');
    figdir = fullfile(nofix, 'figures');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    %% 校正參數 → Ĥ_I（框無關）、ℓ̂、bias 電荷（轉 measure 框）
    S = load(matf, 'ell','gB','Khat','e_hat');
    ell_m  = S.ell * 1e-6;
    Hhat_I = S.gB * S.Khat;                                          % 6×6，mT/A
    tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
    dhat = tip ./ vecnorm(tip);                                     % 3×6 measure 單位方向
    R_act   = [dhat(:,1), dhat(:,3), dhat(:,5)].';
    Pc_base = [1 -1 0 0 0 0; 0 0 1 -1 0 0; 0 0 0 0 1 -1];
    Pc = make_Pc(S.e_hat, Pc_base);                                 % 3×6 actuator 正規化（含 bias）
    dhat_bias = R_act.' * Pc;                                       % 3×6 measure 框 bias 電荷（正規化）

    %% 極座標網格（z=0 平面）：半徑 l_bar、方位 θ
    rr_um = 0:3:150;   thd = 0:5:360;                               % l_bar[µm]、θ[deg]
    [RR, TH] = meshgrid(rr_um, thd*pi/180);                        % size = [nth, nr]
    X = RR .* cos(TH);   Y = RR .* sin(TH);                        % µm（Cartesian 渲染）
    gain = zeros(size(RR));   iso = zeros(size(RR));
    for a = 1:numel(RR)
        r = RR(a)*1e-6;
        p  = r * [cos(TH(a)); sin(TH(a)); 0];                      % z=0 平面物理點 [m]（measure）
        Dk = p/ell_m - dhat_bias;                                  % 3×6（bias 電荷）
        sv = svd((Dk ./ (vecnorm(Dk).^3)) * Hhat_I);              % 3 個奇異值 [mT/A]
        gain(a) = norm(sv);   iso(a) = sv(1)/sv(3);
    end
    fprintf('bias 中心 r=0：gain=%.3f mT/A、iso=%.4f\n', gain(1,1), iso(1,1));

    %% ---- 兩張單一山丘 surf（gain z 標籤 Gain；iso 翻色階）------------------
    render_hill(X, Y, gain, 'Gain (mT/A)', fullfile(figdir,'svd_gain_3d.png'), false);
    render_hill(X, Y, iso,  'iso',       fullfile(figdir,'svd_iso_3d.png'),  true);     % iso 翻色階
end

function render_hill(X, Y, Z, zlab, outpng, flipcmap)
    if nargin < 6, flipcmap = false; end
    zmn = min(Z(:));  zmx = max(Z(:));
    fig = figure('Color','w','Position',[80 80 1000 880]);
    surf(X, Y, Z, 'EdgeColor','none');                            % CData = Z（顏色=gain/iso 大小）
    cmap = jet;  if flipcmap, cmap = flipud(cmap); end            % iso 翻色階
    shading interp; colormap(cmap); caxis([zmn zmx]);

    grid off; box on;                                             % 框體 = box on
    xlim([-150 150]); ylim([-150 150]); zlim([zmn zmx]);
    view(-40, 30);
    pbaspect([1 1 1]);                                            % 立方 plot box（軸異質 → pbaspect）
    set(gca,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    set(gca,'XTick',-150:75:150,'YTick',-150:75:150);
    zt = get(gca,'ZTick');                                        % MATLAB 自動 z 刻度
    tol = 0.055*(zmx - zmn);                                      % 丟掉太靠近端點的自動刻度，避免與 min/max 標籤重疊
    zt = unique([zmn, zt(zt>zmn+tol & zt<zmx-tol), zmx]);        % 保留內部自動刻度 + 多加 min/max
    set(gca,'ZTick', zt, 'ZTickLabel', arrayfun(@(v) sprintf('%.4g',v), zt, 'UniformOutput',false));
    xlabel('x (\mum)','FontWeight','bold');
    ylabel('y (\mum)','FontWeight','bold');
    zlabel(zlab,'FontWeight','bold');

    cb = colorbar; cb.FontSize = 16; cb.FontWeight = 'bold';
    cb.Ticks = cb.Ticks(1:2:end);
    cb.Label.String = zlab; cb.Label.FontWeight = 'bold'; cb.Label.FontSize = 16;

    ax = gca; ax.Toolbar.Visible = 'off';
    exportgraphics(fig, outpng, 'Resolution', 150);
    fprintf('saved %s\n', outpng);
    close(fig);
end
