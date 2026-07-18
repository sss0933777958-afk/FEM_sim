function n_inside = sensing_area_count(varargin)
%SENSING_AREA_COUNT  外圓內可取的小圓中心點數（小圓邊緣不可超出外圓）。
% -------------------------------------------------------------------------
% 問題：外圓半徑 r_max，內含一半徑 r 的小圓；小圓可在外圓內平移，但小圓外邊
%   不可超過外圓。問外圓內可取幾個（小圓中心）點？（sensing area 取點）
%
% 演算法（手稿 step1–4）：
%   step1  小圓邊上點 (r·cosθ, r·sinθ)，θ = 0..360°（小圓上所有點）。
%   step2  位移量 δ = [δ·cosφ ; δ·sinφ]，方向 φ 掃整圈、角度間距 dphi=0.1°
%          （δ = 位移大小；θ 描述小圓邊、φ 描述位移方向，兩者不同）。
%   step3  新座標 [x;y] = [r·cosθ + δ·cosφ ; r·sinθ + δ·sinφ]（位移項用 φ，非 θ）。
%   step4  檢查小圓邊上「所有」點是否都在外圓內；掃完整圈 → 換下一個位移。
%   最後   計數 = 通過 step4（小圓完全在外圓內）的中心點總數。
%
% 說明：候選中心 = (δ·cosφ, δ·sinφ)。小圓完全在外圓內 ⇔ 中心到原點距離 δ ≤ r_max−r
%   （最外側邊點沿徑向 = δ+r）。故有效中心落在半徑 R_valid = r_max−r 的圓盤內。
%
% 參數（手稿只給 dphi=0.1°；徑向步長 ddelta 與 θ 取樣 dtheta 未指定 → 設可調）：
%   'r_max'  外圓半徑 [mm]   預設 2.865
%   'r'      小圓半徑 [mm]   預設 0.15
%   'dphi'   位移方向步長 [deg] 預設 0.1（手稿給定）
%   'ddelta' 位移大小步長 [mm]  預設 0.05（徑向取樣解析度）
%   'dtheta' 小圓邊檢查步長 [deg] 預設 1
%   'plot'   是否輸出圖 [logical] 預設 true → 存到 voltage_base/figures/shared/
% 輸出：n_inside = 外圓內可取的中心點數。
% -------------------------------------------------------------------------

    p = inputParser;
    p.addParameter('r_max', 2.865);
    p.addParameter('r',     0.15);
    p.addParameter('dphi',  0.1);
    p.addParameter('ddelta',0.05);
    p.addParameter('dtheta',1);
    p.addParameter('plot',  true);
    p.parse(varargin{:});
    o = p.Results;
    r_max = o.r_max;  r = o.r;

    % ---- 取樣網格 --------------------------------------------------------
    theta = (0:o.dtheta:360-o.dtheta) * pi/180;   % step1：小圓邊參數 θ
    phi   = (0:o.dphi  :360-o.dphi )  * pi/180;    % step2：位移方向 φ（間距 0.1°）
    delta = 0:o.ddelta:r_max;                       % 位移大小 δ（掃到 r_max）

    ct = cos(theta);  st = sin(theta);              % 小圓邊單位向量

    % ---- step2–4：逐位移檢查小圓是否完全在外圓內，計數有效中心 ----------
    n_inside = 0;
    valid_c  = [];                                  % 保留有效中心（畫圖用，粗取樣）
    keep_every = max(1, round(2/o.dphi));           % 圖上每 ~2° 取一個中心，避免過密
    for d = delta
        for k = 1:numel(phi)
            cx = d*cos(phi(k));  cy = d*sin(phi(k));            % step3：新中心
            ex = cx + r*ct;       ey = cy + r*st;               % 小圓邊上所有點
            if all(ex.^2 + ey.^2 <= r_max^2)                    % step4：全部在外圓內
                n_inside = n_inside + 1;
                if mod(k, keep_every)==0, valid_c(end+1,:) = [cx cy]; end %#ok<AGROW>
            end
            if d==0, break; end                                 % δ=0 中心唯一，不重複掃 φ
        end
    end

    R_valid = r_max - r;
    fprintf('r_max=%.4g mm, r=%.4g mm, dphi=%.3g deg, ddelta=%.3g mm\n', r_max, r, o.dphi, o.ddelta);
    fprintf('有效中心區 = 半徑 R_valid = r_max-r = %.4g mm 的圓盤\n', R_valid);
    fprintf('外圓內可取點數 n_inside = %d\n', n_inside);

    if ~o.plot, return; end

    % ---- 圖：外圓 + 有效中心區 + 有效中心取樣 + 示例小圓 -----------------
    here    = fileparts(mfilename('fullpath'));                 % .../code/function
    outdir  = fullfile(fileparts(fileparts(here)), 'figures', 'shared');
    if ~exist(outdir,'dir'); mkdir(outdir); end
    tt = linspace(0,2*pi,400);

    fig = figure('Color','w','Position',[100 100 720 720]);
    ax  = axes(fig); hold(ax,'on');
    % 外圓（r_max）
    plot(ax, r_max*cos(tt), r_max*sin(tt), 'k-', 'LineWidth', 2.5);
    % 有效中心區邊界（r_max - r）
    plot(ax, R_valid*cos(tt), R_valid*sin(tt), '--', 'Color',[0.2 0.5 0.9], 'LineWidth', 1.8);
    % 有效中心取樣點
    if ~isempty(valid_c)
        plot(ax, valid_c(:,1), valid_c(:,2), '.', 'Color',[0.85 0.33 0.10], 'MarkerSize', 4);
    end

    axis(ax,'equal'); box(ax,'on'); grid(ax,'off');
    lim = r_max*1.08;  xlim(ax,[-lim lim]); ylim(ax,[-lim lim]);
    set(ax,'FontSize',14,'FontWeight','bold','LineWidth',1.5);
    xlabel(ax,'x (mm)','FontWeight','bold'); ylabel(ax,'y (mm)','FontWeight','bold');
    legend(ax, {sprintf('外圓 r_{max}=%.3g mm', r_max), ...
                sprintf('有效中心區 r_{max}-r=%.3g mm', R_valid), ...
                '有效中心（取樣）'}, ...
           'Location','southoutside','FontSize',11);
    text(ax, -lim*0.95, lim*0.9, sprintf('n_{inside} = %d\n(dφ=%.2g°, dδ=%.3g mm)', n_inside, o.dphi, o.ddelta), ...
         'FontSize',13,'FontWeight','bold','BackgroundColor','w','EdgeColor','k','Margin',5);

    outpng = fullfile(outdir, 'sensing_area_count.png');
    exportgraphics(fig, outpng, 'Resolution', 150);
    fprintf('已輸出 %s\n', outpng);
    close(fig);
end
