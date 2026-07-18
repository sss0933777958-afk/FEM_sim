function ntu_draw_plate(x, y, zbot, ztop, col, alpha)
%NTU_DRAW_PLATE  畫 NTU 扁平磁極板：上下面 patch + 側壁 surf + 上下邊線（半透明、無面邊線）。
%   x,y  : Nx1 輪廓（世界座標，mm；通常 = ntu_pole_profile 旋轉後）
%   zbot/ztop : 板底 / 板頂 z（mm）
%   col  : 1x3 顏色（邊線用 col*0.6）；alpha : 面透明度（P1_full 風格 = 0.30）
%
%   [ADDED 2026-07-17] 由 plot/plot_poles_sensors_3d.m 與 plot/plot_p1_sensing_disk.m
%   的重複 local function `draw_plate` 抽出共用。行為與原版逐字相同。
    n=numel(x);
    patch(x,y,zbot*ones(n,1),col,'FaceAlpha',alpha,'EdgeColor','none');
    patch(x,y,ztop*ones(n,1),col,'FaceAlpha',alpha,'EdgeColor','none');
    surf([[x;x(1)] [x;x(1)]],[[y;y(1)] [y;y(1)]],[repmat(zbot,n+1,1) repmat(ztop,n+1,1)], ...
        'FaceColor',col,'FaceAlpha',alpha,'EdgeColor','none');
    plot3([x;x(1)],[y;y(1)],ztop*ones(n+1,1),'-','Color',col*0.6,'LineWidth',1.0);
    plot3([x;x(1)],[y;y(1)],zbot*ones(n+1,1),'-','Color',col*0.6,'LineWidth',1.0);
end
