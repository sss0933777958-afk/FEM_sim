function plot_p1_circuit_3d()
% plot_p1_circuit_3d — NTU P1 整條磁路的 3D B 箭頭圖：含襯套 vs 無襯套（coil1 激發）。
%   目的：解釋「加襯套後 P1 sensor 讀值 495→219 mV（−55.7%），WP |B| 卻 +6.4%」。
%   機制：P1 sensor 距襯套軸 7.416mm < flange 外徑 7.5mm → sensor 正在凸緣外緣正下方 6.14mm。
%
%   資料：ANSYS_data/NTU_hexapole/data/<variant>/coil1（**真實 FEM 節點、不內插**，per plot-real-nodes）。
%   箭頭：air 節點 → 格點抽稀(CELL) → 長度 ∝ |B|（p95 定尺度 + 硬 cap）、AutoScale=0；
%         顏色 log10|B| 分 NB bin turbo；**兩張共用同一 clim**（per figure-style「同類比較圖共用色階」）。
%   幾何輪廓（風格照 figures/shared/P1_full.png）：P1 板（共用 helper ntu_pole_profile/ntu_draw_plate）
%         + 導柱 + 襯套（僅 sleeve 版）+ yoke（只畫邊框，否則屋頂擋住視線）+ sensing area 紅圓柱。
%   出 2 檔 → figures/shared/p1_circuit_3d_{sleeve,nosleeve}.png
%   風格：粗體框選項① + view(-30,40) + daspect([1 1 1])。

    % ======================== 參數（要調就改這裡）========================
    CELL   = 0.40;          % [MODIFIED] 箭頭抽稀格距 [mm]（0.70→0.40：使用者要更密。資料一律真實節點原值、不內插）
    SLAB_Y = 1.20;          % 只取 |y| <= SLAB_Y 的薄片（含極軸 y=0 與 sensor y=-0.99）
    % [MODIFIED] 長度改 **log 壓縮**：|B| 動態範圍 55×（p50=4 / p95=222 mT），線性 ∝|B| 會讓
    %   中位數箭頭只有 0.007mm 長＝看不見。改成長度隨 log10|B| 線性內插到 [LMIN,LMAX]
    %   → 仍是「依強度給長度」，但弱場也看得見；絕對值另由顏色（log10 色階）表達。
    LMIN   = 0.30*CELL;     % 最弱（clim 下限）的箭頭長度 [mm]
    LMAX   = 1.10*CELL;     % 最強（clim 上限）的箭頭長度 [mm]（「適中、不要太長」）
    NB     = 20;            % 顏色 bin 數（log10|B|）
    % [MODIFIED] 稍微 zoom-in 到 sensor 附近：x 從 4 起（切掉尖端外的空白）、z 到 21.5
    %   （仍含板 8.80-9.05、導柱、襯套 flange 15.60-16.60 + tube →20.80、yoke 底面 20.80）
    % [MODIFIED] 左端延伸到 -2.62（原 4）→ 涵蓋極尖端 x=0.353 與 WP x=0，看得到尖端射出的場。
    %   ⚠ x<0 區域在 |y|<SLAB_Y 薄片內會侵入其他極的板（P5@60°、P4@300°）→ 左端的鋼件箭頭不屬於 P1。
    % [MODIFIED] z 上緣 21.5 → 27：yoke 真實範圍是 z[20.80, 26.30]（厚 5.5mm），
    %   zlim=21.5 只露出底部 0.7mm、看起來像薄板（使用者指出「尺寸不對」）。27 才完整含 yoke 頂 26.30。
    BOX    = struct('x',[-2.62 30], 'y',[-8 8], 'z',[8.5 27]);  % P1 磁路盒
    VIEW   = [0 20];        % [MODIFIED] 使用者指定：方位角 0、仰角 20（≈ x–z 側視，最適合看磁路）
    SHOW_IRON = true;       % [MODIFIED] 鋼件內也畫箭頭 —— 磁通主要走鋼裡，濾掉就看不到磁路
                            %   （原本 false：只畫空氣。改 true 後鋼內 |B| 高會撐開色階，屬預期）
    % ====================================================================

    here   = fileparts(mfilename('fullpath'));
    vbase  = fileparts(fileparts(here));                       % .../voltage_base
    figdir = fullfile(vbase,'figures','shared');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    addpath(fullfile(vbase,'code','function'));                % ntu_pole_profile / ntu_draw_plate
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');   % import_ansys_data
    droot = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\NTU_hexapole\data';

    variants = {'full_assembly_sleeve','full_assembly'};
    tags     = {'sleeve','nosleeve'};
    hasSlv   = [true false];

    % ---- Pass 1：兩個 variant 都載入 + 抽稀，並算共用色階 ----
    S = cell(1,2);  allC = [];
    for v = 1:2
        S{v} = prep(fullfile(droot,variants{v},'coil1'), 'coil1', BOX, SLAB_Y, CELL, hasSlv(v), SHOW_IRON);
        allC = [allC; log10(max(S{v}.Bmag,1e-9))];             %#ok<AGROW>
        fprintf('%-22s: %7d 節點在盒內 → 抽稀後 %5d 箭頭 | |B| p50=%.3g p95=%.3g max=%.3g mT\n', ...
                variants{v}, S{v}.nAir, size(S{v}.P,1), median(S{v}.Bmag), prctile(S{v}.Bmag,95), max(S{v}.Bmag));
    end
    CLIM = [prctile(allC,2) prctile(allC,98)];                 % 共用 clim（log10 mT）
    fprintf('共用 clim(log10 mT) = [%.3f %.3f]  → |B| %.3g ~ %.3g mT\n', CLIM, 10^CLIM(1), 10^CLIM(2));

    % ---- Pass 2：逐張渲染（同 clim、同長度對映 → 兩張可直接比）----
    for v = 1:2
        render(S{v}, hasSlv(v), CLIM, LMIN, LMAX, NB, BOX, VIEW, ...
               fullfile(figdir, sprintf('p1_circuit_3d_%s.png', tags{v})));
    end
end

% ============================================================================
function s = prep(dir_in, prefix, BOX, SLAB_Y, CELL, hasSleeve, showIron)
% 載入 → mm/mT → 盒內 + y 薄片 →（可選）去鋼件 → 格點抽稀（每格第一個，比照 long2016 藍本）
    d  = import_ansys_data(dir_in, 'all', prefix);
    P  = [d.x d.y d.z]*1e3;                    % m → mm
    B  = [d.bx d.by d.bz]*1e3;                 % T → mT
    in = P(:,1)>=BOX.x(1) & P(:,1)<=BOX.x(2) & abs(P(:,2))<=SLAB_Y & ...
         P(:,3)>=BOX.z(1) & P(:,3)<=BOX.z(2);
    P=P(in,:); B=B(in,:);
    if ~showIron                               % [MODIFIED] 預設 showIron=true：鋼內也畫（磁通主要走鋼裡）
        air = ~is_iron(P, hasSleeve);
        P=P(air,:); B=B(air,:);
    end
    s.nAir = size(P,1);
    [~,iu] = unique(round(P/CELL),'rows','stable');
    s.P = P(iu,:);  s.B = B(iu,:);  s.Bmag = vecnorm(s.B,2,2);
end

function m = is_iron(P, hasSleeve)
% NTU 鋼件遮罩（mm，world）：6 板 + 6 導柱 + 6 襯套 + yoke。只用來把箭頭排除在鋼外。
    x=P(:,1); y=P(:,2); z=P(:,3);  m=false(size(x));
    prof = ntu_pole_profile();
    ang  = [0 120 240 60 180 300];             % 下 0/120/240（板底 8.80）、上 60/180/300（板底 9.55）
    for k = 1:6
        t=ang(k)*pi/180;  zb = 8.80 + 0.75*(k>3);          % 上層板底 9.55
        xl =  x*cos(t) + y*sin(t);   yl = -x*sin(t) + y*cos(t);   % world → pole-local
        m = m | (inpolygon(xl,yl,prof(:,1),prof(:,2)) & z>=zb & z<=zb+0.25);
    end
    for k = 1:6                                 % 導柱 r2.5 @ r=20.5
        t=(k-1)*60*pi/180;  cx=20.5*cos(t); cy=20.5*sin(t);
        zlo = 9.05; if mod(k-1,2)==1, zlo = 9.80; end        % 60/180/300° 為短柱
        rr = hypot(x-cx,y-cy);
        m = m | (rr<2.5 & z>=zlo & z<=26.80);
        if hasSleeve
            m = m | (rr<7.5 & z>=15.60 & z<=16.60);          % flange
            m = m | (rr<3.5 & z>=16.60 & z<=20.80);          % tube
        end
    end
    % [MODIFIED] yoke = 外盤 − cutout（geom deck：YAX=[29.9,15.5]、YAY=[26.9,10.25]、VSBV 外盤−cutout）
    %   → 是個「框」，中央 |x|<15.5 & |y|<10.25 是洞。原本寫成整片實心 = 錯（把洞當成鋼）。
    inOuter = abs(x)<=29.9 & abs(y)<=26.9;
    inHole  = abs(x)<15.5 & abs(y)<10.25;
    m = m | (inOuter & ~inHole & z>=20.80 & z<=26.30);
end

% ============================================================================
function render(s, hasSleeve, CLIM, LMIN, LMAX, NB, BOX, VIEW, outpng)
    colP=[0.30 0.55 0.95]; colPost=[0.55 0.58 0.62]; colSlv=[0.20 0.65 0.35];
    colYk=[0.35 0.38 0.42]; colS=[0.85 0.10 0.10];

    f = figure('Color','w','Position',[100 100 1400 1000]); hold on;

    % ---- 幾何輪廓（P1_full 風格：半透明面 + col*0.6 邊線）----
    B2 = ntu_pole_profile();                                   % P1 @0° 不旋轉
    ntu_draw_plate(B2(:,1), B2(:,2), 8.80, 9.05, colP, 0.12);  % P1 板（alpha 壓低：鋼內現在有箭頭）
    draw_tube(20.5, 0, 2.5, 9.05, 26.80, colPost, 0.12);       % P1 導柱（長柱）
    if hasSleeve
        draw_tube(20.5, 0, 7.5, 15.60, 16.60, colSlv, 0.16);   % 襯套 flange
        draw_tube(20.5, 0, 3.5, 16.60, 20.80, colSlv, 0.16);   % 襯套 tube
    end
    draw_yoke_slab(BOX, 29.9, 20.80, 26.30, colYk);            % yoke：盒內的板（外緣 x=29.9），淡淡帶過
    draw_cyl_solid(13.15, -0.99, 9.46, 0.15, 0.10, colS);      % sensing area
    text(13.15,-0.99,9.46+1.3,'sensing area','FontSize',14,'FontWeight','bold','Color',colS, ...
         'HorizontalAlignment','center');
    % [ADDED] 元件標名（讓圖自明）
    text(20.5, 0, 25.2,'post','FontSize',13,'FontWeight','bold','Color',colPost*0.7,'HorizontalAlignment','center');
    text(6, 0, 22.0,'yoke','FontSize',13,'FontWeight','bold','Color',colYk,'HorizontalAlignment','center');
    text(2.0, 0, 9.6,'P1 pole','FontSize',13,'FontWeight','bold','Color',colP*0.7,'HorizontalAlignment','center');
    if hasSleeve
        text(20.5, 0, 14.2,'sleeve','FontSize',13,'FontWeight','bold','Color',colSlv*0.8,'HorizontalAlignment','center');
    end

    % ---- 箭頭：長度隨 log10|B| 內插到 [LMIN,LMAX]（共用對映）；顏色 log10|B| 分 bin（共用 clim）----
    Cv = log10(max(s.Bmag,1e-9));
    w  = min(max((Cv-CLIM(1))/(CLIM(2)-CLIM(1)), 0), 1);       % clim 外夾住
    Lw = LMIN + (LMAX-LMIN).*w;                                % 每根箭頭的目標長度 [mm]
    V  = s.B ./ max(s.Bmag,eps) .* Lw;                         % 方向 × 目標長度
    cmap = turbo(NB);  ed = linspace(CLIM(1),CLIM(2),NB+1);
    for k = 1:NB
        if k==1,      sel = Cv<=ed(2);
        elseif k==NB, sel = Cv>ed(NB);
        else,         sel = Cv>ed(k) & Cv<=ed(k+1);
        end
        if ~any(sel); continue; end
        quiver3(s.P(sel,1),s.P(sel,2),s.P(sel,3), V(sel,1),V(sel,2),V(sel,3), 0, ...
                'Color',cmap(k,:),'LineWidth',1.0,'MaxHeadSize',0.5);
    end
    colormap(gca,turbo(NB));  clim(CLIM);
    cb = colorbar;  cb.FontSize=14; cb.FontWeight='bold';
    dt = ceil(CLIM(1)):floor(CLIM(2));
    if numel(dt)>=2
        cb.Ticks = dt;  cb.TickLabels = arrayfun(@(d) sprintf('10^{%d}',d), dt, 'uni',0);
    end
    cb.Label.String = '|B| (mT)';  cb.Label.FontWeight='bold'; cb.Label.FontSize=15;

    % ---- 風格：粗體框選項① ----
    ax = gca;
    xlim(BOX.x); ylim(BOX.y); zlim(BOX.z);
    grid(ax,'off'); box(ax,'on'); daspect(ax,[1 1 1]); view(ax,VIEW(1),VIEW(2));
    set(ax,'FontSize',15,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    xlabel(ax,'x (mm)','FontWeight','bold'); ylabel(ax,'y (mm)','FontWeight','bold');
    zlabel(ax,'z (mm)','FontWeight','bold');
    xt=get(ax,'XTick'); if numel(xt)>3; set(ax,'XTick',xt(1:2:end)); end
    yt=get(ax,'YTick'); if numel(yt)>3; set(ax,'YTick',yt(1:2:end)); end
    ax.Toolbar.Visible='off';
    exportgraphics(f, outpng, 'Resolution', 150);  close(f);
    fprintf('saved %s\n', outpng);
end

% ============================================================================
function draw_tube(cx, cy, R, z0, z1, col, alpha)
% 半透明圓柱（側壁 + 上下蓋）。bore 被導柱填滿 → 畫實心柱在視覺上等價。
    [Xc,Yc,Zc] = cylinder(R,36);  Xc=Xc+cx; Yc=Yc+cy;  Zc = Zc*(z1-z0)+z0;
    surf(Xc,Yc,Zc,'FaceColor',col,'FaceAlpha',alpha,'EdgeColor','none');
    th = linspace(0,2*pi,37);
    for zc = [z0 z1]
        patch(cx+R*cos(th), cy+R*sin(th), zc*ones(size(th)), col, ...
              'FaceAlpha',alpha,'EdgeColor','none');
        plot3(cx+R*cos(th), cy+R*sin(th), zc*ones(size(th)), '-','Color',col*0.6,'LineWidth',1.0);
    end
end

function draw_yoke_slab(BOX, xout, z0, z1, col)
% yoke = 外盤(半寬 X 29.9 / Y 26.9) − cutout(半寬 X 15.5 / Y 10.25) → 是個**框**，中央是洞。
% [MODIFIED ×2] ①改封閉六面盒（原只有上下兩面 + 一片側壁 → 看起來沒填滿）；
%   ②**x 從 15.5 起**（原從 BOX.x(1)=-2.62 起 = 把 cutout 的洞畫成實心，錯）。
%   本圖 y 範圍 [-8,8] 全在 cutout 的 |y|<10.25 內 → 此視野下 yoke 只存在於 x >= 15.5。
%   alpha 與導柱 draw_tube 一致（0.12）+ col*0.6 邊線。
    XHOLE = 15.5;                                   % cutout 半寬 X（YAX(2)）
    x0 = max(BOX.x(1), XHOLE);  x1 = min(BOX.x(2), xout);
    if x1 <= x0, return; end
    y0 = BOX.y(1);  y1 = BOX.y(2);
    a  = 0.12;
    F = { [x0 x1 x1 x0; y0 y0 y1 y1; z0 z0 z0 z0], ...      % 底
          [x0 x1 x1 x0; y0 y0 y1 y1; z1 z1 z1 z1], ...      % 頂
          [x0 x1 x1 x0; y0 y0 y0 y0; z0 z0 z1 z1], ...      % 前 (y=y0)
          [x0 x1 x1 x0; y1 y1 y1 y1; z0 z0 z1 z1], ...      % 後 (y=y1)
          [x1 x1 x1 x1; y0 y1 y1 y0; z0 z0 z1 z1], ...      % 外緣 (x=29.9)
          [x0 x0 x0 x0; y0 y1 y1 y0; z0 z0 z1 z1] };        % 內側 (x=x0)
    for k = 1:numel(F)
        patch(F{k}(1,:), F{k}(2,:), F{k}(3,:), col, ...
              'FaceAlpha',a, 'EdgeColor',col*0.6, 'LineWidth',1.0);
    end
end

% [MODIFIED 2026-07-17] sensor 圓柱改畫成「底面在 cz、往 +n 長 H」，與 build_V_matrix 的取樣一致
%   （a = axial_tol*rand → pts = c + a*n̂，即 c 是**底面中心**、不是幾何中心）。
%   原本畫成以 cz 為中心（cz±H/2）→ 與實際取樣位置差 H/2 = 0.05mm。
function draw_cyl_solid(cx, cy, cz, R, H, col)
    [Xc,Yc,Zc]=cylinder(R,24); Xc=Xc+cx; Yc=Yc+cy; Zc=Zc*H+cz;
    surf(Xc,Yc,Zc,'FaceColor',col,'EdgeColor','none','FaceAlpha',1.0);
    th=linspace(0,2*pi,25);
    for zc=[cz, cz+H]
        patch(cx+R*cos(th),cy+R*sin(th),zc*ones(size(th)),col,'EdgeColor','none');
    end
end
