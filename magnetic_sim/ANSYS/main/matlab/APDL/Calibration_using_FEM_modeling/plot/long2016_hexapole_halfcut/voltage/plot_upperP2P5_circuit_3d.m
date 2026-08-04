function plot_upperP2P5_circuit_3d(EXC, VARIANT, SHOW_FIELD)
% PLOT_UPPERP2P5_CIRCUIT_3D  選定磁極區 3D 磁路箭頭圖（真實 FEM 節點、不內插）。
% -------------------------------------------------------------------------
%   真實 FEM 節點 B 箭頭（長度∝|B|、上限 cap；顏色=|B| log turbo）+ 選定磁極輪廓
%   （per-pole 軸鄰近濾波，上/下極皆可）+ 各極 sensor 取樣圓柱(沿 n+)。**無 WP/極尖點標記。**
%   兩張圖（迴圈）：Fig1 P2/P5/P3 @ az=210；Fig2 P2/P1/P5 @ az=30。
%   EXC    : 激發極 paper 名（預設 'P2'）→ 載對應 coil；all-source（上極不翻、下極 sink 翻）。
%   VARIANT: FEM 變體（預設 'gap_200um'）。 SHOW_FIELD: 預設 true（畫場箭頭）。
% -------------------------------------------------------------------------
    if nargin < 1 || isempty(EXC),        EXC        = 'P2';     end
    if nargin < 2 || isempty(VARIANT),    VARIANT    = 'graded'; end   % [MODIFIED] paper 圖用龍飛 graded
    if nargin < 3 || isempty(SHOW_FIELD), SHOW_FIELD = true;     end
    DPI = 200;

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');   % mt_constants
    % [MODIFIED] 原本 addpath 舊 per-model 樹 voltage_base\code\function（已刪除）→ 改掛共用 function/
    addpath(['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\' ...
             'Calibration_using_FEM_modeling\function']);           % import_ansys_data / filter_iron_nodes
    cnst = mt_constants();
    [sp, sn] = sensor_geometry_local(cnst);                 % 3×6 sensor 中心/法線（WP 框 [m]）
                                                            % [MODIFIED] 原 build_sensor_geometry 隨舊樹刪除，
                                                            %   改用照抄 build_V_matrix.m 之 canonical local（見檔末）
    rr = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';

    % ---- 激發 coil；all-source：只在「激發極是下極 sink」時翻號 ----
    apdl_to_paper_idx = [1,3,6,5,2,4];
    pname = {'P1','P2','P3','P4','P5','P6'};
    pe = find(strcmp(pname, EXC));  kc = find(apdl_to_paper_idx == pe);
    is_lower_exc = ismember(pe, [1 3 6]);
    d = import_ansys_data(fullfile(rr, VARIANT, sprintf('coil%d',kc)),'all',sprintf('coil%d',kc));
    if is_lower_exc, d.bx=-d.bx; d.by=-d.by; d.bz=-d.bz; end
    fprintf('載 coil%d (%s 激發)/%s：%d 節點，|B|max=%.4f T\n', kc, EXC, VARIANT, numel(d.x), max(sqrt(d.bx.^2+d.by.^2+d.bz.^2)));

    airsel = filter_iron_nodes(d.x,d.y,d.z,cnst,struct('visualize',false));
    isair  = false(numel(d.x),1); isair(airsel)=true;
    X = [d.x, d.y, d.z-cnst.SPH_OFST]*1e3;  B = [d.bx,d.by,d.bz];      % WP 框 [mm]
    rxy=cnst.R_norm_xy*1e3; rz=cnst.R_norm_z*1e3; Rnorm=sqrt(rxy^2+rz^2);

    % [MODIFIED] 落點改 paper 圖夾（原 voltage_base\figures\shared 已隨舊樹刪除）
    base = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\paper_fig\Section3_A\';
    if ~exist(base,'dir'), mkdir(base); end

    % ---- paper 圖設定（使用者拍板 2026-07-31：P2+P5 兩極、az=30/el=30、graded、激發 P2）----
    %   舊的三極兩張（P2P3P5 az210 / P1P2P5 az30）輸出到已刪除的舊樹，此處不再產生。
    CFG = { struct('poles',[2 5],'view',[30 30],'name','upperP2P5_circuit_3d') };
    for ci = 1:numel(CFG)
        draw_one(CFG{ci}, X, B, isair, cnst, sp, sn, rxy, rz, Rnorm, SHOW_FIELD, base, DPI, pname);
    end
end

%% ===================== 單張圖 =====================
function draw_one(cfg, X, B, isair, cnst, sp, sn, rxy, rz, Rnorm, SHOW_FIELD, base, DPI, pname)
    poles = cfg.poles;
    % 各選定極：極尖(上極+rz/下極−rz) + sensor 中心/法線（WP 框 mm）
    L_POLE = 4.5;                              % [ADDED] 解析錐體沿極軸畫多長 [mm]（也用於「極內」判準）
    np = numel(poles);
    tip = zeros(np,3); scen = zeros(np,3); nrm = zeros(np,3); pax = zeros(np,3);
    for ii=1:np
        p = poles(ii); th = cnst.pole_angles(p)*pi/180;
        sgn = 1; if cnst.pole_is_lower(p), sgn = -1; end
        tip(ii,:)  = [rxy*cos(th), rxy*sin(th), sgn*rz];
        scen(ii,:) = sp(:,p).'*1e3;  nrm(ii,:) = sn(:,p).';
        if cnst.pole_is_lower(p)                                  % [ADDED] 極軸（尖端→根部），供錐體與選點共用
            pax(ii,:) = [cos(th), sin(th), 0];
        else
            inc = cnst.upper_incline;
            pax(ii,:) = [cos(inc)*cos(th), cos(inc)*sin(th), sin(inc)];
        end
    end

    % box：所選極 tip + 沿極軸外延一段(框入更多錐身) + sensor（含沿 n+ 一小段）+ pad
    tipout = tip + 1.5*(tip./vecnorm(tip,2,2));            % 各極沿自身軸往外 1.5mm
    key = [tip; tipout; scen; scen+0.7*nrm];  pad=0.8; CELL=0.30;
    bx=[min(key(:,1))-pad, max(key(:,1))+pad];
    by=[min(key(:,2))-pad, max(key(:,2))+pad];
    bz=[min(key(:,3))-pad, max(key(:,3))+pad];
    inbox = X(:,1)>=bx(1)&X(:,1)<=bx(2)&X(:,2)>=by(1)&X(:,2)<=by(2)&X(:,3)>=bz(1)&X(:,3)<=bz(2);

    % 場取點 = 盒內空氣節點（原始取樣，格點抽樣）
    %   [MODIFIED 2026-07-31] 曾試加入「磁極內部鐵節點」，但鐵內 |B| 高一個數量級 → 色階上限被拉到
    %   450mT、空氣箭頭全被壓成深藍看不見（使用者決定：只畫空氣）。
    sel = inbox & isair;
    PP=X(sel,:); BB=B(sel,:);
    [~,iu]=unique(round(PP/CELL),'rows','stable'); PP=PP(iu,:); BB=BB(iu,:);
    Bmag=vecnorm(BB,2,2);
    fprintf('  取點：空氣 %d → 抽樣後 %d\n', nnz(sel), size(PP,1));

    % 結構：盒內鐵節點，per-pole 軸鄰近 union（上/下極通用，無 z-floor）
    ironP = X(inbox & ~isair, :);
    RKEEP=3.5; AMIN=0.80*Rnorm; keepP=false(size(ironP,1),1);
    for ii=1:np
        u = tip(ii,:).'/norm(tip(ii,:)); a = ironP*u; perp = vecnorm(ironP - a*u.',2,2);
        keepP = keepP | (perp<RKEEP & a>AMIN);
    end
    ironP = ironP(keepP,:);
    fprintf('[%s] 盒內空氣 %d、pole-only iron %d（極 %s）\n', cfg.name, size(PP,1), size(ironP,1), mat2str(poles));

    %% ---- 圖 ----
    FS = 36;                                   % [MODIFIED] paper 圖字體統一 36 粗體（figure-style.md）
    fig=figure('Position',[40 40 1240 1180],'Color','w');
    ax=axes(fig,'Position',[0.09 0.06 0.67 0.90]); hold(ax,'on');   % 左留 z 刻度數字、右留 colorbar + 36pt 標題

    % [MODIFIED] 磁極輪廓改「解析錐體」（圓角尖端 + 切線直錐 + 圓柱段），只畫選定的極（P2/P5）。
    %   原為 FEM 鐵節點 alphaShape → 形狀糊、且會把鄰極/軛鐵的節點一起包進來。
    %   幾何與 figures/paper_fig_plot/plot_p2_pole_full.m 的 draw_pole_full 同一套。
    rf   = cnst.POLE_TIP_R*1e3;                       % 尖端圓角半徑 [mm]
    beta = atan2(cnst.POLE_R, cnst.POLE_CONE_LEN);    % 半錐角 ≈ 11.31°
    Rcyl = cnst.POLE_R*1e3;                           % 錐→柱轉換半徑 [mm]
    for ii = 1:np
        p = poles(ii);  th = cnst.pole_angles(p)*pi/180;
        if cnst.pole_is_lower(p)
            a = [cos(th); sin(th); 0];                                     % 下極：水平軸
        else
            inc = cnst.upper_incline;                                      % 上極：真實錐軸傾角 ≈ 36.59°
            a = [cos(inc)*cos(th); cos(inc)*sin(th); sin(inc)];
        end
        u = cross(a,[0;0;1]);  u = u/norm(u);  v = cross(a,u);             % 面內 ⊥ 基底
        draw_pole_full(ax, tip(ii,:).', a, u, v, rf, beta, Rcyl, L_POLE, [0.72 0.74 0.78], 0.45);
    end

    if SHOW_FIELD && ~isempty(Bmag)
        bs=sort(Bmag); p95=bs(max(1,round(0.95*numel(bs)))); SCALE=0.6/p95;
        Vq=BB*SCALE; LMAX=0.8; Lq=vecnorm(Vq,2,2); over=Lq>LMAX; Vq(over,:)=Vq(over,:).*(LMAX./Lq(over));
        % [MODIFIED] 色階改「|B| 場 colorbar 標準樣式」(figure-style.md)：turbo + 線性 clim[0 CLIM]、
        %   CLIM 進位到 50mT、28 個 bin 逐 bin 上色；單位 Tesla → mT。(原為 log10 + Tesla)
        Bm_mT = Bmag*1e3;
        CLIM  = ceil(max(Bm_mT)/50)*50;
        nb=28; ed=linspace(0,CLIM,nb+1); cmap=turbo(nb);
        for k=1:nb
            in=Bm_mT>=ed(k)&Bm_mT<ed(k+1); if k==nb, in=in|(Bm_mT>=ed(end)); end
            if any(in)
                quiver3(PP(in,1),PP(in,2),PP(in,3), Vq(in,1),Vq(in,2),Vq(in,3), 0, ...
                        'Color',cmap(k,:),'LineWidth',1.0,'MaxHeadSize',0.5);
            end
        end
        colormap(ax,turbo); clim(ax,[0 CLIM]); cb=colorbar(ax);
        style_cbar(cb, FS);
    end

    % sensor 圓柱 + n+ 箭頭（保留實體；[MODIFIED] 使用者要求拿掉「P2 sensor / P5 sensor」文字標籤）
    R_SENS=0.15; H_SENS=0.10;
    for ii=1:np
        c=scen(ii,:); nh=nrm(ii,:);
        draw_sensor_cyl(ax, c, nh, R_SENS, H_SENS, [.1 .8 .25]);
        quiver3(c(1),c(2),c(3),nh(1),nh(2),nh(3),0.7,'Color',[.9 .1 .1],'LineWidth',3,'MaxHeadSize',0.9);
    end

    % [MODIFIED] 風格照 figure-style.md：3D 變體 A（x,y,z 同單位同尺度 → daspect + 手動框邊、
    %   省「離相機最近角」相連 3 邊）；字體 36 粗體；**移除三軸標題**（使用者要求）；刻度不含端點、奇數個。
    grid(ax,'off'); box(ax,'off'); daspect(ax,[1 1 1]);
    xlim(bx); ylim(by); zlim(bz); view(ax,cfg.view(1),cfg.view(2));
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',3.0,'TickLength',[.018 .018]);
    set(ax,'XTick',ticks3(bx),'YTick',ticks3(by),'ZTick',ticks3(bz));
    draw_box_edges3(ax, bx, by, bz, 3.0);          % 手動框邊（黑、等粗、省最近角 3 邊）
    hold(ax,'off');                                % ⚠ 必須在框邊畫完之後才 hold off，
    ax.Toolbar.Visible='off';                      %   否則 plot3 會逐次清空座標軸內容

    out=[base cfg.name '.png']; exportgraphics(fig,out,'Resolution',DPI);
    fprintf('saved: %s\n', out); close(fig);
end

%% ---- local：解析磁極整根（圓角尖端 + 切線直錐 + 圓柱段）；單位 mm ----
%   照抄 figures/paper_fig_plot/plot_p2_pole_full.m 的 draw_pole_full；改動要兩邊同步。
function draw_pole_full(ax, tip, a, u, v, rf, beta, Rcyl, L, col, alp)
    ts = rf*(1+sin(beta));  rt = rf*cos(beta);
    L_cone = ts + (Rcyl - rt)/tan(beta);            % 錐→柱轉換軸向位置
    psi = linspace(0, pi/2+beta, 16);               % 圓角尖端
    axoff = rf*(1-cos(psi));   radm = rf*sin(psi);
    Lc = min(L, L_cone);                            % 直錐
    tc = linspace(ts, Lc, 44);
    axoff = [axoff, tc(2:end)];   radm = [radm, rt + (tc(2:end)-ts)*tan(beta)];
    if L > L_cone                                   % 圓柱段
        tcyl = linspace(L_cone, L, 20);
        axoff = [axoff, tcyl(2:end)];   radm = [radm, Rcyl*ones(1,numel(tcyl)-1)];
    end
    phi = linspace(0, 2*pi, 72);
    N = numel(axoff);  X=zeros(N,numel(phi)); Y=X; Z=X;
    for j = 1:numel(phi)
        rad = u*cos(phi(j)) + v*sin(phi(j));   P = tip + a*axoff + rad*radm;
        X(:,j)=P(1,:).'; Y(:,j)=P(2,:).'; Z(:,j)=P(3,:).';
    end
    surf(ax, X,Y,Z, 'FaceColor',col, 'FaceAlpha',alp, 'EdgeColor','none');
    rim = tip + a*axoff(end) + (u*cos(phi)+v*sin(phi))*radm(end);  cen = tip + a*axoff(end);
    patch(ax, [cen(1) rim(1,:)],[cen(2) rim(2,:)],[cen(3) rim(3,:)], col, 'EdgeColor','none','FaceAlpha',alp);
end

%% ---- local：|B| colorbar 標準樣式（figure-style.md；同 plot_circuit_side.m 的 style_cbar）----
function style_cbar(cb, FS)
    cb.FontSize = FS;  cb.FontWeight = 'bold';
    cb.Label.Interpreter = 'latex';             % ⚠ 必須先設 Interpreter 再設 String，
    cb.Label.String = '$\mathbf{|B|\;(mT)}$';   %   否則預設 tex 解譯不懂 \mathbf 會報錯
    cb.Label.FontSize = FS;
    cb.Units = 'normalized';
    CBW_RATIO = 0.009;                          % [MODIFIED] 照 figure-style.md 定案：cbw/figW = 0.009（細長條）
    cb.Position = [0.800 0.16 CBW_RATIO 0.68];  % 手動定位：緊靠圖框、標題不被裁
end

%% ---- local：3D 刻度（3 個等距、不含端點；figure-style「角落不放 tick、奇數個」）----
function tk = ticks3(lim)
    lo = lim(1);  hi = lim(2);  mid = (lo+hi)/2;
    s  = nice_step((hi-lo)/5);                  % 等距 nice 步長
    tk = round(mid/s)*s + [-1 0 1]*s;           % 3 個等距、對稱於中心
    tk = tk(tk > lo + 0.06*(hi-lo) & tk < hi - 0.06*(hi-lo));   % 不含端點/角落
end

function s = nice_step(x)
    k = floor(log10(x));  m = x/10^k;
    cand = [1 2 2.5 5 10];
    [~,i] = min(abs(cand - m));
    s = cand(i)*10^k;
end

%% ---- local：手動框邊（黑、等粗；省「離相機最近角」相連的 3 邊）----
function draw_box_edges3(ax, bx, by, bz, lw)
    C = [bx([1 1 1 1 2 2 2 2]).', by([1 1 2 2 1 1 2 2]).', bz([1 2 1 2 1 2 1 2]).'];   % 8 角
    E = [1 2;3 4;5 6;7 8; 1 3;2 4;5 7;6 8; 1 5;2 6;3 7;4 8];                            % 12 邊
    cp = campos(ax);
    [~, near] = min(sum((C - cp).^2, 2));       % 離相機最近的角
    for e = 1:size(E,1)
        if any(E(e,:) == near), continue; end   % 省掉與最近角相連的 3 邊
        plot3(ax, C(E(e,:),1), C(E(e,:),2), C(E(e,:),3), 'k-', 'LineWidth', lw);
    end
end

%% ---- local：sensor 幾何（照抄 function/build_V_matrix.m 的 canonical local sensor_geometry）----
%   ⚠ 與 build_V_matrix.m 為同一份定義，改動要兩邊同步。
function [sensor_pos, sensor_n] = sensor_geometry_local(cfg)
    SOFF_upper = 4.572e-3;
    beta   = atan2(cfg.POLE_R, cfg.POLE_CONE_LEN);    % 半錐角 ≈ 11.31°
    psi0   = atan2(cfg.R_norm_z, cfg.R_norm_xy);      % 仰角 ≈ 35.26°（magic-angle）
    inc_up = cfg.upper_incline;                       % 上極真實錐軸傾角 ≈ 36.59°
    ell    = cfg.R_norm;
    SOFF_lower = 4.572e-3;  AIR = 0.41e-3;
    dir = @(el,az) [cos(el)*cos(az); cos(el)*sin(az); sin(el)];
    sensor_pos = zeros(3,6);  sensor_n = zeros(3,6);
    for i = 1:6
        th = cfg.pole_angles(i)*pi/180;
        if cfg.pole_is_lower(i), psi = -psi0; else, psi = +psi0; end
        e1 = dir(psi, th);
        if cfg.pole_is_lower(i)
            e2 = dir(-beta, th);   nhat = dir(-beta-pi/2, th);   soff = SOFF_lower;
        else
            e2 = dir(inc_up+beta, th);  nhat = dir(inc_up+beta+pi/2, th);  soff = SOFF_upper;
        end
        sensor_pos(:,i) = ell*e1 + soff*e2 + AIR*nhat;
        sensor_n(:,i)   = nhat;
    end
end

%% ---- local：sensor 取樣圓柱（軸=n̂、半徑 R、高 H）----
function draw_sensor_cyl(ax, c, nhat, R, H, col)
    nhat=nhat(:)/norm(nhat); c=c(:).';
    t1=[-nhat(2);nhat(1);0]; if norm(t1)<1e-9, t1=[1;0;0]; end
    t1=t1/norm(t1); t2=cross(nhat,t1);
    th=linspace(0,2*pi,28).';
    ring = R*(cos(th)*t1.' + sin(th)*t2.');
    bot = c + ring;  top = c + H*nhat.' + ring;
    surf(ax, [bot(:,1) top(:,1)], [bot(:,2) top(:,2)], [bot(:,3) top(:,3)], 'FaceColor',col,'FaceAlpha',0.92,'EdgeColor','none');
    fill3(ax, bot(:,1),bot(:,2),bot(:,3), col,'FaceAlpha',0.92,'EdgeColor','none');
    fill3(ax, top(:,1),top(:,2),top(:,3), col,'FaceAlpha',0.92,'EdgeColor','none');
end
