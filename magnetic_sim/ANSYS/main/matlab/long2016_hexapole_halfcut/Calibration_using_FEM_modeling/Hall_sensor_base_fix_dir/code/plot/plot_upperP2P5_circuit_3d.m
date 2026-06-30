function plot_upperP2P5_circuit_3d(EXC, VARIANT, SHOW_FIELD)
% PLOT_UPPERP2P5_CIRCUIT_3D  選定磁極區 3D 磁路箭頭圖（真實 FEM 節點、不內插）。
% -------------------------------------------------------------------------
%   真實 FEM 節點 B 箭頭（長度∝|B|、上限 cap；顏色=|B| log turbo）+ 選定磁極輪廓
%   （per-pole 軸鄰近濾波，上/下極皆可）+ 各極 sensor 取樣圓柱(沿 n+)。**無 WP/極尖點標記。**
%   兩張圖（迴圈）：Fig1 P2/P5/P3 @ az=210；Fig2 P2/P1/P5 @ az=30。
%   EXC    : 激發極 paper 名（預設 'P2'）→ 載對應 coil；all-source（上極不翻、下極 sink 翻）。
%   VARIANT: FEM 變體（預設 'gap200um_mueq'）。 SHOW_FIELD: 預設 true（畫場箭頭）。
% -------------------------------------------------------------------------
    if nargin < 1 || isempty(EXC),        EXC        = 'P2';            end
    if nargin < 2 || isempty(VARIANT),    VARIANT    = 'gap200um_mueq'; end
    if nargin < 3 || isempty(SHOW_FIELD), SHOW_FIELD = true;            end
    DPI = 200;

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath(['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
             'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir\code\function']);
    cnst = mt_constants();
    [sp, sn] = build_sensor_geometry(cnst);                 % 3×6 sensor 中心/法線（WP 框 [m]）
    rr = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';

    % ---- 激發 coil；all-source：只在「激發極是下極 sink」時翻號 ----
    apdl_to_paper_idx = [1,3,6,5,2,4];
    pname = {'P1','P2','P3','P4','P5','P6'};
    pe = find(strcmp(pname, EXC));  kc = find(apdl_to_paper_idx == pe);
    is_lower_exc = ismember(pe, [1 3 6]);
    d = import_ansys_data(fullfile(rr,sprintf('coil%d',kc),VARIANT),'all',sprintf('coil%d',kc));
    if is_lower_exc, d.bx=-d.bx; d.by=-d.by; d.bz=-d.bz; end
    fprintf('載 coil%d (%s 激發)/%s：%d 節點，|B|max=%.4f T\n', kc, EXC, VARIANT, numel(d.x), max(sqrt(d.bx.^2+d.by.^2+d.bz.^2)));

    airsel = filter_iron_nodes(d.x,d.y,d.z,cnst,struct('visualize',false));
    isair  = false(numel(d.x),1); isair(airsel)=true;
    X = [d.x, d.y, d.z-cnst.SPH_OFST]*1e3;  B = [d.bx,d.by,d.bz];      % WP 框 [mm]
    rxy=cnst.R_norm_xy*1e3; rz=cnst.R_norm_z*1e3; Rnorm=sqrt(rxy^2+rz^2);

    base = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
            'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir\figures\'];

    % ---- 兩張圖設定 ----
    CFG = { struct('poles',[2 5 3],'view',[210 36.59],'name','circuit_3d_P2P3P5_az210'), ...
            struct('poles',[2 1 5],'view',[ 30 36.59],'name','circuit_3d_P1P2P5_az30') };
    for ci = 1:numel(CFG)
        draw_one(CFG{ci}, X, B, isair, cnst, sp, sn, rxy, rz, Rnorm, SHOW_FIELD, base, DPI, pname);
    end
end

%% ===================== 單張圖 =====================
function draw_one(cfg, X, B, isair, cnst, sp, sn, rxy, rz, Rnorm, SHOW_FIELD, base, DPI, pname)
    poles = cfg.poles;
    % 各選定極：極尖(上極+rz/下極−rz) + sensor 中心/法線（WP 框 mm）
    np = numel(poles);
    tip = zeros(np,3); scen = zeros(np,3); nrm = zeros(np,3);
    for ii=1:np
        p = poles(ii); th = cnst.pole_angles(p)*pi/180;
        sgn = 1; if cnst.pole_is_lower(p), sgn = -1; end
        tip(ii,:)  = [rxy*cos(th), rxy*sin(th), sgn*rz];
        scen(ii,:) = sp(:,p).'*1e3;  nrm(ii,:) = sn(:,p).';
    end

    % box：所選極 tip + 沿極軸外延一段(框入更多錐身) + sensor（含沿 n+ 一小段）+ pad
    tipout = tip + 1.5*(tip./vecnorm(tip,2,2));            % 各極沿自身軸往外 1.5mm
    key = [tip; tipout; scen; scen+0.7*nrm];  pad=0.8; CELL=0.30;
    bx=[min(key(:,1))-pad, max(key(:,1))+pad];
    by=[min(key(:,2))-pad, max(key(:,2))+pad];
    bz=[min(key(:,3))-pad, max(key(:,3))+pad];
    inbox = X(:,1)>=bx(1)&X(:,1)<=bx(2)&X(:,2)>=by(1)&X(:,2)<=by(2)&X(:,3)>=bz(1)&X(:,3)<=bz(2);

    % 場取點（盒內空氣，格點抽樣）
    PP=X(inbox&isair,:); BB=B(inbox&isair,:);
    [~,iu]=unique(round(PP/CELL),'rows','stable'); PP=PP(iu,:); BB=BB(iu,:);
    Bmag=vecnorm(BB,2,2);

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
    fig=figure('Position',[60 60 1180 900],'Color','w');
    ax=axes(fig,'Position',[0.06 0.06 0.86 0.86]); hold(ax,'on');

    if size(ironP,1) >= 4
        try
            shp=alphaShape(ironP(:,1),ironP(:,2),ironP(:,3),0.9); hs=plot(shp);
            if SHOW_FIELD, set(hs,'FaceColor',[.62 .66 .72],'FaceAlpha',0.30,'EdgeColor','none');
            else,         set(hs,'FaceColor',[.62 .66 .72],'FaceAlpha',0.55,'EdgeColor',[.35 .38 .42],'EdgeAlpha',0.15); end
        catch
        end
    end

    if SHOW_FIELD && ~isempty(Bmag)
        bs=sort(Bmag); p95=bs(max(1,round(0.95*numel(bs)))); SCALE=0.6/p95;
        Vq=BB*SCALE; LMAX=0.8; Lq=vecnorm(Vq,2,2); over=Lq>LMAX; Vq(over,:)=Vq(over,:).*(LMAX./Lq(over));
        Cv=log10(max(Bmag,1e-9)); nb=20; ed=linspace(min(Cv),max(Cv),nb+1); cmap=turbo(nb);
        for k=1:nb
            in=Cv>=ed(k)&Cv<ed(k+1); if k==nb, in=in|(Cv>=ed(end)); end
            if any(in)
                quiver3(PP(in,1),PP(in,2),PP(in,3), Vq(in,1),Vq(in,2),Vq(in,3), 0, ...
                        'Color',cmap(k,:),'LineWidth',1.0,'MaxHeadSize',0.5);
            end
        end
        colormap(ax,turbo); clim([min(Cv) max(Cv)]); cb=colorbar(ax);
        lo=floor(min(Cv)); hi=ceil(max(Cv)); tk=lo:hi;
        cb.Ticks=tk; cb.TickLabels=arrayfun(@(t)sprintf('10^{%d}',t),tk,'UniformOutput',false);
        ylabel(cb,'|B| (3D) [Tesla]'); cb.FontSize=13; cb.FontWeight='bold';
    end

    % sensor 圓柱 + n+ + 標籤（每選定極；無 WP、無極尖點）
    R_SENS=0.15; H_SENS=0.10;
    for ii=1:np
        c=scen(ii,:); nh=nrm(ii,:);
        draw_sensor_cyl(ax, c, nh, R_SENS, H_SENS, [.1 .8 .25]);
        text(c(1),c(2),c(3),sprintf('  %s sensor',pname{poles(ii)}),'FontSize',12,'FontWeight','bold','Color',[.05 .5 .15]);
        quiver3(c(1),c(2),c(3),nh(1),nh(2),nh(3),0.7,'Color',[.9 .1 .1],'LineWidth',3,'MaxHeadSize',0.9);
    end

    hold(ax,'off'); grid(ax,'off'); box(ax,'on'); daspect(ax,[1 1 1]);
    xlim(bx); ylim(by); zlim(bz); view(ax,cfg.view(1),cfg.view(2));
    set(ax,'FontSize',13,'FontWeight','bold','LineWidth',1.5);
    xlabel('x (mm)','FontWeight','bold'); ylabel('y (mm)','FontWeight','bold'); zlabel('z (mm)','FontWeight','bold');
    ax.Toolbar.Visible='off';

    out=[base cfg.name '.png']; exportgraphics(fig,out,'Resolution',DPI);
    fprintf('saved: %s\n', out); close(fig);
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
