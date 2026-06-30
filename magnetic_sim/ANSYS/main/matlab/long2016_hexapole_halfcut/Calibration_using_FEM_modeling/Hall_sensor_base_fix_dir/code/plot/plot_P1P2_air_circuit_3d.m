function plot_P1P2_air_circuit_3d(PREVIEW, FOCUS, VIEW, VARIANT, DATASET, FLIP)
% PLOT_P1P2_AIR_CIRCUIT_3D  baseline(完美接觸) P1 激發、all-source 的 3D 磁路箭頭圖。
% -------------------------------------------------------------------------
%   真實 FEM 節點 B 箭頭（不內插，all-source 翻 P1）。箭頭=單位方向、顏色=|B|(log)。
%   FOCUS：
%     'sensor'(預設)：P2 sensor + 鄰近 cone flank（匯進↔漏出轉換）。只空氣節點。
%     'full'        ：P1+P2 對極近 WP 全景。只空氣節點。
%     'p2pole'      ：**P2 整根磁極 + 支撐座(protrusion) + yoke 接點**（WP→錐→座→yoke）。
%                     含**鋼體內磁通**（支撐座的磁路在鋼裡）→ 看 flux 經支撐座灌進 yoke（完美接觸的捷徑）。
%   VIEW=[az el]（預設依 FOCUS）。PREVIEW=true → 暫存。
% -------------------------------------------------------------------------
    if nargin < 1 || isempty(PREVIEW), PREVIEW = false;   end
    if nargin < 2 || isempty(FOCUS),   FOCUS   = 'sensor'; end
    if nargin < 3 || isempty(VIEW)
        switch FOCUS
            case 'sensor', VIEW=[-32 16];
            case 'p2pole', VIEW=[-12 12];
            otherwise,     VIEW=[-25 12];
        end
    end
    if nargin < 4 || isempty(VARIANT), VARIANT = 'standard'; end   % FEM 變體子夾（'standard'/'graded_p2'…）
    if nargin < 5 || isempty(DATASET), DATASET = 'all';      end   % 匯出 dataset（'all'/'p2reg'…）
    if nargin < 6 || isempty(FLIP),    FLIP    = 'flip';     end   % 'flip'=all-source 翻 P1 / 'none'=不翻（資料已 all-source）
    DPI = 200;

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath(['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
             'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir\code\function']);
    cnst = mt_constants();
    [sp, sn] = build_sensor_geometry(cnst);
    rr = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';

    % ---- 載 coil1（P1 激發）→ 視 FLIP 決定是否 all-source 翻號 ----
    d = import_ansys_data(fullfile(rr,'coil1',VARIANT),DATASET,'coil1');
    if strcmp(FLIP,'flip'), d.bx=-d.bx; d.by=-d.by; d.bz=-d.bz; end  % all-source（P1 下極）；'none'=資料已 all-source(如 graded CURR=-1)
    fprintf('載入 coil1/%s [%s]：%d 節點，|B|max=%.4f T（FLIP=%s）\n', ...
            VARIANT, DATASET, numel(d.x), max(sqrt(d.bx.^2+d.by.^2+d.bz.^2)), FLIP);

    % ---- air / iron 分類 ----
    airsel = filter_iron_nodes(d.x,d.y,d.z,cnst,struct('visualize',false));
    isair  = false(numel(d.x),1); isair(airsel)=true;

    % ---- WP 框 [mm] ----
    X = [d.x, d.y, d.z-cnst.SPH_OFST]*1e3;  B = [d.bx,d.by,d.bz];

    % ---- 關鍵點（WP 框 mm）----
    th1=cnst.pole_angles(1)*pi/180; th2=cnst.pole_angles(2)*pi/180;
    rxy=cnst.R_norm_xy*1e3; rz=cnst.R_norm_z*1e3;
    P1tip=[rxy*cos(th1), rxy*sin(th1), -rz];
    P2tip=[rxy*cos(th2), rxy*sin(th2), +rz];
    scm=(sp(:,2).')*1e3;  nmm=sn(:,2).';
    Pprot=[-47.5, 0, 5.5];                                        % P2 支撐座(protrusion)中心(azimuth180) [mm]

    incw=cnst.upper_incline; betaw=atan2(3.0,15.0);
    rotw=@(v,a)[cos(a)*v(1)-sin(a)*v(2); sin(a)*v(1)+cos(a)*v(2)];
    fdir2=rotw([-cos(incw);sin(incw)],betaw);  fdir=[fdir2(1);0;fdir2(2)];
    Pflip = P2tip(:) + 3.0*fdir;

    % ---- 盒 + 取點策略 ----
    if strcmp(FOCUS,'sensor')
        Pa=P2tip(:)+1.5*fdir; Pb=P2tip(:)+5.5*fdir;
        key=[Pa.'; Pb.'; scm; scm+0.7*nmm]; padxz=0.7; pady=1.3; CELL=0.18; ARROW=0.17; useAll=false;
        bx=[min(key(:,1))-padxz, max(key(:,1))+padxz]; by=[-pady,pady]; bz=[min(key(:,3))-padxz, max(key(:,3))+padxz];
    elseif strcmp(FOCUS,'p2pole')
        % P2 整根極 + 支撐座 + yoke 接點（azimuth 180、y≈0 平面）：固定盒涵蓋整根 + 薄 y 片
        %   graded_p2 資料已細解析到 yoke(~-47mm)；baseline 只到 ~-16mm（盒大但後段空）。
        bx=[-52, 4]; by=[-3, 3]; bz=[-3, 19];
        CELL=1.2; ARROW=1.6; useAll=true;                       % 含鋼體內磁通（薄片去 3D 雜訊）
    else
        key=[P1tip; P2tip; 0 0 0; scm]; padxz=0.9; pady=2.0; CELL=0.32; ARROW=0.26; useAll=false;
        bx=[min(key(:,1))-padxz, max(key(:,1))+padxz]; by=[-pady,pady]; bz=[min(key(:,3))-padxz, max(key(:,3))+padxz];
    end
    inbox = X(:,1)>=bx(1)&X(:,1)<=bx(2)&X(:,2)>=by(1)&X(:,2)<=by(2)&X(:,3)>=bz(1)&X(:,3)<=bz(2);

    % ---- 取節點（p2pole=全節點含鋼，其餘=只空氣）+ 格點抽樣 ----
    sel = inbox; if ~useAll, sel = inbox & isair; end
    PP=X(sel,:); BB=B(sel,:);
    [~,iu]=unique(round(PP/CELL),'rows','stable');  PP=PP(iu,:); BB=BB(iu,:);
    Bmag=vecnorm(BB,2,2); uvec=BB./max(Bmag,eps);
    ironP = X(inbox & ~isair, :);
    fprintf('盒內：取點 %d（useAll=%d）；iron 節點 %d\n', size(PP,1), useAll, size(ironP,1));

    %% ---- 圖 ----
    fig=figure('Position',[60 60 1180 900],'Color','w');
    ax=axes(fig,'Position',[0.08 0.08 0.85 0.82]); hold(ax,'on');

    % 輪廓：sensor/full 用 alphaShape；p2pole 範圍大改用淡 iron 點雲襯底
    if ~useAll && size(ironP,1) >= 4
        try
            shp=alphaShape(ironP(:,1),ironP(:,2),ironP(:,3),0.5);
            hs=plot(shp); set(hs,'FaceColor',[.55 .55 .58],'FaceAlpha',0.13,'EdgeColor','none');
        catch
            scatter3(ironP(:,1),ironP(:,2),ironP(:,3),3,[.6 .6 .6],'filled','MarkerFaceAlpha',0.2);
        end
    end

    % 箭頭：log 色階分箱
    Cv=log10(max(Bmag,1e-9)); nb=20; ed=linspace(min(Cv),max(Cv),nb+1); cmap=turbo(nb);
    for k=1:nb
        in=Cv>=ed(k)&Cv<ed(k+1); if k==nb, in=in|(Cv>=ed(end)); end
        if any(in)
            quiver3(PP(in,1),PP(in,2),PP(in,3), uvec(in,1)*ARROW,uvec(in,2)*ARROW,uvec(in,3)*ARROW, 0, ...
                    'Color',cmap(k,:),'LineWidth',1.0,'MaxHeadSize',0.6);
        end
    end
    colormap(ax,turbo); clim([min(Cv) max(Cv)]); cb=colorbar(ax);
    lo=floor(min(Cv)); hi=ceil(max(Cv)); tk=lo:hi;
    cb.Ticks=tk; cb.TickLabels=arrayfun(@(t)sprintf('10^{%d}',t),tk,'UniformOutput',false);
    ylabel(cb,'|B| (3D) [Tesla]');

    % 標記
    plot3(0,0,0,'p','MarkerSize',16,'MarkerFaceColor',[1 .84 0],'MarkerEdgeColor','k','LineWidth',1);
    text(0,0,0,'  WP','FontSize',11,'FontWeight','bold');
    plot3(P2tip(1),P2tip(2),P2tip(3),'o','MarkerSize',9,'MarkerFaceColor',[.1 .5 .95],'MarkerEdgeColor','k');
    text(P2tip(1),P2tip(2),P2tip(3),'  P2 tip','FontSize',11,'FontWeight','bold');
    if strcmp(FOCUS,'p2pole')
        plot3(Pprot(1),Pprot(2),Pprot(3),'s','MarkerSize',13,'MarkerFaceColor',[.95 .55 .1],'MarkerEdgeColor','k','LineWidth',1);
        text(Pprot(1),Pprot(2),Pprot(3),'  支撐座 protrusion→yoke','FontSize',12,'FontWeight','bold','Color',[.6 .3 0]);
    end
    if strcmp(FOCUS,'full')
        plot3(P1tip(1),P1tip(2),P1tip(3),'o','MarkerSize',9,'MarkerFaceColor',[.95 .4 .1],'MarkerEdgeColor','k');
        text(P1tip(1),P1tip(2),P1tip(3),'  P1 tip','FontSize',11,'FontWeight','bold');
    end
    if ~useAll
        plot3(scm(1),scm(2),scm(3),'p','MarkerSize',20,'MarkerFaceColor',[.1 .8 .25],'MarkerEdgeColor','k','LineWidth',1.2);
        text(scm(1),scm(2),scm(3),'  P2 sensor','FontSize',12,'FontWeight','bold','Color',[.05 .5 .15]);
        quiver3(scm(1),scm(2),scm(3),nmm(1),nmm(2),nmm(3),0.6,'Color',[.9 .1 .1],'LineWidth',3,'MaxHeadSize',0.9);
    end

    hold(ax,'off'); grid(ax,'on'); box(ax,'on'); daspect(ax,[1 1 1]);
    xlim(bx); ylim(by); zlim(bz); view(ax,VIEW(1),VIEW(2));
    xlabel('x [mm] (WP frame)'); ylabel('y [mm] (WP frame)'); zlabel('z [mm] (WP frame)');
    if strcmp(FOCUS,'p2pole')
        title({'P2 整根磁極 + 支撐座 3D 磁路 ｜ P1 激發 (all-source) ｜ baseline 完美接觸 ｜ 真實 FEM 節點(含鋼體內磁通)', ...
               '看磁通 WP→錐→支撐座→yoke：完美接觸下大量 flux 經支撐座灌回 yoke（捷徑）'},'Interpreter','none');
    else
        title({'P2 sensor 區 3D 空氣磁路 ｜ P1 激發 (all-source) ｜ baseline ｜ 真實 FEM 節點(未內插)', ...
               '近尖端匯進 P2，過了~3mm 轉成漏出/回 yoke，sensor 落在漏出段→朝外'},'Interpreter','none');
    end
    ax.Toolbar.Visible='off';

    %% ---- 輸出 ----
    base = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
            'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir\figures\'];
    switch FOCUS
        case 'sensor', nm='P2sensor_air_circuit_3d'; defv=[-32 16];
        case 'p2pole', nm='P2pole_circuit_3d';        defv=[-12 12];
        otherwise,     nm='P1P2_air_circuit_3d';      defv=[-25 12];
    end
    vtag=''; if ~strcmp(VARIANT,'standard'), vtag=['_' VARIANT]; end   % graded_p2 等變體加後綴
    suf=vtag; if ~isequal(VIEW, defv), suf=[vtag sprintf('_az%d_el%d',VIEW(1),VIEW(2))]; end
    if PREVIEW, out=fullfile(tempdir,[nm '_preview.png']); else, out=[base nm suf '.png']; end
    exportgraphics(fig,out,'Resolution',DPI);
    fprintf('saved: %s\n', out);
end
