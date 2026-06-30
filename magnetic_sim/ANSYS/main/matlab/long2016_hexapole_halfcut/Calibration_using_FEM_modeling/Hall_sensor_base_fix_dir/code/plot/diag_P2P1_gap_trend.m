function diag_P2P1_gap_trend()
% DIAG_P2P1_GAP_TREND  P1 激發下，P2 sensor 處 all-source B·n+ 隨 gap 的趨勢（+解釋圖）。
%   迴圈 gap∈{0,50,100,150,200}µm，在 P2 sensor 中心用 sensor_local 精細網格單點內插 B，
%   算 all-source B·n+（下極 P1 ×-1）+ B 與 n+ 夾角 → 證實「gap 越大、B 越倒向 n−」。
%   畫 2-panel 解釋圖：(A) 兩條磁通路徑機制示意；(B) B·n+ @P2 vs gap 量化曲線。
%   注意：sensor B·n+ 為「局部單點內插」估計（baseline 粗網格 sensor 圓柱內節點少）。

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath(['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
             'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir\code\function']);
    cnst = mt_constants();
    [sp, sn] = build_sensor_geometry(cnst);
    rr  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
    mcd = fullfile(rr,'mesh','standard','csv');
    figdir = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
              'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir\figures'];

    %% ---- 標準局部網格 → triangulation（恆用 baseline；gap 同 node_id）----
    N=readmatrix(fullfile(mcd,'sensor_local_nodes.csv'));  E=readmatrix(fullfile(mcd,'sensor_local_elems.csv'));
    nid=N(:,1); P=N(:,2:4); mxid=max(nid); g2l=zeros(mxid,1); g2l(nid)=1:numel(nid);
    Sl=E(:,2:9); tets=zeros(size(Sl,1),4); kk=0;
    for r=1:size(Sl,1), u=unique(Sl(r,:),'stable'); if numel(u)==4, kk=kk+1; tets(kk,:)=g2l(u); end, end
    tets=tets(1:kk,:);
    v1=P(tets(:,2),:)-P(tets(:,1),:); v2=P(tets(:,3),:)-P(tets(:,1),:); v3=P(tets(:,4),:)-P(tets(:,1),:);
    vol=dot(v1,cross(v2,v3,2),2); bad=vol<0; tets(bad,[3 4])=tets(bad,[4 3]);
    TR=triangulation(tets,P);
    ctr=(sp(:,2)+[0;0;cnst.SPH_OFST]).';  ni=sn(:,2);          % P2 sensor 中心(ANSYS frame) + n+
    ti=pointLocation(TR,ctr); bc=cartesianToBarycentric(TR,ti,ctr); conn=TR.ConnectivityList(ti,:);

    %% ---- 迴圈 gap：載 coil1/gapXXum_mueq、內插 B、算 all-source B·n+ + 夾角 ----
    gaps=[0 50 100 150 200];  Bn=zeros(size(gaps)); ang=zeros(size(gaps)); Bmag=zeros(size(gaps));
    muUp=[280 137.2 95.4 73.1 59.3]; muLo=[280 164.6 114.4 87.6 71.0];
    fprintf('\n gap(um)  mu_up/lo   B.n+ (all-src,T)   |B|(T)    angle(B,n+)\n');
    for k=1:numel(gaps)
        ds=import_ansys_data(fullfile(rr,'coil1',sprintf('gap%dum_mueq',gaps(k))),'all','coil1');
        m2=max(max(nid),max(ds.node_id)); id2=zeros(m2,1); id2(ds.node_id)=1:numel(ds.node_id);
        li=zeros(numel(nid),1); inb=nid<=m2; li(inb)=id2(nid(inb));
        Bnode=[ds.bx(li), ds.by(li), ds.bz(li)];
        Bp=bc*Bnode(conn,:);                 % 內插 B @ P2 sensor (raw FEM)
        Bas=-Bp;                             % all-source（P1 下極 ×-1）
        Bn(k)=Bas*ni; Bmag(k)=norm(Bas); ang(k)=acosd(max(-1,min(1,(Bas*ni)/norm(Bas))));
        fprintf(' %5d   %5.0f/%-5.0f   %+.4e      %.4e   %6.1f deg\n', gaps(k),muUp(k),muLo(k),Bn(k),Bmag(k),ang(k));
    end
    % 過零點（線性內插）
    zc=NaN; for k=1:numel(gaps)-1, if Bn(k)>0 && Bn(k+1)<=0, zc=gaps(k)+(0-Bn(k))/(Bn(k+1)-Bn(k))*(gaps(k+1)-gaps(k)); break; end, end
    fprintf('過零點 gap ≈ %.0f µm（B·n+ 由 +n+ 翻 −n+）\n\n', zc);

    %% ================= 解釋圖 =================
    fig=figure('Position',[60 60 1500 680],'Color','w');

    % ---------- Panel A：機制示意（horseshoe：P2 腳 + yoke + P1 腳）----------
    axA=subplot(1,2,1); hold(axA,'on'); axis(axA,'equal','off');
    title(axA,'機制：P1 磁通到 P2 尖端的兩條路徑','FontSize',14,'FontWeight','bold');
    gray=[.74 .74 .74]; LWfe=12; blue=[.1 .35 .9]; red=[.85 .12 .12];
    yL=[-3 3.4]; yR=[3 3.4]; P2t=[-0.55 0.15]; P1t=[0.55 0.15];      % yoke 兩端 + 兩尖端(近 WP)
    plot(axA,[yL(1) yR(1)],[yL(2) yR(2)],'-','Color',gray,'LineWidth',LWfe);      % yoke
    plot(axA,[yL(1) P2t(1)],[yL(2) P2t(2)],'-','Color',gray,'LineWidth',LWfe);    % P2 腳
    plot(axA,[yR(1) P1t(1)],[yR(2) P1t(2)],'-','Color',gray,'LineWidth',LWfe);    % P1 腳
    % 支撐座 gap 缺口（白）
    g2=[-1.75 1.95]; g1=[1.75 1.95];
    plot(axA,g2(1),g2(2),'s','MarkerSize',16,'MarkerFaceColor','w','MarkerEdgeColor',[.4 .4 .4],'LineWidth',1.2);
    plot(axA,g1(1),g1(2),'s','MarkerSize',16,'MarkerFaceColor','w','MarkerEdgeColor',[.4 .4 .4],'LineWidth',1.2);
    text(axA,g2(1)-1.05,g2(2),'gap','FontSize',10); text(axA,g1(1)+0.25,g1(2),'gap','FontSize',10);
    % coil on P1 腳
    text(axA,2.05,2.7,'coil','FontSize',11,'FontWeight','bold');
    % WP + 標記
    plot(axA,0,-0.45,'p','MarkerSize',20,'MarkerFaceColor',[1 .85 .1],'MarkerEdgeColor','k');
    text(axA,0.15,-0.75,'WP','FontSize',12,'FontWeight','bold');
    text(axA,P2t(1)-1.55,P2t(2)+0.05,'P2 (sensor)','FontSize',11,'FontWeight','bold','Color',[0 .4 0]);
    text(axA,P1t(1)+0.25,P1t(2)+0.0,'P1 (driven)','FontSize',11,'FontWeight','bold');
    % ① 藍：鐵內 flux 方向(P1腳↑ / yoke← / P2腳↓) + 出 P2 尖端(+n+)
    quiver(axA,1.95,1.0,(yR(1)-P1t(1))*0.22,(yR(2)-P1t(2))*0.22,0,'Color',blue,'LineWidth',2.4,'MaxHeadSize',3);
    quiver(axA,1.4,3.4,-1.3,0,0,'Color',blue,'LineWidth',2.4,'MaxHeadSize',1.4);
    quiver(axA,-1.95,2.55,(P2t(1)-yL(1))*0.22,(P2t(2)-yL(2))*0.22,0,'Color',blue,'LineWidth',2.4,'MaxHeadSize',3);
    npl=(P1t-P2t)/norm(P1t-P2t);                                     % +n+ 朝 WP 方向
    quiver(axA,P2t(1),P2t(2),npl(1)*0.85,npl(2)*0.85,0,'Color',blue,'LineWidth',3.2,'MaxHeadSize',2.5);
    text(axA,-2.95,1.35,'① 繞 yoke → P2 射出 +n_+','Color',blue,'FontSize',12,'FontWeight','bold','Interpreter','tex');
    % ② 紅：P1 尖端 → 直接灌進 P2 尖端(−n+)
    quiver(axA,P1t(1),P1t(2),(P2t(1)-P1t(1))*0.88,(P2t(2)-P1t(2))*0.88,0,'Color',red,'LineWidth',3.2,'MaxHeadSize',1.8);
    text(axA,-0.65,0.62,'② 直灌進 P2 −n_+','Color',red,'FontSize',12,'FontWeight','bold','Interpreter','tex');
    % 註解框
    text(axA,-3.4,-1.75,{'gap\uparrow → \mu_{eff}\downarrow → 藍路(支撐座)磁阻\uparrow', ...
        '→ 藍弱、紅(直灌)主導 → sensor B 倒向 n_-'}, ...
        'FontSize',11.5,'Color',[0 0 0],'Interpreter','tex','BackgroundColor',[1 1 .88],'EdgeColor',[.6 .6 .3]);
    xlim(axA,[-3.9 3.9]); ylim(axA,[-2.2 4]);

    % ---------- Panel B：B·n+ vs gap ----------
    axB=subplot(1,2,2); hold(axB,'on'); box(axB,'on');
    xl=[-8 208];
    fill(axB,[xl fliplr(xl)],[0 0 max(Bn)*1.25 max(Bn)*1.25],[.85 .9 1],'EdgeColor','none');   % +n+ 區
    fill(axB,[xl fliplr(xl)],[0 0 min(Bn)*1.25 min(Bn)*1.25],[1 .88 .88],'EdgeColor','none');  % -n+ 區
    plot(axB,xl,[0 0],'k-','LineWidth',1);
    plot(axB,gaps,Bn,'-o','Color',[.1 .1 .1],'LineWidth',2,'MarkerSize',8,'MarkerFaceColor',[.95 .75 .1]);
    if ~isnan(zc), plot(axB,zc,0,'kv','MarkerSize',11,'MarkerFaceColor','k'); text(axB,zc,max(Bn)*0.12,sprintf(' 過零 ~%.0fµm',zc),'FontSize',11,'FontWeight','bold'); end
    text(axB,150,max(Bn)*0.6,'+n_+（P2 射出 / 回流）','Color',[.1 .35 .9],'FontSize',12,'FontWeight','bold','Interpreter','tex');
    text(axB,150,min(Bn)*0.6,'-n_+（灌進 P2）','Color',[.85 .12 .12],'FontSize',12,'FontWeight','bold','Interpreter','tex');
    xlabel(axB,'gap (µm)','FontSize',13); ylabel(axB,'B·n_+ @ P2 sensor  (all-source, T)','FontSize',13,'Interpreter','tex');
    title(axB,'P2 sensor B·n_+ 隨 gap 翻向 n_-（局部內插估計）','FontSize',14,'FontWeight','bold','Interpreter','tex');
    set(axB,'FontSize',12); xlim(axB,xl); ylim(axB,[min(Bn)*1.25 max(Bn)*1.25]);

    if ~exist(figdir,'dir'); mkdir(figdir); end
    out=fullfile(figdir,'P2P1_signflip_gap_explain.png');
    exportgraphics(fig,out,'Resolution',150);
    fprintf('saved: %s\n', out);
end
