function plot_circuit_side(USE_BIAS)
% plot_circuit_side -- long2016 半切六極 P1|P2 磁路箭頭「側視合併圖」(x-z, y=0)+ 電荷位置
% =========================================================================
%   左右合併:左 P1(coil1)、右 P2(coil5),側視(x-z, y=0)畫尖端區磁路箭頭 + 粉色電荷點。
%   raw 座標(WP 在 raw z=-12.711mm);source 號誌 s 讓場從尖端射出(全 source)。
%   場:'all' dataset、每格取最近 y=0 真實節點(不內插)、turbo 依 |B|(mT) 疊 quiver。
%   電荷:single = l_hat·dhat；eighteen = l_hat·R_actᵀ·(Pc_base+E(e))(讀 R150 .mat)。
%   **P1/P2 共用單一 colorbar**(同 flux 合併圖:shared clim=max,弱場 panel 顯冷色);兩 panel 共用 box 高度 H → 邊框對齊。
%   colorbar 寬度 = CBW_RATIO·figW(定案 0.009,與 flux 圖共用同值)。
%   USE_BIAS=false → circuit_side_merged.png、true → circuit_side_merged_bias.png。
%   nargin<1 → 兩張都產。輸出 → figures/paper_fig/Section2_E/(覆蓋迭代)。
% =========================================================================
    if nargin < 1
        plot_circuit_side(false);  plot_circuit_side(true);  return;
    end
    clc;
    here   = fileparts(mfilename('fullpath'));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CAL,'function'));  addpath(fullfile(CAL,'common_path'));
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    c = mt_constants();
    FS = 30;
    bstr = ''; if USE_BIAS, bstr = '_bias'; end
    % ---- 佈局（manual pixel；兩 panel 共用 box 高度 H → 邊框對齊；各自 colorbar）----
    H=680; y0=130; leftm=110; cbgap=22; cblab=140; midgap=120; rightm=30;
    CBW_RATIO = 0.009;                                  % colorbar 寬度佔比(cbw/figW)，定案；與 flux 共用同值

    S1 = load_panel(1, USE_BIAS, c, CAL);
    S2 = load_panel(2, USE_BIAS, c, CAL);
    CMAX = max(S1.CLIM, S2.CLIM);                       % P1/P2 共用色階(同 flux 合併圖;弱場 panel 顯冷色)
    w1 = H*diff(S1.XL)/diff(S1.ZL);   w2 = H*diff(S2.XL)/diff(S2.ZL);
    base = leftm + w1 + midgap + w2 + cbgap + cblab + rightm;
    cbw  = CBW_RATIO*base/(1-CBW_RATIO);               % 單一 colorbar，解 cbw/figW = ratio
    figW = base + cbw;   figH = y0 + H + 60;
    x1 = leftm;   x2 = x1 + w1 + midgap;

    fig = figure('Color','w','Units','pixels','Position',[20 40 figW figH]);
    ax1 = axes(fig,'Units','pixels');  render_panel_into(ax1, S1, CMAX, FS);  ax1.Position=[x1 y0 w1 H];
    ax2 = axes(fig,'Units','pixels');  render_panel_into(ax2, S2, CMAX, FS);  ax2.Position=[x2 y0 w2 H];
    cb  = colorbar(ax2,'Units','pixels');  cb.Position=[x2+w2+cbgap y0 cbw H];  ax2.Position=[x2 y0 w2 H];  style_cbar(cb, FS);

    outp = fullfile(figdir, sprintf('circuit_side_merged%s.png', bstr));
    exportgraphics(fig, outp, 'Resolution', 150);
    fprintf('wrote %s\n', outp);
end

% ============================================================================
function [XL,ZL,XT,ZT] = win(pidx)
% per-pole 視野 + 刻度（WP frame；z 已平移）。x 起終點都標；z 只留內縮 tick。
    switch pidx
        case 1, XL=[0.2 0.9]; ZL=[-0.6 0];  XT=[0.2 0.4 0.6 0.9];    ZT=[-0.5 -0.3 -0.1];
        case 2, XL=[-1.2 0];  ZL=[0 1];     XT=[-1.2 -0.8 -0.4 0];   ZT=[0.25 0.5 0.75];
    end
end

% ============================================================================
function S = load_panel(pidx, USE_BIAS, c, CAL)
% 載入單極場 + grid-sample(不內插) + 電荷位置 + 磁極輪廓；回傳繪圖用結構(含自然 CLIM)。
    switch pidx, case 1, cname='coil1'; case 2, cname='coil5'; end
    [XL,ZL,XT,ZT] = win(pidx);

    % ---- 載入場(graded / 'all')----
    ddir = ansys_path('long2016_hexapole_halfcut','data','graded',cname);
    d = import_ansys_data(ddir, 'all', cname);
    fprintf('P%d = %s (graded/all): matched=%d, |B|max=%.4f T\n', pidx, cname, numel(d.x), max(d.bsum));
    s = 1 - 2*c.pole_is_lower(pidx);                    % P1 -1 / P2 +1（全 source）

    % ---- y=0 薄板 + grid-sample 最近 y=0 真實節點(不內插)----
    zoff = -c.SPH_OFST*1e3;                             % raw z → WP frame（縱軸平移到 WP）
    x=d.x*1e3; y=d.y*1e3; z=d.z*1e3+zoff;               % mm；z 已平移到 WP frame
    slab=100; cell=0.016;                               % y 投影全取(不限薄板) / grid 格邊 ~0.016mm
    nx=round(diff(XL)/cell); nz=round(diff(ZL)/cell);
    inwin = x>=XL(1)&x<=XL(2) & z>=ZL(1)&z<=ZL(2) & abs(y)<slab;
    gi=find(inwin);
    xe=linspace(XL(1),XL(2),nx+1); ze=linspace(ZL(1),ZL(2),nz+1);
    ix=discretize(x(gi),xe); iz=discretize(z(gi),ze);
    ok=~isnan(ix)&~isnan(iz); gi=gi(ok); ix=ix(ok); iz=iz(ok);
    cid=(ix-1)*nz+iz; ay=abs(y(gi));
    [uc,~,g]=unique(cid); sel=zeros(numel(uc),1);
    for k=1:numel(uc), idx=find(g==k); [~,b]=min(ay(idx)); sel(k)=gi(idx(b)); end
    keep=d.bsum(sel)>1e-4; sel=sel(keep);

    Xs=x(sel); Zs=z(sel);
    Bx=s*d.bx(sel); Bz=s*d.bz(sel); Bm=d.bsum(sel);

    % ---- 箭頭長度(|B|^0.25 縮放,單位 mm)----
    arrow_max=0.020; bmax=max(Bm);
    bxz=hypot(Bx,Bz); bxz(bxz==0)=1e-12;
    scl=arrow_max.*(Bm./bmax).^0.25./bxz;

    % ---- 電荷位置(single = l_hat·dhat；eighteen = l_hat·R_actᵀ·(Pc_base+E))----
    tip=[c.pole_tip_x; c.pole_tip_y; c.pole_tip_z_wp]; dhat=tip./vecnorm(tip);   % WP frame
    if USE_BIAS
        M = load(fullfile(CAL,'data','long2016_hexapole_halfcut','.mat','calib_current_graded_R150_eighteen.mat'),'l_hat','e');
        R_act = [dhat(:,1),dhat(:,3),dhat(:,5)].';       % actuator frame 旋轉
        Pc_base = R_act*dhat;                            % 理想電荷格(actuator, 正規化)
        Pc = make_Pc(M.e, Pc_base);                      % + 18-param 偏移
        qc = M.l_hat * (R_act.' * Pc(:,pidx));           % 轉回 measure/WP frame(m)
    else
        M = load(fullfile(CAL,'data','long2016_hexapole_halfcut','.mat','calib_current_graded_R150_single.mat'),'l_hat');
        qc = M.l_hat*dhat(:,pidx);
    end
    fprintf('  charge WP = (%.3f, %.3f) mm   l_hat=%.1f um\n', qc(1)*1e3, qc(3)*1e3, M.l_hat*1e6);

    % ---- 磁極 y=0 截面輪廓(下極水平半切平頂、上極傾斜全錐;龍飛模型)----
    th=c.pole_angles(pidx)*pi/180;
    if c.pole_is_lower(pidx)
        pa=[cos(th);sin(th);0]; pv=[0;0;-1]; phalf=true;
    else
        inc=c.upper_incline; pa=[cos(inc)*cos(th);cos(inc)*sin(th);sin(inc)];
        pu=cross(pa,[0;0;1]); pu=pu/norm(pu); pv=cross(pa,pu); phalf=false;
    end
    ptip=[c.pole_tip_x(pidx); c.pole_tip_y(pidx); c.pole_tip_z_wp(pidx)]*1e3;  % mm, WP frame
    prf=c.POLE_TIP_R*1e3; pbeta=atan2(c.POLE_R,c.POLE_CONE_LEN); pL=3.5;
    [pox,poz]=pole_outline_xz(ptip,pa,pv,prf,pbeta,pL,phalf);

    % ---- 打包 ----
    S.XL=XL; S.ZL=ZL; S.XT=XT; S.ZT=ZT;
    S.Xs=Xs; S.Zs=Zs; S.Uq=Bx.*scl; S.Wq=Bz.*scl; S.Bm_mT=Bm*1e3;
    S.qcx=qc(1)*1e3; S.qcz=qc(3)*1e3;  S.pox=pox; S.poz=poz;
    S.CLIM=ceil(max(S.Bm_mT)/50)*50;
end

% ============================================================================
function render_panel_into(ax, S, CLIM, FS)
% 把 load_panel 的結構畫進 ax（quiver 依給定 CLIM 分 bin）。
    hold(ax,'on');
    patch(ax, S.pox, S.poz, [0.82 0.84 0.88], 'FaceAlpha',0.30, 'EdgeColor',[0.28 0.30 0.36], 'LineWidth',2.2);
    nb=28; edges=linspace(0,CLIM,nb+1); cmap=turbo(nb); lw=linspace(0.5,2.2,nb);
    for k=1:nb
        if k<nb, m=S.Bm_mT>=edges(k)&S.Bm_mT<edges(k+1); else, m=S.Bm_mT>=edges(k); end
        if any(m), quiver(ax,S.Xs(m),S.Zs(m),S.Uq(m),S.Wq(m),0,'Color',cmap(k,:),'LineWidth',lw(k),'MaxHeadSize',0.35); end
    end
    plot(ax, S.qcx, S.qcz, 'o', 'MarkerFaceColor',[1 0.30 0.65], 'MarkerEdgeColor','k', 'MarkerSize',22, 'LineWidth',1.6);
    axis(ax,'equal'); xlim(ax,S.XL); ylim(ax,S.ZL);
    box(ax,'on'); grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',3.0,'TickLength',[.018 .018]);
    set(ax,'XTick',S.XT,'YTick',S.ZT);                 % 刻度數字維持 Helvetica(不改字體)
    colormap(ax,turbo); clim(ax,[0 CLIM]);
    ax.Toolbar.Visible='off'; hold(ax,'off');
end

% ============================================================================
function style_cbar(cb, FS)
% colorbar 樣式：軸標題標準數學字體(\mathbf CM)、刻度數字 Helvetica 粗體。
    cb.FontSize=FS; cb.FontWeight='bold';
    cb.Label.Interpreter='latex'; cb.Label.String='$\mathbf{|B|\;(mT)}$'; cb.Label.FontSize=FS;
end

% ============================================================================
function Pc = make_Pc(e17, Pc_base)
% 電荷格 Pc = Pc_base + E(e)（含 e6z 約束；與 solve_current/plot_err_hist 一致）。
    if isempty(e17) || all(e17(:)==0)
        Pc = Pc_base;  return;
    end
    E = zeros(3, 6);
    E(:,1) = e17(1:3);     E(:,2) = e17(4:6);
    E(:,3) = e17(7:9);     E(:,4) = e17(10:12);
    E(:,5) = e17(13:15);
    E(1,6) = e17(16);      E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);
    Pc = Pc_base + E;
end

% ============================================================================
function [px, pz] = pole_outline_xz(tip, a, v, rf, beta, L, half)
% 磁極 y=0 截面外框(x,z);a=極軸、v=面內⊥、rf 圓角、beta 半錐角、L 錐長;half=下極半切平頂。
    ts=rf*(1+sin(beta)); rt=rf*cos(beta);
    psi=linspace(0,pi/2+beta,16); nt=24; t=linspace(ts,L,nt);
    axoff=[rf*(1-cos(psi)),  t(2:end)];
    radm =[rf*sin(psi),      rt+(t(2:end)-ts)*tan(beta)];
    edgeP = tip + a.*axoff + v.*radm;        % +v 側錐邊
    if half
        axline = tip + a.*axoff;             % 半切平頂 = 軸線(y=0 截面)
        poly = [axline, fliplr(edgeP)];
    else
        edgeM = tip + a.*axoff - v.*radm;    % -v 側錐邊(全錐)
        poly = [edgeP, fliplr(edgeM)];
    end
    px = poly(1,:); pz = poly(3,:);
end
