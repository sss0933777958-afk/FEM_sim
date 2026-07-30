function plot_circuit_side(USE_BIAS, WITH_DIST)
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
        plot_circuit_side(false,false);  plot_circuit_side(true,false);      % 原圖(不動)
        plot_circuit_side(false,true);   plot_circuit_side(true,true);       % 距離標註版(另存)
        return;
    end
    if nargin < 2, WITH_DIST = false; end
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

    S1 = load_panel(1, USE_BIAS, c, CAL, WITH_DIST);
    S2 = load_panel(2, USE_BIAS, c, CAL, WITH_DIST);
    CMAX = max(S1.CLIM, S2.CLIM);                       % P1/P2 共用色階(同 flux 合併圖;弱場 panel 顯冷色)
    if WITH_DIST, SL = S2;  SR = S1;                    % 距離版:P2 在左、P1 在右(左右對調)
    else,         SL = S1;  SR = S2;  end               % 原圖:P1 在左、P2 在右
    w1 = H*diff(SL.XL)/diff(SL.ZL);   w2 = H*diff(SR.XL)/diff(SR.ZL);
    base = leftm + w1 + midgap + w2 + cbgap + cblab + rightm;
    cbw  = CBW_RATIO*base/(1-CBW_RATIO);               % 單一 colorbar，解 cbw/figW = ratio
    figW = base + cbw;   figH = y0 + H + 60;
    x1 = leftm;   x2 = x1 + w1 + midgap;

    fig = figure('Color','w','Units','pixels','Position',[20 40 figW figH]);
    ax1 = axes(fig,'Units','pixels');  render_panel_into(ax1, SL, CMAX, FS, WITH_DIST);  ax1.Position=[x1 y0 w1 H];
    ax2 = axes(fig,'Units','pixels');  render_panel_into(ax2, SR, CMAX, FS, WITH_DIST);  ax2.Position=[x2 y0 w2 H];
    cb  = colorbar(ax2,'Units','pixels');  cb.Position=[x2+w2+cbgap y0 cbw H];  ax2.Position=[x2 y0 w2 H];  style_cbar(cb, FS);

    if WITH_DIST
        outp = fullfile(figdir, sprintf('circuit_side_dist%s.png', bstr));      % 距離標註版(另存,不覆蓋原圖)
    else
        outp = fullfile(figdir, sprintf('circuit_side_merged%s.png', bstr));    % 原圖
    end
    exportgraphics(fig, outp, 'Resolution', 150);
    fprintf('wrote %s\n', outp);
end

% ============================================================================
function [XL,ZL,XT,ZT] = win(pidx, WITH_DIST)
% per-pole 視野 + 刻度（WP frame；z 已平移）。x 起終點都標；z 只留內縮 tick。
% WITH_DIST=true：視野必須含 WP 中心(0,0)。P1 左界由 0.2 → 0（中心落在左上角）；
%   P2 原本 XL=[-1.2 0]、ZL=[0 1]，中心已在右下角，不必改。
    if nargin < 2, WITH_DIST = false; end
    switch pidx
        case 1
            if WITH_DIST                                 % 原點在左上角 → x 左、z 上各留 0.1mm 讓「+」完整
                XL=[-0.1 1.2]; ZL=[-0.6 0.1]; XT=[0 0.4 0.8 1.2];    ZT=[-0.5 -0.3 -0.1];
            else
                XL=[0.2 0.9]; ZL=[-0.6 0];  XT=[0.2 0.4 0.6 0.9];    ZT=[-0.5 -0.3 -0.1];
            end
        case 2
            if WITH_DIST                                 % 原點在右下角 → x 右、z 下各留 0.1mm
                XL=[-1.5 0.1]; ZL=[-0.1 1];  XT=[-1.5 -1 -0.5 0];    ZT=[0.25 0.5 0.75];
            else
                XL=[-1.2 0];  ZL=[0 1];     XT=[-1.2 -0.8 -0.4 0];   ZT=[0.25 0.5 0.75];
            end
    end
end

% ============================================================================
function S = load_panel(pidx, USE_BIAS, c, CAL, WITH_DIST)
% 載入單極場 + grid-sample(不內插) + 電荷位置 + 磁極輪廓；回傳繪圖用結構(含自然 CLIM)。
% [ADDED] 抽稀後的繪圖結構快取:載一份 .dat 要數分鐘,純改樣式不需重載 → 有快取就用。
%   要強制重算把 circuit_side_panel*.mat 刪掉。
    switch pidx, case 1, cname='coil1'; case 2, cname='coil5'; end
    [XL,ZL,XT,ZT] = win(pidx, WITH_DIST);

    here   = fileparts(mfilename('fullpath'));   % 快取鍵含視野模式（grid-sample 隨 XL/ZL 而異）
    cachef = fullfile(here, sprintf('circuit_side_panel%d_%s_%s.mat', pidx, ...
                      ternary(USE_BIAS,'bias','fix'), ternary(WITH_DIST,'dist','base')));
    if exist(cachef,'file')
        L = load(cachef,'S');  S = L.S;
        if isequal(S.XL,XL) && isequal(S.ZL,ZL)          % 視野一致才可用(grid-sample 綁 XL/ZL)
            S.XT=XT; S.ZT=ZT;                            % 刻度可直接更新(不影響取樣)
            fprintf('P%d loaded cache %s  (charge r=%.1f um)\n', pidx, cachef, S.qcr_um);
            return;
        end
        fprintf('P%d cache 視野不符(舊 x[%g %g] vs 新 x[%g %g]) → 重算\n', pidx, S.XL, XL);
    end

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
    qcr_um = norm(qc)*1e6;                              % [ADDED] 電荷到 WP 中心的「3D 總距離」[µm]（含 y 分量）
    fprintf('  charge WP = (%.3f, %.3f) mm   l_hat=%.1f um   |r|=%.1f um\n', ...
            qc(1)*1e3, qc(3)*1e3, M.l_hat*1e6, qcr_um);

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
    S.qcx=qc(1)*1e3; S.qcz=qc(3)*1e3;  S.qcr_um=qcr_um;  S.pox=pox; S.poz=poz;
    S.CLIM=ceil(max(S.Bm_mT)/50)*50;
    save(cachef, 'S');   fprintf('P%d saved cache %s\n', pidx, cachef);
end

% ============================================================================
function o = ternary(cond, a, b)
    if cond, o = a; else, o = b; end
end

% ============================================================================
function render_panel_into(ax, S, CLIM, FS, WITH_DIST)
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
    if nargin >= 5 && WITH_DIST
        plot(ax, 0, 0, '+', 'Color','k', 'MarkerSize',26, 'LineWidth',3.0);   % WP 中心(視野內)
        draw_center_dist(ax, S, FS);                     % 虛線指向 WP 中心 + 3D 總距離數值
    end
    box(ax,'on'); grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',3.0,'TickLength',[.018 .018]);
    set(ax,'XTick',S.XT,'YTick',S.ZT);                 % 刻度數字維持 Helvetica(不改字體)
    colormap(ax,turbo); clim(ax,[0 CLIM]);
    ax.Toolbar.Visible='off'; hold(ax,'off');
end

% ============================================================================
function draw_center_dist(ax, S, FS)
% [ADDED] 由電荷點朝 WP 中心(原點)畫虛線 + 箭頭 + 標「3D 總距離」[µm]。
%   中心多半在視野外 → 虛線裁到框邊(或到原點為止,取較近者)。
%   數值擺線段中點、沿法線偏移;偏移方向自動選「遠離磁極輪廓」那側,避免壓到極面。
    p0 = [S.qcx; S.qcz];
    d  = -p0;   L0 = norm(d);   if L0 < 1e-9, return; end
    d  = d / L0;

    tc = [];                                             % 裁到框邊的參數 t
    if d(1) > 0, tc(end+1) = (S.XL(2)-p0(1))/d(1); elseif d(1) < 0, tc(end+1) = (S.XL(1)-p0(1))/d(1); end
    if d(2) > 0, tc(end+1) = (S.ZL(2)-p0(2))/d(2); elseif d(2) < 0, tc(end+1) = (S.ZL(1)-p0(2))/d(2); end
    t  = min([tc, L0]);                                  % 不超過原點本身
    pe = p0 + t*d;

    plot(ax, [p0(1) pe(1)], [p0(2) pe(2)], 'k--', 'LineWidth', 2.5);
    hl = 0.055*diff(S.XL);  a = 26*pi/180;  b = -d;      % 箭頭(指向中心)兩翼
    R  = @(th) [cos(th) -sin(th); sin(th) cos(th)];
    for w = [R(a)*b, R(-a)*b]
        plot(ax, [pe(1) pe(1)+w(1)*hl], [pe(2) pe(2)+w(2)*hl], 'k-', 'LineWidth', 2.5);
    end

    pm  = p0 + 0.5*t*d;   n = [-d(2); d(1)];   off = 0.085*diff(S.XL);
    cen = [mean(S.pox); mean(S.poz)];                    % 磁極輪廓形心
    if norm(pm + n*off - cen) < norm(pm - n*off - cen), n = -n; end   % 選遠離磁極那側
    q = pm + n*off;
    q(1) = min(max(q(1), S.XL(1)+0.10*diff(S.XL)), S.XL(2)-0.10*diff(S.XL));   % 夾在框內
    q(2) = min(max(q(2), S.ZL(1)+0.07*diff(S.ZL)), S.ZL(2)-0.07*diff(S.ZL));
    text(ax, q(1), q(2), sprintf('$\\mathbf{%.0f\\;\\mu m}$', S.qcr_um), ...
         'Interpreter','latex', 'FontSize',FS, 'Color','k', ...
         'HorizontalAlignment','center', 'VerticalAlignment','middle');
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
