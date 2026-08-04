function plot_circuit_side(USE_BIAS, WITH_DIST, CHARGE_SRC, DUAL)
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
    % 用法:plot_circuit_side(USE_BIAS, WITH_DIST, 'maxwell') → 電荷改 Maxwell、場仍 APDL
    if nargin < 2, WITH_DIST = false; end
    if nargin < 3, CHARGE_SRC = 'apdl'; end
    if nargin < 4, DUAL = false; end
    % [ADDED] DUAL=true：一張圖同時畫「有 e(bias,粉)」+「沒 e(fix,深藍)」電荷 + x_a 虛線軸、視野含 WP 中心、無文字標籤。
    % [ADDED] CHARGE_SRC 只換「電荷位置(粉紅點)」的來源校正 .mat:'apdl' | 'maxwell'。
    %   ⚠ **場(箭頭/輪廓/colorbar)一律用 APDL 的 graded .dat,不隨 CHARGE_SRC 改**——
    %   使用者拍板:磁路箭頭要 APDL 資料,只有磁荷位置改用 Maxwell 校正結果。
    clc;
    here   = fileparts(mfilename('fullpath'));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';   % 場來源固定 APDL
    addpath(fullfile(CAL,'function'));  addpath(fullfile(CAL,'common_path'));
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    c = mt_constants();
    FS = 30;
    bstr = ''; if USE_BIAS, bstr = '_bias'; end
    cstr = ''; if strcmpi(CHARGE_SRC,'maxwell'), cstr = '_maxwell'; end   % [ADDED] 電荷來源後綴
    % ---- 佈局（manual pixel；兩 panel 共用 box 高度 H → 邊框對齊；各自 colorbar）----
    H=680; y0=130; leftm=110; cbgap=22; cblab=140; midgap=120; rightm=30;
    CBW_RATIO = 0.009;                                  % colorbar 寬度佔比(cbw/figW)，定案；與 flux 共用同值

    % [MODIFIED 2026-08-03] 視野模式三選一：base / dist(含 WP 中心) / dual(DUAL 專用，見 win())
    MODE = 'base';  if WITH_DIST, MODE = 'dist'; end;  if DUAL, MODE = 'dual'; end
    S1 = load_panel(1, USE_BIAS, c, CAL, MODE, CHARGE_SRC);
    S2 = load_panel(2, USE_BIAS, c, CAL, MODE, CHARGE_SRC);
    if DUAL                                             % [ADDED] 同時算 有e(bias)+沒e(fix) 兩電荷位置(Maxwell)
        qb1=charge_xz(1,true ,c,CAL,CHARGE_SRC); qf1=charge_xz(1,false,c,CAL,CHARGE_SRC);
        qb2=charge_xz(2,true ,c,CAL,CHARGE_SRC); qf2=charge_xz(2,false,c,CAL,CHARGE_SRC);
        S1.qbx=qb1(1)*1e3; S1.qbz=qb1(3)*1e3; S1.qfx=qf1(1)*1e3; S1.qfz=qf1(3)*1e3;
        S2.qbx=qb2(1)*1e3; S2.qbz=qb2(3)*1e3; S2.qfx=qf2(1)*1e3; S2.qfz=qf2(3)*1e3;
        tip=[c.pole_tip_x;c.pole_tip_y;c.pole_tip_z_wp]; dh=tip(:,1)/norm(tip(:,1));  % dhat(P1)=+x_a
        axd=[dh(1);dh(3)]; axd=axd/norm(axd);          % actuator x_a 軸方向(x-z 投影,單位)
        S1.axd=axd; S2.axd=axd;
        % 視野/刻度已由 win(pidx,'dual') 給定（兩 panel 等框寬 + 對稱留白，見該函式）。
    end
    CMAX = max(S1.CLIM, S2.CLIM);                       % P1/P2 共用色階(同 flux 合併圖;弱場 panel 顯冷色)
    SL = S2;  SR = S1;                                  % [MODIFIED] 左右對調：P2 在左、P1 在右（merged/dist 一致）
    w1 = H*diff(SL.XL)/diff(SL.ZL);   w2 = H*diff(SR.XL)/diff(SR.ZL);
    base = leftm + w1 + midgap + w2 + cbgap + cblab + rightm;
    cbw  = CBW_RATIO*base/(1-CBW_RATIO);               % 單一 colorbar，解 cbw/figW = ratio
    figW = base + cbw;   figH = y0 + H + 60;
    x1 = leftm;   x2 = x1 + w1 + midgap;

    fig = figure('Color','w','Units','pixels','Position',[20 40 figW figH]);
    ax1 = axes(fig,'Units','pixels');  render_panel_into(ax1, SL, CMAX, FS, WITH_DIST, DUAL);  ax1.Position=[x1 y0 w1 H];
    ax2 = axes(fig,'Units','pixels');  render_panel_into(ax2, SR, CMAX, FS, WITH_DIST, DUAL);  ax2.Position=[x2 y0 w2 H];
    cb  = colorbar(ax2,'Units','pixels');  cb.Position=[x2+w2+cbgap y0 cbw H];  ax2.Position=[x2 y0 w2 H];  style_cbar(cb, FS);
    if DUAL, cb.Label.String = ''; end   % [MODIFIED 2026-08-03] DUAL 拿掉 colorbar 標題(只留色階+數字)

    % [MODIFIED 2026-08-03] 檔名明確標「電荷來源_模型」（場固定 APDL → 只標電荷來源；使用者拍板）。
    %   dualcharge 同時畫有 e / 沒 e 兩顆電荷 → 不標模型。
    mstr = 'single'; if USE_BIAS, mstr = 'eighteen'; end
    csrc = lower(CHARGE_SRC);
    if DUAL
        outp = fullfile(figdir, sprintf('circuit_side_dualcharge_charge-%s.png', csrc));
    elseif WITH_DIST
        outp = fullfile(figdir, sprintf('circuit_side_dist_charge-%s_%s.png', csrc, mstr));
    else
        outp = fullfile(figdir, sprintf('circuit_side_merged_charge-%s_%s.png', csrc, mstr));
    end
    exportgraphics(fig, outp, 'Resolution', 150);
    fprintf('wrote %s\n', outp);
end

% ============================================================================
function [XL,ZL,XT,ZT] = win(pidx, MODE)
% per-pole 視野 + 刻度（WP frame；z 已平移）。x 起終點都標；z 只留內縮 tick。
% MODE = 'base' | 'dist' | 'dual'
%   'dist'：視野必須含 WP 中心(0,0)。P1 左界由 0.2 → 0（中心落在左上角）；
%           P2 原本 XL=[-1.2 0]、ZL=[0 1]，中心已在右下角，不必改。
%   'dual'：[ADDED 2026-08-03] DUAL 專用。兩 panel **等框寬**（w=H·Δx/Δz，z 各自不動）
%           + **兩端留白對稱**；tick 3/3 奇數、步長 0.4。P2 的 Δx 由 P1 反推，
%           右界因此推到 x≈+0.39（超出 'dist' 取樣範圍 → 自己一份快取）。
    if nargin < 2, MODE = 'base'; end
    switch pidx
        case 1
            switch MODE
                case 'dual'                              % P1 右邊內縮(1.2→0.9)：留白 0.1/0.1
                    % [MODIFIED 2026-08-03] z 刻度改含 0：-0.5/-0.25/0（3 個奇數、步長 0.25、上下留白 0.1/0.1）
                    XL=[-0.1 0.9]; ZL=[-0.6 0.1]; XT=[0 0.4 0.8];        ZT=[-0.5 -0.25 0];
                case 'dist'                              % 原點在左上角 → x 左、z 上各留 0.1mm 讓「+」完整
                    XL=[-0.1 1.2]; ZL=[-0.6 0.1]; XT=[0 0.4 0.8 1.2];    ZT=[-0.5 -0.3 -0.1];
                otherwise
                    XL=[0.2 0.9]; ZL=[-0.6 0];  XT=[0.2 0.4 0.6 0.9];    ZT=[-0.5 -0.3 -0.1];
            end
        case 2
            switch MODE
                case 'dual'
                    ZL=[-0.1 1];  XT=[-0.8 -0.4 0];  ZT=[0 0.5 1];
                    [X1,Z1] = win(1,'dual');             % 框寬對齊 P1：Δx = (Δx₁/Δz₁)·Δz
                    dx = diff(X1)/diff(Z1)*diff(ZL);
                    m  = (dx - (XT(end)-XT(1)))/2;       % 兩端留白對稱(≈0.386，接近 tick 步長 0.4)
                    XL = [XT(1)-m, XT(end)+m];
                case 'dist'                              % 原點在右下角 → x 右、z 下各留 0.1mm
                    XL=[-1.5 0.1]; ZL=[-0.1 1];  XT=[-1.5 -1 -0.5 0];    ZT=[0.25 0.5 0.75];
                otherwise
                    XL=[-1.2 0];  ZL=[0 1];     XT=[-1.2 -0.8 -0.4 0];   ZT=[0.25 0.5 0.75];
            end
    end
end

% ============================================================================
function S = load_panel(pidx, USE_BIAS, c, CAL, MODE, CHARGE_SRC)
% 載入單極場 + grid-sample(不內插) + 電荷位置 + 磁極輪廓；回傳繪圖用結構(含自然 CLIM)。
% [ADDED] 抽稀後的繪圖結構快取:載一份 .dat 要數分鐘,純改樣式不需重載 → 有快取就用。
%   要強制重算把 circuit_side_panel*.mat 刪掉。
    switch pidx, case 1, cname='coil1'; case 2, cname='coil5'; end
    [XL,ZL,XT,ZT] = win(pidx, MODE);

    here   = fileparts(mfilename('fullpath'));   % 快取鍵含視野模式（grid-sample 隨 XL/ZL 而異）
    cachef = fullfile(here, sprintf('circuit_side_panel%d_%s_%s.mat', pidx, ...
                      ternary(USE_BIAS,'bias','fix'), MODE));
                      % 快取只存「場」(APDL graded,與 CHARGE_SRC 無關);電荷每次現算後覆寫,
                      % 這樣換電荷來源不必重讀大 .dat。
    if exist(cachef,'file')
        L = load(cachef,'S');  S = L.S;
        if isequal(S.XL,XL) && isequal(S.ZL,ZL)          % 視野一致才可用(grid-sample 綁 XL/ZL)
            S.XT=XT; S.ZT=ZT;                            % 刻度可直接更新(不影響取樣)
            [qc, qcr_um] = charge_xz(pidx, USE_BIAS, c, CAL, CHARGE_SRC);   % [ADDED] 電荷現算(可切來源)
            S.qcx=qc(1)*1e3; S.qcz=qc(3)*1e3; S.qcr_um=qcr_um;
            fprintf('P%d loaded cache %s  (charge src=%s, r=%.1f um)\n', pidx, cachef, CHARGE_SRC, S.qcr_um);
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
    [qc, qcr_um] = charge_xz(pidx, USE_BIAS, c, CAL, CHARGE_SRC);   % [MODIFIED] 抽成 local(快取命中時也要算)

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
function [qc, qcr_um] = charge_xz(pidx, USE_BIAS, c, CAL, CHARGE_SRC)
% [ADDED] 電荷位置(single = l̂·d̂；eighteen = l̂·R_actᵀ·(Pc_base+E))。
%   CHARGE_SRC='apdl'|'maxwell' 只換讀哪一份 R150 校正 .mat——**場永遠是 APDL,不受此影響**。
%   兩分支幾何(tip40um)完全相同 → d̂ / R_act / Pc_base 共用;e 的欄一律 paper 序 P1..P6。
    if strcmpi(CHARGE_SRC,'maxwell')
        CDIR  = fullfile('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell', ...
                         'data','long2016_hexapole_halfcut','.mat');
        cfile = @(tag) fullfile(CDIR, sprintf('calib_current_maxwell_R150_%s.mat', tag));
    else
        cfile = @(tag) fullfile(CAL,'data','long2016_hexapole_halfcut','.mat', ...
                         sprintf('calib_current_graded_R150_%s.mat', tag));
    end
    tip = [c.pole_tip_x; c.pole_tip_y; c.pole_tip_z_wp];  dhat = tip./vecnorm(tip);   % WP frame
    if USE_BIAS
        M = load(cfile('eighteen'),'l_hat','e');
        R_act   = [dhat(:,1),dhat(:,3),dhat(:,5)].';     % actuator frame 旋轉
        Pc_base = R_act*dhat;                            % 理想電荷格(actuator, 正規化)
        Pc      = make_Pc(M.e, Pc_base);                 % + 18-param 偏移
        qc      = M.l_hat * (R_act.' * Pc(:,pidx));      % 轉回 measure/WP frame(m)
    else
        M  = load(cfile('single'),'l_hat');
        qc = M.l_hat*dhat(:,pidx);
    end
    qcr_um = norm(qc)*1e6;                               % 電荷到 WP 中心的 3D 總距離 [µm]（含 y 分量）
    fprintf('  [charge src=%s] WP = (%.3f, %.3f) mm   l_hat=%.1f um   |r|=%.1f um\n', ...
            CHARGE_SRC, qc(1)*1e3, qc(3)*1e3, M.l_hat*1e6, qcr_um);
end

% ============================================================================
function o = ternary(cond, a, b)
    if cond, o = a; else, o = b; end
end

% ============================================================================
function render_panel_into(ax, S, CLIM, FS, WITH_DIST, DUAL)
% 把 load_panel 的結構畫進 ax（quiver 依給定 CLIM 分 bin）。
    if nargin < 6, DUAL = false; end
    hold(ax,'on');
    patch(ax, S.pox, S.poz, [0.82 0.84 0.88], 'FaceAlpha',0.30, 'EdgeColor',[0.28 0.30 0.36], 'LineWidth',2.2);
    nb=28; edges=linspace(0,CLIM,nb+1); cmap=turbo(nb); lw=linspace(0.5,2.2,nb);
    for k=1:nb
        if k<nb, m=S.Bm_mT>=edges(k)&S.Bm_mT<edges(k+1); else, m=S.Bm_mT>=edges(k); end
        if any(m), quiver(ax,S.Xs(m),S.Zs(m),S.Uq(m),S.Wq(m),0,'Color',cmap(k,:),'LineWidth',lw(k),'MaxHeadSize',0.35); end
    end
    if DUAL
        draw_axis_line(ax, S);                           % x_a 虛線軸(過 WP、+x_a 端箭頭；無文字)
        plot(ax, 0, 0, '+', 'Color','k', 'MarkerSize',24, 'LineWidth',3.0);   % WP 中心十字
        plot(ax, S.qfx, S.qfz, 'o', 'MarkerFaceColor',[0.12 0.24 0.45], 'MarkerEdgeColor','k', 'MarkerSize',20, 'LineWidth',1.6);  % 沒 e (fix) 深藍
        plot(ax, S.qbx, S.qbz, 'o', 'MarkerFaceColor',[1 0.30 0.65],   'MarkerEdgeColor','k', 'MarkerSize',20, 'LineWidth',1.6);  % 有 e (bias) 粉
    else
        plot(ax, S.qcx, S.qcz, 'o', 'MarkerFaceColor',[1 0.30 0.65], 'MarkerEdgeColor','k', 'MarkerSize',22, 'LineWidth',1.6);
    end
    axis(ax,'equal'); xlim(ax,S.XL); ylim(ax,S.ZL);
    if nargin >= 5 && WITH_DIST && ~DUAL
        plot(ax, 0, 0, '+', 'Color','k', 'MarkerSize',26, 'LineWidth',3.0);   % WP 中心(視野內)
        draw_center_dist(ax, S, FS);                     % 虛線指向 WP 中心 + 3D 總距離數值
    end
    box(ax,'on'); grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',3.0,'TickLength',[.018 .018]);
    set(ax,'XTick',S.XT,'YTick',S.ZT);                 % 刻度數字 Helvetica 粗體(維持原字體)
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
function draw_axis_line(ax, S)
% [ADDED] actuator x_a 軸虛線：過 WP 原點、沿 dhat(P1) 的 x-z 投影,裁到框內;+x_a 端畫箭頭;無文字。
    d = S.axd(:);                                        % 單位方向(+x_a)
    ts = [];
    if abs(d(1))>1e-9, ts=[ts S.XL(1)/d(1) S.XL(2)/d(1)]; end
    if abs(d(2))>1e-9, ts=[ts S.ZL(1)/d(2) S.ZL(2)/d(2)]; end
    tin = [];
    for t = ts
        p = t*d;
        if p(1)>=S.XL(1)-1e-6 && p(1)<=S.XL(2)+1e-6 && p(2)>=S.ZL(1)-1e-6 && p(2)<=S.ZL(2)+1e-6
            tin(end+1) = t; %#ok<AGROW>
        end
    end
    if numel(tin) < 2, return; end
    tlo=min(tin); thi=max(tin);
    % [MODIFIED 2026-08-03 依參考圖] 虛線**畫滿整個視野**：從 -x_a 端邊界(tlo)一路拉到 +x_a 端；
    %   箭頭固定放 **+x_a 端(thi，右下角附近)**、朝 +x_a。
    %   原作法(頭放「近 WP 原點」那端 + 0.62 內縮)會讓 P1 的頭跑到左上、線在該端被截短 → 與參考圖相反。
    tn = thi;  tf = tlo;                                 % 頭在 +x_a 端；線起於 -x_a 端邊界
    plo = tf*d;  phi = tn*d;                             % 先放到邊界，下方防裁迴圈再往內退
    col=[0.20 0.30 0.45];
    pe = phi;  hl=0.06*diff(S.XL); a=26*pi/180; b=-d;   % [MODIFIED] 固定 -d：wings 沿 b 張開 ⇒ 視覺箭頭指 -b = +x_a
    Rm=@(th)[cos(th) -sin(th); sin(th) cos(th)];
    Wg = [Rm(a)*b, Rm(-a)*b];                            % 2×2：每欄一支翼的方向
    % [MODIFIED 2026-08-03] 箭頭整顆(頂點 + 兩翼端)都不可被框裁掉：頭起於 +x_a 邊界，沿 -x_a 往框內退到全露。
    mg = 0.012*diff(S.XL);                               % 離框的最小留白
    for it = 1:400
        tp = [pe, pe + Wg*hl];                           % 2×3：頂點 + 兩支翼端
        if all(tp(1,:) >= S.XL(1)+mg) && all(tp(1,:) <= S.XL(2)-mg) && ...
           all(tp(2,:) >= S.ZL(1)+mg) && all(tp(2,:) <= S.ZL(2)-mg), break; end
        pe = pe - 0.01*hl*d;                             % 沿 -x_a 往框內退
    end
    plot(ax,[plo(1) pe(1)],[plo(2) pe(2)],'--','Color',col,'LineWidth',2.5);   % 虛線收在(調整後的)箭頭底
    for k = 1:2
        w = Wg(:,k);
        plot(ax,[pe(1) pe(1)+w(1)*hl],[pe(2) pe(2)+w(2)*hl],'-','Color',col,'LineWidth',2.5);
    end
end

% ============================================================================
function style_cbar(cb, FS)
% colorbar 樣式：軸標題標準數學字體(\mathbf CM)、刻度數字 Helvetica 粗體。
    cb.FontSize=FS; cb.FontWeight='bold';   % colorbar 數字 Helvetica 粗體(維持原字體)
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
