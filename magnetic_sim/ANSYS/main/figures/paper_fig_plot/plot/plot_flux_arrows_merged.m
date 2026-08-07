function plot_flux_arrows_merged()
% plot_flux_arrows_merged -- long2016 半切六極 coil5=P2 俯視磁通:巨觀 + WP zoom 併一張
% =========================================================================
%   左 panel = 巨觀俯視(裁到 ±60mm、'all' dataset)、右 panel = WP zoom(±1mm、'wp')。
%   同一激發(coil5=P2、上極 source s=+1)、**共用 |B| 色階 + 單一 colorbar**(shared clim)。
%   畫法沿用龍飛 fig(a):全節點、格取 |B|max 真實節點(不內插)、|B|^0.25 長度、turbo 依 |B|(mT)。
%   軸:**單位 mm**、LaTeX 數學字型(對齊 3D 幾何圖 x_a/y_a/z_a 風格);tick 奇數含原點、
%   去端點(角落不標數字);字體 36。左 panel 留 6 黑圈+軛環+P1..P6;右 panel **不標 P**。
%   輸出 → figures/paper_fig/Section2_A/flux_arrows_merged.png(覆蓋迭代)。
% =========================================================================
    clc;
    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CAL,'function'));  addpath(fullfile(CAL,'common_path'));
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    cnst = mt_constants();

    cname='coil5';  pidx=2;                 % coil5 = P2(上極 source）
    s = 1 - 2*cnst.pole_is_lower(pidx);     % = +1(不翻)
    ddir = ansys_path('long2016_hexapole_halfcut','data','graded',cname);

    % ---- 取樣兩 panel（各自比例的降取樣，回傳仍在 µm；後面統一轉 mm） ----
    dA = import_ansys_data(ddir,'all',cname);
    fprintf('macro: matched=%d, |B|max=%.4f T\n', numel(dA.x), max(dA.bsum));
    A  = sample_macro(dA, s);
    dB = import_ansys_data(ddir,'wp',cname);
    fprintf('zoom : matched=%d, |B|max=%.4f T\n', numel(dB.x), max(dB.bsum));
    B  = sample_zoom(dB, s);

    % ---- µm → mm（位置+箭頭同乘，比例不變） ----
    for f = {'xq','yq','uq','vq'}, A.(f{1})=A.(f{1})*1e-3;  B.(f{1})=B.(f{1})*1e-3;  end

    % ---- 共用色階（兩 panel 取最大） ----
    CLIM = ceil(max([max(A.bmag_mT), max(B.bmag_mT)])/50)*50;
    fprintf('shared CLIM = %d mT\n', CLIM);

    % ---- 手動 pixel 佈局：兩 panel 共用 box 高度 H → 上下邊框對齊（右對齊左）----
    FS = 36;
    BHX = 68;  BHY = 60;                                % [mm] 左半框:x 稍寬、y 較窄(去外圈空白)
    H = 680;  y0 = 120;  leftm = 130;  midgap = 150;  cbgap = 22;  cblab = 150;
    CBW_RATIO = 0.009;                                 % colorbar 寬度佔比(cbw/figW)；與 circuit 圖共用同值
    w1 = H*(2*BHX)/(2*BHY);                             % 左 aspect = 136/120 → box 寬(共用 H)
    w2 = H*(2)/(2);                                     % 右 ±1 → 正方(box 寬 = H)
    x1 = leftm;  x2 = x1 + w1 + midgap;
    base = leftm + w1 + midgap + w2 + cbgap + cblab + 30;
    cbw  = CBW_RATIO*base/(1-CBW_RATIO);   figW = base + cbw;  figH = y0 + H + 60;
    fig = figure('Color','w','Units','pixels','Position',[20 30 figW figH]);

    % ---- 左:巨觀 裁到 ±60mm ----
    ax1 = axes(fig,'Units','pixels');  hold(ax1,'on');
    render_panel(ax1, A, CLIM);
    draw_structure(ax1, cnst);                          % mm 幾何
    axis(ax1,'equal');  xlim(ax1,[-BHX BHX]);  ylim(ax1,[-BHY BHY]);   % axis equal 保圓
    box(ax1,'on');  grid(ax1,'off');
    set(ax1,'FontSize',FS,'FontWeight','bold','LineWidth',3.0,'TickLength',[.012 .012]);
    set(ax1,'XTick',[-40 0 40],'YTick',[-40 0 40]);    % 3 個(奇數)含 0、去端點(無軸標題)
    ax1.Toolbar.Visible='off';  hold(ax1,'off');
    ax1.Position = [x1 y0 w1 H];                        % 共用 H → 上下邊框對齊

    % ---- 右:WP zoom ±1mm（不標 P） ----
    ax2 = axes(fig,'Units','pixels');  hold(ax2,'on');
    render_panel(ax2, B, CLIM);
    plot(ax2, 0,0,'k+','MarkerSize',16,'LineWidth',2.6);   % WP 十字
    axis(ax2,'equal');  xlim(ax2,[-1 1]);  ylim(ax2,[-1 1]);
    box(ax2,'on');  grid(ax2,'off');
    set(ax2,'FontSize',FS,'FontWeight','bold','LineWidth',3.0,'TickLength',[.012 .012]);
    set(ax2,'XTick',[-0.5 0 0.5],'YTick',[-0.5 0 0.5]); % 3 個(奇數)含 0、去端點(無軸標題)
    ax2.Toolbar.Visible='off';  hold(ax2,'off');
    ax2.Position = [x2 y0 w2 H];                        % 同 H、同 y0 → 與左 panel 邊框齊平

    % ---- 單一共用 colorbar（pixel 定位，貼右 panel、同高 H）----
    colormap(fig,turbo);  clim(ax1,[0 CLIM]);  clim(ax2,[0 CLIM]);
    cb = colorbar(ax2,'Units','pixels');
    cb.Position = [x2+w2+cbgap y0 cbw H];  ax2.Position = [x2 y0 w2 H];
    cb.FontSize=30;  cb.FontWeight='bold';                                  % tick 數字:一般字型、粗體
    cb.Label.Interpreter='latex';  cb.Label.String='$\mathbf{|B|\;(mT)}$';  % 只有標籤用 LaTeX 數學字型 + 粗體

    out = fullfile(figdir,'flux_arrows_merged.png');
    exportgraphics(fig, out, 'Resolution', 110);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function A = sample_macro(d, s)
% 巨觀 ±80000µm:每 2000µm 格取 |B|max 真實節點;|B|^0.25 長度(arrow_max=4000µm @1T)。
    bx=s*d.bx; by=s*d.by; x_um=d.x*1e6; y_um=d.y*1e6; bsum=d.bsum;
    RNG=80000; gc=80; xe=linspace(-RNG,RNG,gc+1); ye=linspace(-RNG,RNG,gc+1);
    ix=discretize(x_um,xe); iy=discretize(y_um,ye);
    vb=~isnan(ix)&~isnan(iy);
    cid=zeros(numel(x_um),1); cid(vb)=(ix(vb)-1)*gc+iy(vb);
    uc=unique(cid(cid>0)); idx=zeros(numel(uc),1);
    for k=1:numel(uc), cn=find(cid==uc(k)); [~,b]=max(bsum(cn)); idx(k)=cn(b); end
    keep=bsum(idx)>1e-4; idx=idx(keep);
    bmag=bsum(idx);
    arrow_max=4000; bxy=hypot(bx(idx),by(idx)); bxy(bxy==0)=1e-10;
    scl=arrow_max.*bmag.^0.25./bxy;
    A.xq=x_um(idx); A.yq=y_um(idx); A.uq=bx(idx).*scl; A.vq=by(idx).*scl;
    A.bmag_mT=bmag*1e3; A.lw=linspace(0.40,1.75,28);
end

% ============================================================================
function B = sample_zoom(d, s)
% WP zoom ±1000µm、|z_wp|<1mm:每 ~23µm 格取 |B|max 真實節點;|B|^0.25 長度(arrow_max=40µm)。
    bx=s*d.bx; by=s*d.by; bmag=d.bsum; x=d.x*1e6; y=d.y*1e6; zwp=d.z*1e3+12.711;
    XR=[-1000 1000]; YR=[-1000 1000];
    inbox = x>=XR(1)&x<=XR(2)&y>=YR(1)&y<=YR(2)&abs(zwp)<1.0;
    ncx=88; ncy=88; xe=linspace(XR(1),XR(2),ncx+1); ye=linspace(YR(1),YR(2),ncy+1);
    ix=discretize(x,xe); iy=discretize(y,ye);
    good=inbox&~isnan(ix)&~isnan(iy);
    cid=(ix-1)*ncy+iy; vi=find(good);
    [~,~,g]=unique(cid(vi)); sel=zeros(max(g),1);
    for k=1:max(g), cn=vi(g==k); [~,b]=max(bmag(cn)); sel(k)=cn(b); end
    keep=bmag(sel)>1e-4; sel=sel(keep);
    Bm=bmag(sel); bmax=max(Bm);
    arrow_max=40; bxy=hypot(bx(sel),by(sel)); bxy(bxy==0)=1e-12;
    scl=arrow_max.*(Bm./bmax).^0.25./bxy;
    B.xq=x(sel); B.yq=y(sel); B.uq=bx(sel).*scl; B.vq=by(sel).*scl;
    B.bmag_mT=Bm*1e3; B.lw=linspace(0.5,1.9,28);
end

% ============================================================================
function render_panel(ax, S, CLIM)
% turbo 依 |B|(mT)分 28 bin 疊 quiver、linewidth ramp。
    nb=28; edges=linspace(0,CLIM,nb+1); cmap=turbo(nb);
    for k=1:nb
        if k<nb, m=S.bmag_mT>=edges(k)&S.bmag_mT<edges(k+1); else, m=S.bmag_mT>=edges(k); end
        if any(m), quiver(ax,S.xq(m),S.yq(m),S.uq(m),S.vq(m),0, ...
                          'Color',cmap(k,:),'LineWidth',S.lw(k),'MaxHeadSize',0.3); end
    end
end

% ============================================================================
function draw_structure(ax, cnst)
% 巨觀(mm):6 磁極黑圓圈(@ YOKE_MID_R、半徑 PROT_R)+ P1..P6 + 紅虛線軛環 + WP 十字。
    Rin=cnst.YOKE_IN_R*1e3; Rout=cnst.YOKE_OUT_R*1e3; Rmid=cnst.YOKE_MID_R*1e3; PR=cnst.PROT_R*1e3;
    th=linspace(0,2*pi,240);
    plot(ax, Rin*cos(th),  Rin*sin(th),  '--','Color',[.85 .1 .1],'LineWidth',1.6);
    plot(ax, Rout*cos(th), Rout*sin(th), '--','Color',[.85 .1 .1],'LineWidth',1.6);
    for i=1:6
        a=cnst.pole_angles(i)*pi/180; cx=Rmid*cos(a); cy=Rmid*sin(a);
        plot(ax, cx+PR*cos(th), cy+PR*sin(th), 'k-','LineWidth',2.0);   % 6 磁極黑圈(不標 P)
    end
    plot(ax, 0,0,'k+','MarkerSize',16,'LineWidth',2.6);
end
