function plot_p1p2_poles_3d(SOURCE, FLD, XLin, ZLin, YHW, EXC)
% [ADDED 2026-08-14] EXC = 激發哪根極（paper 編號，預設 1）。
%   'apdl' 分支據此挑 graded coil（long2016 map [1,3,6,5,2,4]：P1->coil1、P2->coil5），
%   輸出檔名 p<EXC>_apdl_side.png。'maxwell' 分支則由 .fld 檔名的 p<N> 決定。
% plot_p1p2_poles_3d -- long2016 P1+P2 磁路「側視圖」(x-z, y=0) + sensor + n+ + WP 十字
% =========================================================================
%   側視(x-z, y=0)：P1(下極半切、水平+x)與 P2(上極傾斜全錐)共尖端於 WP(原點)。
%   每 (x,z) 格取最近 y=0 真實節點(不內插)、turbo 依 |B|(mT)。磁極 y=0 截面輪廓灰填+深邊。
%   colorbar 樣式同 circuit_side（style_cbar：粗體 + $\mathbf{|B|\;(mT)}$）。字體 36、無軸標題。
%   場源 SOURCE：
%     'apdl'（預設）— coil1(=P1) graded .dat；限 R=7mm 球內 → p1p2_poles_3d.png
%     'maxwell'     — P1 激發的 Maxwell .fld（給 FLD 路徑）；矩形視野 XLin/ZLin、y 帶 |y|<=YHW
%                     → p1_maxwell_side.png（用 import_maxwell_fld 讀 .fld）
% =========================================================================
    clc;
    if nargin<1 || isempty(SOURCE), SOURCE='apdl'; end
    if nargin<2, FLD=''; end
    if nargin<6 || isempty(EXC), EXC=1; end
    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    % [MODIFIED 2026-08-08] 脫離 backup（規則 no-backup-data）→ live 樹。
    CALROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    % [MODIFIED 2026-08-14] ⚠ Maxwell\function 必須**先**加：addpath 預設 prepend，後加的會蓋前面。
    %   原本順序讓 Maxwell 的 model_config 蓋掉 APDL 版 → 拿到 identity 的 coil→pole map，
    %   P2 被解成 coil2（實際應為 coil5），載到錯的激發資料。
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell\function');   % import_maxwell_fld
    addpath(fullfile(CALROOT,'function'), fullfile(CALROOT,'utils'), fullfile(CALROOT,'common_path'));
    c = model_config('long2016_hexapole_halfcut','tip40um');

    rf = c.POLE_TIP_R*1e3;  beta = atan2(c.POLE_R, c.POLE_CONE_LEN);
    L  = 8;
    switch lower(SOURCE)
        case 'apdl'
            XL=[-6 6]; ZL=[-3 5]; YHW=7; R_SPH=7;                   % [MODIFIED] 視野同 maxwell 版，便於對照
            XT=[-6 -3 0 3 6]; ZT=[-2 0 2 4];
            outname = sprintf('p%d_apdl_side.png', EXC);
        case 'maxwell'
            if isempty(FLD), error('SOURCE=maxwell 需給 .fld 路徑 FLD'); end
            if nargin<3||isempty(XLin), XLin=[-6 6]; end
            if nargin<4||isempty(ZLin), ZLin=[-3 5]; end
            if nargin<5||isempty(YHW),  YHW=4;       end
            XL=XLin; ZL=ZLin; R_SPH=Inf;                            % 矩形視野、不限球
            XT=[-6 -3 0 3 6]; ZT=[-2 0 2 4];
            % [MODIFIED 2026-08-14] 輸出檔名由 .fld 檔名裡的 p<N> 決定（激發哪根極），
            %   不再寫死 p1 —— 例：B_voltage_p2.fld -> p2_maxwell_side.png
            [~, fbn] = fileparts(FLD);
            tk = regexp(fbn, 'p(\d+)', 'tokens', 'once');
            iexc = 1;  if ~isempty(tk), iexc = str2double(tk{1}); end
            outname = sprintf('p%d_maxwell_side.png', iexc);
        otherwise, error('SOURCE 必為 ''apdl'' | ''maxwell''');
    end
    FS = 36;  GREEN=[0.10 0.60 0.20];  RED=[0.85 0.10 0.10];

    % ---- P1 y=0 切面場（每 (x,z) 格取最近 y=0 真實節點 + 縮放好的箭頭分量）----
    [Xs,Zs,Uq,Wq,Bm_mT,CLIM,SL] = field_slice(c, XL, ZL, R_SPH, SOURCE, FLD, YHW, EXC);

    % ---- sensor / n+ / 錐體幾何：一律取自 pole_sensor_geometry（唯一來源）----
    % [MODIFIED 2026-08-14] 原本 sensor 位置是寫死的舊值（P1 4.8047/-1.6189、P2 -3.1454/3.9774），
    %   與正規函式差 15/25 um。改為呼叫函式；同時輪廓改用 CAD 實測的 per-pole 半錐角與倒圓半徑
    %   （原本兩根極共用 mt_constants 的名目值，對兩層都不準）。
    [spos, snor, geo] = pole_sensor_geometry(c);
    to2  = @(v) v([1 3])*1e3;
    spP1 = to2(spos(:,1));   nP1 = snor([1 3],1);
    spP2 = to2(spos(:,2));   nP2 = snor([1 3],2);

    % ---- P1 / P2 的 y=0 截面輪廓（per-pole CAD 幾何）----
    [p1x,p1z] = pole_outline(c, 1, geo.r_tip(1)*1e3, geo.beta(1), L);
    [p2x,p2z] = pole_outline(c, 2, geo.r_tip(2)*1e3, geo.beta(2), L);

    fig = figure('Color','w','Position',[60 60 1280 760]);
    ax  = axes(fig);  hold(ax,'on');

    % 1) 磁極輪廓（灰填 + 深邊，畫在最底）
    patch(ax, p1x, p1z, [0.82 0.84 0.88], 'FaceAlpha',0.30, 'EdgeColor',[0.28 0.30 0.36], 'LineWidth',2.2);
    patch(ax, p2x, p2z, [0.82 0.84 0.88], 'FaceAlpha',0.30, 'EdgeColor',[0.28 0.30 0.36], 'LineWidth',2.2);

    % 2) 磁路箭頭（turbo 依 |B|mT 分 bin；線寬隨場漸粗）
    nb=28; edges=linspace(0,CLIM,nb+1); cmap=turbo(nb); lw=linspace(0.5,2.2,nb);
    drawq = @(X,Z,U,W,M) arrayfun(@(k) qbin(ax,X,Z,U,W,M,edges,cmap,lw,nb,k), 1:nb);
    drawq(Xs,Zs,Uq,Wq,Bm_mT);

    % 2b) [ADDED 2026-08-14] 磁極表面帶（**內插**）：讓磁通看得出從激發極表面轉出、
    %     從接收極表面轉入。離面 0.10/0.24/0.40 mm 三層、沿輪廓每 0.18mm 一點。
    for pp = 1:2
        if pp==1, bx_=p1x; bz_=p1z; else, bx_=p2x; bz_=p2z; end
        [Qx,Qz,Ub,Wb,Bb] = surface_band(bx_, bz_, SL, [0.10 0.24 0.40], 0.18);
        fprintf('  P%d 表面帶箭頭 %d 支（內插）\n', pp, numel(Qx));
        drawq(Qx,Qz,Ub,Wb,Bb);
    end

    % 3) sensor 感測區域（實際尺寸圓柱側視矩形 2R×H，綠）+ n+ 紅箭頭
    draw_sensor2d(ax, spP1, nP1, 0.15, 0.10, GREEN);
    draw_sensor2d(ax, spP2, nP2, 0.15, 0.10, GREEN);
    draw_narrow2d(ax, spP1, nP1, 1.0, RED, 3.5);
    draw_narrow2d(ax, spP2, nP2, 1.0, RED, 3.5);

    % 4) WP（原點）十字：screen-space marker（MarkerSize22 / LineWidth3）
    plot(ax, 0, 0, '+', 'MarkerSize',22, 'Color','k', 'LineWidth',3);

    axis(ax,'equal');  xlim(ax,XL);  ylim(ax,ZL);
    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',3.0,'TickLength',[.018 .018]);
    set(ax,'XTick',XT,'YTick',ZT);                       % 奇數等距、無軸標題（不設 xlabel/ylabel）
    colormap(ax,turbo);  clim(ax,[0 CLIM]);
    cb = colorbar(ax);  style_cbar(cb, FS);
    % --- 手動版面：留右側空間給 colorbar 旋轉標題，避免被 figure 右緣裁掉 ---
    set(ax,'Units','normalized','Position',[0.075 0.11 0.70 0.86]);
    set(cb,'Units','normalized','Position',[0.795 0.13 0.028 0.82]);
    ax.Toolbar.Visible='off';  hold(ax,'off');

    out = fullfile(figdir, outname);
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('wrote %s  (CLIM=%d mT)\n', out, CLIM);
end

% ============================================================================
function [Xs,Zs,Uq,Wq,Bm_mT,CLIM,SL] = field_slice(c, XL, ZL, R_SPH, SOURCE, FLD, YHW, EXC)
% 激發極 y=0 切面：每 (x,z) 格取最近 y=0 的真實節點(不內插)。回箭頭分量(mm)。
%   SOURCE='apdl' → graded coil<j> .dat（j 由 EXC 經 apdl_to_paper_idx 反查）；
%   'maxwell' → import_maxwell_fld(FLD)。
%   SL = 薄板內的原始節點（供磁極表面帶內插用）。
    switch lower(SOURCE)
        case 'apdl'
            j = find(c.apdl_to_paper_idx == EXC, 1);
            assert(~isempty(j), 'apdl_to_paper_idx 找不到 P%d', EXC);
            cn = sprintf('coil%d', j);
            d = import_ansys_data(ansys_path('long2016_hexapole_halfcut','data','graded',cn), 'all', cn);
            fprintf('P%d = %s (graded/all): matched=%d, |B|max=%.4f T\n', EXC, cn, numel(d.x), max(d.bsum));
            s = 1 - 2*c.pole_is_lower(EXC);                  % 下極 -1 / 上極 +1（全 source 慣例）
        case 'maxwell'
            d = import_maxwell_fld(FLD);
            fprintf('P%d = Maxwell .fld: %d pts, |B|max=%.4f T\n', EXC, numel(d.x), max(d.bsum));
            s = 1;                                           % Maxwell 場方向（待驗；若與 APDL 反向改 -1）
    end
    zoff = -c.SPH_OFST*1e3;
    x=d.x*1e3; y=d.y*1e3; z=d.z*1e3+zoff;                 % mm, WP frame
    rWP = sqrt(x.^2+y.^2+z.^2);

    cell = 0.20;                                          % (x,z) 格邊 (mm)
    inwin = x>=XL(1)&x<=XL(2) & z>=ZL(1)&z<=ZL(2) & abs(y)<=YHW & rWP<R_SPH;
    gi = find(inwin);
    xe = XL(1):cell:XL(2);  ze = ZL(1):cell:ZL(2);
    ix = discretize(x(gi),xe);  iz = discretize(z(gi),ze);
    ok = ~isnan(ix)&~isnan(iz);  gi=gi(ok); ix=ix(ok); iz=iz(ok);
    nz = numel(ze)-1;  cid=(ix-1)*nz+iz;  ay=abs(y(gi));
    [uc,~,g] = unique(cid);  sel=zeros(numel(uc),1);
    for k=1:numel(uc), idx=find(g==k); [~,b]=min(ay(idx)); sel(k)=gi(idx(b)); end
    keep = d.bsum(sel)>1e-4;  sel=sel(keep);

    Xs=x(sel);  Zs=z(sel);
    if strcmpi(SOURCE,'maxwell')                         % 打散規則格「格子感」：位置微抖動(場向量仍為真實格點值、不內插)
        rng(7);  J=0.45*cell;
        Xs = Xs + (2*rand(size(Xs))-1)*J;
        Zs = Zs + (2*rand(size(Zs))-1)*J;
    end
    Bx=s*d.bx(sel);  Bz=s*d.bz(sel);  Bm=d.bsum(sel);
    arrow_max=0.30;  bmax=max(Bm);                        % 箭頭長度 ∝ |B|^0.25 (mm)
    bxz=hypot(Bx,Bz);  bxz(bxz==0)=1e-12;
    scl=arrow_max.*(Bm./bmax).^0.25./bxz;
    Uq=Bx.*scl;  Wq=Bz.*scl;  Bm_mT=Bm*1e3;
    CLIM=ceil(max(Bm_mT)/50)*50;

    % [ADDED 2026-08-14] 薄板原始節點（|y| <= YSL），供磁極表面帶做 2D 內插用。
    YSL = 0.30;
    q = x>=XL(1)-0.5 & x<=XL(2)+0.5 & z>=ZL(1)-0.5 & z<=ZL(2)+0.5 & abs(y)<=YSL & rWP<R_SPH;
    SL = struct('x',x(q),'z',z(q),'bx',s*d.bx(q),'bz',s*d.bz(q),'bm',d.bsum(q), ...
                'bmax',bmax,'arrow_max',arrow_max);
    fprintf('  表面帶內插源：|y|<=%.2fmm 薄板節點 %d 個\n', YSL, nnz(q));
end

% ============================================================================
function qbin(ax, X, Z, U, W, M, edges, cmap, lw, nb, k)
% 依 |B| bin 畫一批箭頭（turbo 分色、線寬隨場漸粗）
    if k < nb, m = M>=edges(k) & M<edges(k+1); else, m = M>=edges(k); end
    if any(m)
        quiver(ax, X(m), Z(m), U(m), W(m), 0, 'Color',cmap(k,:), ...
               'LineWidth',lw(k), 'MaxHeadSize',0.35);
    end
end

% ============================================================================
function [Qx,Qz,Uq,Wq,Bm_mT] = surface_band(px, pz, SL, offs, ds)
% [ADDED 2026-08-14] 沿磁極 y=0 截面輪廓「外側」鋪一帶查詢點，內插磁場後回箭頭分量。
%   目的（使用者要求）：讓磁通看起來是從激發極表面轉出來、從接收極表面轉進去。
%   ⚠ 這一帶是 **scatteredInterpolant 線性內插**（不是節點原值）——依 plot-real-nodes
%     規則，使用者明確要求時才可用，且必須在回覆/說明中標示為內插。
%   px,pz 磁極截面多邊形；offs 離面距離陣列 [mm]；ds 沿邊界取樣間距 [mm]。
    P = [px(:) pz(:)];
    if ~isequal(P(1,:), P(end,:)), P = [P; P(1,:)]; end
    seg = diff(P);  Lg = hypot(seg(:,1), seg(:,2));
    Q = zeros(0,2);
    for k = 1:size(seg,1)
        if Lg(k) < 1e-9, continue; end
        n  = max(1, round(Lg(k)/ds));
        t  = ((0:n-1)+0.5)/n;
        pts= P(k,:) + t(:)*seg(k,:);
        nk = [seg(k,2), -seg(k,1)]/Lg(k);                    % 邊的一側法線
        pr = pts(1,:) + 0.02*nk;                             % 探測：若落在多邊形內就翻向
        if inpolygon(pr(1), pr(2), px, pz), nk = -nk; end
        for dd = offs(:).'
            Q = [Q; pts + dd*nk];                            %#ok<AGROW>
        end
    end
    Fx = scatteredInterpolant(SL.x, SL.z, SL.bx, 'linear','none');
    Fz = scatteredInterpolant(SL.x, SL.z, SL.bz, 'linear','none');
    Fm = scatteredInterpolant(SL.x, SL.z, SL.bm, 'linear','none');
    bx = Fx(Q(:,1), Q(:,2));  bz = Fz(Q(:,1), Q(:,2));  bm = Fm(Q(:,1), Q(:,2));
    ok = ~isnan(bx) & ~isnan(bz) & ~isnan(bm) & bm > 1e-4;
    % 去掉落在磁極內部的查詢點（內插會拿到鋼內場）
    ok = ok & ~inpolygon(Q(:,1), Q(:,2), px, pz);
    Qx = Q(ok,1);  Qz = Q(ok,2);  bx = bx(ok);  bz = bz(ok);  bm = bm(ok);
    bxz = hypot(bx,bz);  bxz(bxz==0) = 1e-12;
    scl = SL.arrow_max .* (bm./SL.bmax).^0.25 ./ bxz;
    Uq = bx.*scl;  Wq = bz.*scl;  Bm_mT = bm*1e3;
end

% ============================================================================
function [px, pz] = pole_outline(c, pidx, rf, beta, L)
% 磁極 y=0 截面外框(x,z)：下極水平半切平頂 / 上極傾斜全錐（龍飛模型）。
    th = c.pole_angles(pidx)*pi/180;
    if c.pole_is_lower(pidx)
        a=[cos(th);sin(th);0];  v=[0;0;-1];  half=true;
    else
        inc=c.upper_incline;  a=[cos(inc)*cos(th);cos(inc)*sin(th);sin(inc)];
        u=cross(a,[0;0;1]); u=u/norm(u);  v=cross(a,u);  half=false;
    end
    tip=[c.pole_tip_x(pidx); c.pole_tip_y(pidx); c.pole_tip_z_wp(pidx)]*1e3;
    ts=rf*(1+sin(beta));  rt=rf*cos(beta);
    psi=linspace(0,pi/2+beta,16);  nt=30;  t=linspace(ts,L,nt);
    axoff=[rf*(1-cos(psi)),  t(2:end)];
    radm =[rf*sin(psi),      rt+(t(2:end)-ts)*tan(beta)];
    edgeP = tip + a.*axoff + v.*radm;
    if half
        axline = tip + a.*axoff;  poly=[axline, fliplr(edgeP)];
    else
        edgeM = tip + a.*axoff - v.*radm;  poly=[edgeP, fliplr(edgeM)];
    end
    px=poly(1,:);  pz=poly(3,:);
end

% ============================================================================
function draw_sensor2d(ax, sp, n, R, H, col)
% sensor 圓柱側視矩形(x-z)：底面在 sp、沿 +n 長 H、半寬 R。綠填黑邊。
    sp=sp(:);  n=n(:)/norm(n);  td=[-n(2);n(1)];
    base=sp;  top=sp+H*n;
    rc=[base+R*td, base-R*td, top-R*td, top+R*td];
    patch(ax, rc(1,:), rc(2,:), col, 'EdgeColor','k', 'LineWidth',1.6);
end

% ============================================================================
function draw_narrow2d(ax, sp, n, Ln, col, lw)
% 紅 n+ 箭頭(x-z)：桿 + 小三角頭（對齊 n）。
    sp=sp(:);  n=n(:)/norm(n);  tip=sp+Ln*n;  td=[-n(2);n(1)];
    hl=0.32*Ln;  a=24*pi/180;
    d1=cos(a)*(-n)+sin(a)*td;  d2=cos(a)*(-n)-sin(a)*td;
    plot(ax,[sp(1) tip(1)],[sp(2) tip(2)],'-','Color',col,'LineWidth',lw);
    plot(ax,[tip(1) tip(1)+d1(1)*hl],[tip(2) tip(2)+d1(2)*hl],'-','Color',col,'LineWidth',lw);
    plot(ax,[tip(1) tip(1)+d2(1)*hl],[tip(2) tip(2)+d2(2)*hl],'-','Color',col,'LineWidth',lw);
    % [ADDED 2026-08-14] 箭頭旁標 n+：沿桿 0.62Ln、再往側向 (td) offset 0.30Ln，
    %   放在桿旁而不是頂端外側 —— 兩根極的箭頭都指向框邊，標在外側會貼到框線。
    lab = sp + 0.35*Ln*n + 0.60*Ln*td;
    text(ax, lab(1), lab(2), '$\mathbf{n^{+}}$', 'Interpreter','latex', ...
         'FontSize',30, 'Color',col, 'HorizontalAlignment','center', ...
         'VerticalAlignment','middle', 'Clipping','off');
end

% ============================================================================
function style_cbar(cb, FS)
% colorbar 標準樣式（同 circuit_side）：軸標題 LaTeX 數學粗體、刻度數字粗體。
    cb.FontSize=FS;  cb.FontWeight='bold';
    % [MODIFIED 2026-08-14] 場量標籤改 ‖b‖（figure-style 2026-08-05 拍板；source of truth
    %   = plot_circuit_side.m 的 style_cbar）。原本是 |B|。
    cb.Label.Interpreter='latex';  cb.Label.String='$\mathbf{\|b\|\;(mT)}$';  cb.Label.FontSize=FS;
end
