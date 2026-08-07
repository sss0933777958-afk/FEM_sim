function out = plot_surface_flux(SURF, VARIANT, delta, h_cell, s_sec, POLE, COIL, SRC, WANT_CELLS)
%   [ADDED 2026-08-05] WANT_CELLS=true → out 另回 .P（格心）、.N（外法線）、.corners（n×3×4 格子四角，
%     順序可直接餵 patch）。給「取樣示意圖」用，確保示意圖與實際積分是**同一份幾何**。
%   [ADDED 2026-08-05] SRC = 'apdl'（預設，graded 貼合網格）| 'maxwell'
%     'maxwell' 改讀 D:\Maxwell_sim\...\export\B_voltage_p<POLE>.fld（**唯一涵蓋整根磁極**的匯出，
%     0.1mm 規則格；WP 細格框只有 ±0.6mm）。Maxwell 分支 coil→pole 是 identity → 直接用 POLE 編號。
%     ⚠⚠ 依 .claude/rules/surface-integral-vs-point-sampling.md：**規則直角格不適合做表面積分**
%        （格點不落在傾斜錐面上、每個格心都要跨鐵/空氣界面內插，B 切向分量在界面不連續 → 抹平）。
%        此分支僅供對照，數值有系統性偏差，引用前必須與 APDL 版比對。
%PLOT_SURFACE_FLUX  六極單一磁極側表面的漏磁通量 Phi = ∮ B·n̂ dA（正方形格積分）。
%   out = plot_surface_flux('flat',      'graded')                 % P1 削平面（coil1 自激）
%   out = plot_surface_flux('cone_upper','graded',[],[],[],2,'coil5')  % P2 上半錐面（coil5 = P2 自激）
%   out = plot_surface_flux(SURF, VARIANT, delta, h_cell, s_sec, POLE, COIL)
%
%   SURF：'flat'（削平面，僅半切下極有）| 'cone_upper' | 'cone_lower' | 'section'
%   POLE：1 = P1（下極、方位 0°、錐軸水平）| 2 = P2（上極、方位 180°、錐軸傾 36.59°）
%   COIL：資料夾/檔名前綴（自激發請用該極對應的 coil；long2016 map [1,3,6,5,2,4] → P1=coil1、P2=coil5）
%
%   局部基底（每極自己的）：axis = 錐軸（指向根部）、w = ŷ、up = 垂直 axis 且朝上。
%     表面點 = tip + s·axis + R(s)·(cos φ·w + sin φ·up)；φ>0 = 上半、φ<0 = 下半。
%     R(s) = RT + SLP·s（至 s = S_CONE 後封頂 RMAX，轉等徑圓柱）；鋼件至 s = S_STEEL。
%     ⚠ 幾何常數（RT/SLP/S_CONE/RMAX）兩極共用——已用節點雲分別驗證。
%     半切（VARIANT='graded' 且 POLE=1）：鋼件只在 (P−tip)·up ≤ 0 那半。
%
%   正方形格：
%     'flat'        平面 (s,u)，格心 (k+0.5)h → 格線落在 h 整數倍（切面精確切齊）
%     'cone_*'      圓錐側面是**可展面** → 展開成平面扇形後鋪正方形（每格面積精確 = h²、
%                   曲面上局部為正方形），再映射回 3D。半角 alpha = atan(SLP)。
%     'section'     垂直錐軸的截面 @ s_sec，n̂ = +axis（供散度定理檢查）
%
%   ⚠ 內插聲明：格心的 B 由 FEM **鋼件節點** scatteredInterpolant('linear','nearest') 內插，
%     非節點原值。delta = 沿內法線退入鋼件的距離，**預設 0**（面上本就有節點層；
%     實測往內縮是系統性偏高：0.05mm → +7.8%、0.10mm → +12.9%）。
%
%   單位：Phi[uWb] = B[T]·A[mm²]。號誌：n̂ 取離開鋼件向外為正（回傳為 raw FEM 號誌，
%   下游要 all-source 需自行乘 s_source）。

    if nargin < 1 || isempty(SURF),    SURF    = 'flat';    end
    if nargin < 2 || isempty(VARIANT), VARIANT = 'graded';  end
    if nargin < 3 || isempty(delta),   delta   = 0;         end
    if nargin < 4 || isempty(h_cell),  h_cell  = 0.02;      end
    if nargin < 5 || isempty(s_sec),   s_sec   = 14.827;    end   % 錐底（局部 s）
    if nargin < 6 || isempty(POLE),    POLE    = 1;         end
    if nargin < 7 || isempty(COIL),    COIL    = 'coil1';   end
    if nargin < 8 || isempty(SRC),     SRC     = 'apdl';    end   % [ADDED 2026-08-05]
    if nargin < 9 || isempty(WANT_CELLS), WANT_CELLS = false; end  % [ADDED] 回傳格子四角(給示意圖用)

    ROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
    g = geo_pole(POLE);
    alpha = atan(g.SLP);
    isHalf = (POLE == 1) && ~strcmp(VARIANT,'lower_filled');   % 半切只發生在下極、且非 filled 變體

    %% ---- 載入 + 建鋼件內插器（只用該極錐體內節點）----
    if strcmpi(SRC,'maxwell')
        [x,y,z,B] = load_maxwell_B(POLE);                       % [ADDED] 0.1mm 規則格（見檔頭警告）
        srctag = sprintf('maxwell p%d', POLE);
    else
        [x,y,z,B] = load_node_B(fullfile(ROOT, VARIANT, COIL), COIL);
        srctag = sprintf('%s %s', VARIANT, COIL);
    end
    [s, r, up_c] = local_coord([x y z], g);
    in = s >= 0 & s <= g.S_STEEL & r <= Rloc(s,g)*1.02;
    if isHalf, in = in & up_c <= 0.001; end
    fprintf('%-20s P%d | 節點/格點 %d, 該極鋼件 %d\n', srctag, POLE, numel(x), sum(in));
    assert(sum(in) > 1000, '錐體鋼件節點過少（%d）—— 幾何/variant/coil 可能不對', sum(in));

    F  = scatteredInterpolant(x(in),y(in),z(in),B(in,1),'linear','nearest');
    Fy = F;  Fy.Values = B(in,2);
    Fz = F;  Fz.Values = B(in,3);

    %% ---- 建面 + 格心 + 外法線（全部在局部座標，最後映射回 3D）----
    switch SURF
        case 'flat'      % 削平面：含錐軸的平面，域 = {|u| <= R(s), 0<=s<=S_CONE}，n̂ = up
            k = @(a,b) (floor(a/h_cell)+0.5 : 1 : ceil(b/h_cell)-0.5)*h_cell;
            [Sg, Ug] = meshgrid(k(0,g.S_CONE), k(-g.RMAX-0.1, g.RMAX+0.1));
            keep = abs(Ug) <= Rloc(Sg,g);
            P = g.tip + Sg(keep)*g.ax + Ug(keep)*g.w;
            N = repmat(g.up, sum(keep(:)), 1);
            A_ana = 2*integral(@(t) Rloc(t,g), 0, g.S_CONE);
            if WANT_CELLS                                   % [ADDED] 平面格四角（(s,u) 各 ±h/2）
                h2 = h_cell/2;  sc0 = Sg(keep);  uc0 = Ug(keep);
                cn = zeros(numel(sc0),3,4);  dd = [ +1 +1; -1 +1; -1 -1; +1 -1 ];
                for q = 1:4
                    cn(:,:,q) = g.tip + (sc0+dd(q,1)*h2)*g.ax + (uc0+dd(q,2)*h2)*g.w;
                end
            end

        case {'cone_upper','cone_lower'}   % 半個錐側面（可展面展開鋪正方形）
            upflag = strcmp(SURF,'cone_upper');
            s1 = g.S_CONE/cos(alpha);  psi_max = pi*sin(alpha);
            k = @(a,b) (floor(a/h_cell)+0.5 : 1 : ceil(b/h_cell)-0.5)*h_cell;
            gg = k(-s1, s1);
            [Pq,Qq] = meshgrid(gg,gg);
            sl = hypot(Pq,Qq);  psi = atan2(Qq,Pq);
            if upflag, keep = sl<=s1 & psi>=0 & psi<=psi_max;
            else,      keep = sl<=s1 & psi<=0 & psi>=-psi_max;  end
            sv = sl(keep);  phv = psi(keep)/sin(alpha);
            sc_ = sv*cos(alpha);   rv = sv*sin(alpha);
            rad = cos(phv)*g.w + sin(phv)*g.up;                    % 徑向單位向量
            P = g.tip + sc_*g.ax + rv.*rad;
            N = -sin(alpha)*repmat(g.ax,numel(sv),1) + cos(alpha)*rad;
            A_ana = 0.5*psi_max*s1^2;
            if WANT_CELLS                                   % [ADDED] 展開面上的正方形四角 → 映射回 3D
                h2 = h_cell/2;  pq0 = Pq(keep);  qq0 = Qq(keep);
                cn = zeros(numel(pq0),3,4);  dd = [ +1 +1; -1 +1; -1 -1; +1 -1 ];
                for q = 1:4
                    slq = hypot(pq0+dd(q,1)*h2, qq0+dd(q,2)*h2);
                    phq = atan2(qq0+dd(q,2)*h2, pq0+dd(q,1)*h2)/sin(alpha);
                    cn(:,:,q) = g.tip + (slq*cos(alpha))*g.ax + (slq*sin(alpha)).*(cos(phq)*g.w + sin(phq)*g.up);
                end
            end

        case 'section'   % 垂直錐軸截面 @ s_sec，n̂ = +axis
            Rc = Rloc(s_sec,g);
            k = @(a,b) (floor(a/h_cell)+0.5 : 1 : ceil(b/h_cell)-0.5)*h_cell;
            gg = k(-Rc,Rc);
            [Wg,Ug] = meshgrid(gg,gg);
            keep = hypot(Wg,Ug) <= Rc;
            if isHalf, keep = keep & Ug < 0; end
            P = g.tip + s_sec*g.ax + Wg(keep)*g.w + Ug(keep)*g.up;
            N = repmat(g.ax, sum(keep(:)), 1);
            A_ana = pi*Rc^2;  if isHalf, A_ana = A_ana/2; end

        otherwise
            error('SURF 必為 ''flat'' | ''cone_upper'' | ''cone_lower'' | ''section''');
    end

    A_grid = size(P,1)*h_cell^2;
    Q  = P - delta*N;
    bx = F(Q(:,1),Q(:,2),Q(:,3));  by = Fy(Q(:,1),Q(:,2),Q(:,3));  bz = Fz(Q(:,1),Q(:,2),Q(:,3));
    Bn  = sum([bx by bz].*N, 2);
    Phi = sum(Bn)*h_cell^2;

    sloc = (P - g.tip) * g.ax.';          % 每格的局部軸向位置（供累積曲線）
    out = struct('SURF',SURF,'VARIANT',VARIANT,'POLE',POLE,'COIL',COIL, ...
                 'delta',delta,'h',h_cell,'Phi',Phi,'ncell',size(P,1), ...
                 'A_grid',A_grid,'A_ana',A_ana,'Bn_mean',mean(Bn), ...
                 's_cell',sloc,'dPhi',Bn*h_cell^2,'S_CONE',g.S_CONE);
    if WANT_CELLS                                          % [ADDED] 示意圖用（格心/法線/四角）
        out.P = P;  out.N = N;  out.Bn = Bn;  out.Bvec = [bx by bz];   % [ADDED] 供 |B|/夾角拆解
        if exist('cn','var'), out.corners = cn; end
        out.tip = g.tip;  out.ax = g.ax;  out.up = g.up;  out.w = g.w;
    end

    fprintf(['  %-10s delta=%.3f h=%.3f | cells=%7d  A=%8.4f/%8.4f mm² (%+.3f%%) | ' ...
             'Phi=%+.5f uWb  <B·n>=%+.4f mT\n'], SURF, delta, h_cell, out.ncell, ...
            A_grid, A_ana, (A_grid/A_ana-1)*100, Phi, mean(Bn)*1e3);
end

% ================= local =================

function g = geo_pole(POLE)
% 每極的局部基底與錐體常數（tip/axis 取自 mt_constants；已用節點雲驗證）
    g.RT = 0.04;  g.SLP = 0.2028;  g.S_CONE = 14.827;  g.S_STEEL = 17.5;
    g.RMAX = g.RT + g.SLP*g.S_CONE;                      % 3.047 mm
    if POLE == 1
        g.tip = [ 0.4082 0 -13.0000];   g.ax = [ 1 0 0];
    else
        inc = 36.5895*pi/180;
        g.tip = [-0.4082 0 -12.4226];   g.ax = [-cos(inc) 0 sin(inc)];
    end
    g.w  = [0 1 0];                                       % 水平橫向（兩極皆在 xz 平面）
    g.up = cross(g.ax, g.w);  g.up = g.up/norm(g.up);     % 垂直 axis、朝上
    if g.up(3) < 0, g.up = -g.up; end
end

function [s, r, up_c] = local_coord(P, g)
    d = P - g.tip;
    s = d*g.ax.';
    up_c = d*g.up.';
    r = vecnorm(d - s*g.ax, 2, 2);
end

function R = Rloc(s, g)
    R = min(g.RT + g.SLP*s, g.RMAX);   R = max(R, 0);
end

function [x,y,z,B] = load_maxwell_B(POLE)
% [ADDED 2026-08-05] Maxwell sensor 粗格匯出（0.1mm、x±15.23 y±14 z[-16.1,0] mm、13.88M 格點）——
%   **唯一涵蓋整根磁極**的 Maxwell 匯出（WP 細格框只有 ±0.6mm）。coil→pole = identity。
%   同一顆極會被呼叫多次（不同面）→ persistent 快取，避免重讀 1.97 GB。
%   回傳與 load_node_B 同介面：座標 mm（raw/Maxwell frame，與 APDL 同框）、B [T]。
    persistent CACHE
    if isempty(CACHE), CACHE = containers.Map('KeyType','double','ValueType','any'); end
    if ~CACHE.isKey(POLE)
        addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell\function');
        f = sprintf('D:\\Maxwell_sim\\long2016_hexapole_halfcut\\export\\B_voltage_p%d.fld', POLE);
        d = import_maxwell_fld(f);
        CACHE(POLE) = struct('x',d.x*1e3, 'y',d.y*1e3, 'z',d.z*1e3, 'B',[d.bx d.by d.bz]);
    end
    S = CACHE(POLE);  x = S.x;  y = S.y;  z = S.z;  B = S.B;
end

function [x,y,z,B] = load_node_B(dir_coil, cname)
    C = readmatrix(fullfile(dir_coil,[cname '_coord_all.dat']),'FileType','text');
    C = C(~any(isnan(C(:,1:4)),2),:);
    t = fileread(fullfile(dir_coil,[cname '_bfield_all.dat']));
    t = regexprep(t,'(\d)([-+])','$1 $2');
    tf=[tempname '.dat']; f=fopen(tf,'w'); fprintf(f,'%s',t); fclose(f);
    M = readmatrix(tf,'FileType','text'); delete(tf);
    M = M(~any(isnan(M(:,1:5)),2),:);
    [~,ic,ib] = intersect(C(:,1),M(:,1));
    x=C(ic,2)*1e3; y=C(ic,3)*1e3; z=C(ic,4)*1e3;  B=M(ib,2:4);
end
