%% calib_fem.m
%  ============================================================================
%  依「Calibration using FEM modeling_V2.pdf」實作【步驟 1–7】。
%  步驟 1–7 的目標:建立純量成本函數 J(ell_hat)。其中每次模擬 j 的 6x1 電荷向量
%  g_j 已用最小二乘「剖面化(profile out)」閉形式解掉,使成本只剩單一非線性參數
%  ell_hat(特徵長度)。步驟 8(最小化 J)與步驟 9–11(回推 g_j、G、g_B*K_I)本檔尚未做。
%
%  文件對應到本專案資料的地方以 >>確認<< 標出,請使用者拍板。每行皆有中文註解。
%  ============================================================================
clear; clc;                                                  % 清空工作區與命令視窗

% ---- 路徑 ------------------------------------------------------------------
% mt_constants / import_ansys_data / filter_iron_nodes 放在 long2016 的 analysis 夾
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');   % 加入共用工具路徑
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';  % FEM 結果(coil1..coil6)根目錄

% ---- 頂層參數 --------------------------------------------------------------
cnst       = mt_constants();        % 幾何常數結構(R_norm、極尖位置、SPH_OFST 等)
R_select   = 150e-6;                % >>確認<< 步驟1 取點半徑 [m]:只保留 |p|<=R_select 的 air 節點
ell_design = cnst.R_norm;           % >>確認<< 步驟5 ell_hat 初始猜 [m] = 尖端到 WP 距離(500 um)
N_I        = 6;                     % FEM 模擬次數 = 6 個單線圈解(見步驟2)
apdl_to_paper_idx = [1,3,6,5,2,4];  % APDL coil j → paper pole 索引(coil1→P1, coil2→P3, ...)

%% ===========================================================================
%  步驟 1 - 在工作空間選 N_p 個點 p_i(3x1)。
%  本實作:用落在 WP 周圍 R_select 球內的 FEM air 節點當這些點。6 顆 coil 解共用
%  同一個 mesh,所以節點集合與順序在 6 顆之間完全一致(既有腳本已驗證)。
%  ===========================================================================
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];  % 3x6 極尖位置(WP 座標)[m]
dhat = tip ./ vecnorm(tip);                                    % 3x6 極尖「單位方向」(每欄各自正規化)

d1   = import_ansys_data(fullfile(results_root, 'coil1', 'standard'),'wp','coil1');  % 載 coil1 的 'wp' 資料(x,y,z,bx,by,bz)
air1 = filter_iron_nodes(d1.x,d1.y,d1.z,cnst,struct('visualize',false));% 邏輯遮罩:TRUE=air 節點(已排除鐵芯)
zwp1 = d1.z - cnst.SPH_OFST;                                    % 把 ANSYS z 平移到 WP 為原點的 z
P_all = [d1.x(air1), d1.y(air1), zwp1(air1)];                   % (Nair x 3) 所有 air 節點座標(WP 座標)[m]

insel = sum(P_all.^2,2) <= R_select^2;                         % 邏輯:半徑 |p| <= R_select 的節點
P     = P_all(insel,:);                                        % (N_p x 3) 選出來的工作點 p_i
Np    = size(P,1);                                             % N_p = 選出的點數
fprintf('步驟1: 在 R_select = %g um 內選出 N_p = %d 個點\n', R_select*1e6, Np);  % 回報點數

%% ===========================================================================
%  步驟 2 - 做 N_I 次 FEM,輸入電流向量 I_j(6x1);合併電流矩陣
%  F = [I_1 ... I_{N_I}] 須 rank = 6。
%  本實作:N_I = 6 個單線圈解就是這些模擬。第 j 次模擬「依序激發第 j 顆 APDL coil」
%  通 1A,對應 paper pole = apdl_to_paper_idx(j),故電流向量(pole 空間)
%  I_j = e_{pj(j)}。F = [I_1 ... I_6] 為置換矩陣(rank 6)。 >>確認<<
%  ===========================================================================
F = zeros(6, N_I);                                            % 6x6 電流矩陣(pole 空間,P1..P6)
for j = 1:N_I                                                 % 逐次模擬
    F(apdl_to_paper_idx(j), j) = 1;                           % 第 j 次:paper pole pj(j) 通 1A,其餘 0
end
assert(rank(F) == 6, '步驟2: 電流矩陣 F 必須 rank = 6');        % 強制滿足 rank 6 的要求

%% ===========================================================================
%  步驟 3 - 記錄磁通密度 b_ij = b^FEM(p_i, I_j)(3x1)。
%  單線圈模擬下,b_ij = 第 j 顆 coil 在點 i 的場。存成 B(Np,3,N_I)。
%  符號慣例:取 Bn = -B^FEM(整體變號),與本專案 all-source 流程(產生
%  fit_KI_R150.tex 的 sweep)一致,如此步驟 11 的 gB 為正、K̄_I gauge 對齊、
%  可逐數比較。J(ell) 與 ell_hat 對此整體變號「不變」(b 出現兩次),只影響
%  步驟 9–11 的 g_j/K̄_I 符號。
%  ===========================================================================
B = zeros(Np,3,N_I);                                          % 預配置 b_ij 容器(點, xyz, 模擬)
for j = 1:N_I                                                  % 逐一掃 6 個單線圈模擬
    if j == 1                                                 % coil1 在步驟1 已載入 → 直接重用,不重載
        dj = d1;  airj = air1;                                % 重用步驟1 的 coil1 資料與 air 遮罩
    else                                                      % coil2..6 才新載入
        cn = sprintf('coil%d', j);                            % 目錄/工作名 'coil2'..'coil6'
        dj = import_ansys_data(fullfile(results_root, cn, 'standard'),'wp',cn);  % 載 coil j 的 'wp' 資料
        airj = filter_iron_nodes(dj.x,dj.y,dj.z,cnst,struct('visualize',false));  % air 遮罩(幾何相同→遮罩相同)
    end
    Bj_all = -[dj.bx(airj), dj.by(airj), dj.bz(airj)];        % (Nair x 3) coil j 在所有 air 節點的場 [T];取負(all-source 慣例)
    B(:,:,j) = Bj_all(insel,:);                               % 只留 R_select 內選出的點 → b_ij
end
fprintf('步驟3: 已記錄 b_ij(N_I = %d 個模擬,每個 %d 點)\n', N_I, Np);  % 回報

%% ===========================================================================
%  步驟 4 - 記錄感測器電壓 v_jk(k=1..6, j=1..N_I)。
%  步驟 7 的成本 J(ell) 用不到 v_jk:它只在 page2 的霍爾感測器模型/模型參數估計
%  (步驟 7 以後)才會用到。這裡先留占位,等做步驟 9+ 再載入。 >>確認<<
%  ===========================================================================
v = [];                                                       % v_jk(6 x N_I)— 暫緩;J(ell) 不需要

%% ===========================================================================
%  步驟 5 - 引入特徵長度 ell_hat(初始猜 = ell_design)。
%  ===========================================================================
ell = ell_design;                                             % 純量 ell_hat [m];唯一的非線性未知數

%% ===========================================================================
%  步驟 6 - 非量綱 3x6 空間函數矩陣 S(pbar_i; ell)。
%  電荷位置矩陣 Pc(3x6):文件理想型把 6 顆電荷放在 ±x,±y,±z(致動器座標)。
%  這裡 Pc = dhat = FEM/量測座標下「實際的 6 極尖端單位方向」(3 組互相正交的
%  對極 = 理想型的旋轉版),如此點可留在 FEM 座標、免做旋轉。 >>確認<<
%  S 第 k 欄 = (pbar_i - pbar_ck)/||pbar_i - pbar_ck||^3,pbar = p/ell。
%  (build_S 定義在檔案最下方的本地函式。)
%  ===========================================================================
Pc = dhat;                                                    % 3x6 電荷(極)位置,放在單位距離

%% ===========================================================================
%  步驟 7 - 成本函數 J(ell):
%    J(ell) = sum_j { sum_i b_ij^T b_ij
%                     - (sum_i b_ij^T S_i)(sum_i S_i^T S_i)^{-1}(sum_i S_i^T b_ij) }
%  中括號 = 模擬 j 把它的最佳電荷向量 g_j = (sum S^T S)^{-1}(sum S^T b_ij)
%  剖面化解掉後的殘差平方和。sum_i S_i^T S_i(6x6)與 j 無關,每個 ell 只算一次。
%  (costJ 為檔案最下方的本地函式。)
%  ===========================================================================
J0 = costJ(ell, P, B, Pc);                                    % 在初始猜評估一次 J(sanity check)
fprintf('步驟7: J(ell_design = %.1f um) = %.6e  (單位 T^2,預期有限正值)\n', ell*1e6, J0);  % 回報

%% ===========================================================================
%  步驟 8 - 最小化 J(ell) ⟹ ell_hat(初始猜 = ell_design)。
%  J 只是單一純量 ell 的函數(g_j 已被剖面化解掉),用一維 fminbnd 找最小值。
%  ===========================================================================
obj    = @(x) costJ(x, P, B, Pc);                            % 目標函式:成本對 ell 的一維函數
ell_lo = 0.2e-3;  ell_hi = 2.0e-3;                           % ell 搜尋區間 [0.2, 2.0] mm(實體合理範圍)
fopt   = optimset('TolX',1e-9,'Display','iter');             % 收斂容差 1 nm,顯示每次迭代
[ell_hat, Jmin] = fminbnd(obj, ell_lo, ell_hi, fopt);        % 一維最小化 → 最佳 ell_hat 與最小成本 Jmin
fprintf('步驟8: ell_hat = %.4f mm  (Jmin = %.6e),初始猜 ell_design = %.3f mm\n', ...
        ell_hat*1e3, Jmin, ell_design*1e3);                 % 回報 ell_hat(預期 ≈ fit_KI_R150.tex 的 0.856 mm)

%% ===========================================================================
%  步驟 9 - 在最佳 ell_hat 算各模擬電荷向量
%           g_j = (Σ_i S^T S)^{-1}(Σ_i S^T b_ij)(6x1)。
%  在 ell_hat 下重建法矩陣 M = Σ S_i^T S_i 與右端 c_j = Σ S_i^T b_ij(同步驟7)。
%  ===========================================================================
M = zeros(6,6);                                              % Σ_i S_i^T S_i(6x6,與 j 無關)
c = zeros(6,N_I);                                            % Σ_i S_i^T b_ij(每模擬 6x1)
for i = 1:Np                                                 % 逐點累加
    Si = build_S(P(i,:), ell_hat, Pc);                      % 該點 3x6(用最佳 ell_hat)
    M  = M + Si.' * Si;                                      % 累加 S_i^T S_i
    for j = 1:N_I                                            % 逐模擬
        c(:,j) = c(:,j) + Si.' * squeeze(B(i,:,j)).';       % 累加 S_i^T b_ij
    end
end
G = M \ c;                                                  % 6xN_I:各欄 g_j = M^{-1} c_j(步驟9 的 g_j)

%% ===========================================================================
%  步驟 10 - G = [g_1 ... g_{N_I}](6 x N_I)。上式 G 即是。
%  ===========================================================================
% (G 已於步驟9 取得)

%% ===========================================================================
%  步驟 11 - 求 g_B*K̄_I:
%    H = G F^T (F F^T)^{-1}        (6x6,= gB_hat * K̄_I)
%    gB_hat = (6/5) h11            (因 gauge k̄_I(1,1)=5/6 → H(1,1)=gB*5/6)
%    K̄_I    = 5/(6 h11) H          (= H/gB_hat,其 (1,1) 自動 = 5/6)
%  註:文件第二式原寫 "5/(6 g11) H" 為筆誤,正解用 H 的 (1,1) = h11(非 G 的 g11);
%      兩式一致才會使 K̄_I(1,1)=5/6。本例因 coil1 激發 P1(置換不動第1欄),h11=g11。
%  ===========================================================================
H      = G * F.' / (F * F.');                               % 6x6:gB_hat * K̄_I
h11    = H(1,1);                                             % H 的 (1,1) 元素 = h11
gB_hat = (6/5) * h11;                                        % ĝ_B = (6/5) h11
KbarI  = (5/(6*h11)) * H;                                    % K̄_I = 5/(6 h11) H(K̄_I(1,1)=5/6)

% all-source 顯示:翻上極 P2/P4/P5 欄(純符號重貼,不改場;對齊 fit_KI_R150.tex)
coil_sign  = [1 -1 1 -1 -1 1];                              % paper P1..P6(下極+1、上極-1)
KbarI_disp = KbarI .* coil_sign;                           % 顯示用 K̄_I(對角全正)

%% ===========================================================================
%  與 fit_KI_R150.tex 比較(該檔由昨天全節點 fit 產生:ell=0.856mm, gB=8.429e-3)
%  ===========================================================================
Kref = [ 0.8333 -0.0822 -0.1857 -0.1721 -0.1751 -0.2005;    % fit_KI_R150.tex 的 K̂_I(all-source)
        -0.1947  0.9782 -0.2254 -0.1460 -0.1476 -0.2353;
        -0.1896 -0.1707  0.8318 -0.0801 -0.1712 -0.2002;
        -0.2226 -0.1476 -0.1911  0.9800 -0.1490 -0.2321;
        -0.2247 -0.1474 -0.2214 -0.1463  0.9818 -0.2048;
        -0.1892 -0.1724 -0.1860 -0.1714 -0.0839  0.8311];
ell_ref = 0.856e-3;  gB_ref = 8.429e-3;                     % 參考的 ell、gB

fprintf('\n=============== 步驟 9–11 結果 vs fit_KI_R150.tex ===============\n');
fprintf('  ell_hat = %.4f mm   (參考 %.3f mm,差 %.4f mm)\n', ell_hat*1e3, ell_ref*1e3, (ell_hat-ell_ref)*1e3);
fprintf('  gB_hat  = %.4e      (參考 %.4e,相對差 %.2f%%)\n', gB_hat, gB_ref, (gB_hat-gB_ref)/gB_ref*100);
fprintf('  K̄_I (all-source) =\n');
for i=1:6, fprintf('    % .4f % .4f % .4f % .4f % .4f % .4f\n', KbarI_disp(i,:)); end
fprintf('  max|K̄_I - K_ref| = %.4f\n', max(abs(KbarI_disp(:)-Kref(:))));  % 與參考矩陣的最大逐元素差
fprintf('================================================================\n');

%% ===========================================================================
%  PAGE 2 - Hall-sensor-based 模型 + 模型參數估計 d
%  模型:b_ij = Σ_k d_k v_jk s_ik = S_i V_j d(電荷綁到 sensor 讀數 v_jk)。
%  解  :d = (Σ_j V_j M V_j)^{-1}(Σ_j V_j cc_j),M=Σ_i S_i^T S_i、cc_j=Σ_i S_i^T b_ij
%        皆重用步驟 9 在 ell_hat 下算好的 M 與 c(:,j)。ell_hat 沿用步驟 8(R<=150um)。
%  v_jk:依兩張圖的 sensor 擺放(= 既有 B_bar 慣例)抽取。
%  符號(per feedback_sensor_sign_convention_toward_wp):取「n+ 方向上的分量」,
%        磁場朝 n+ → 正、朝 n− → 負,即物理 signed B·n+(不壓 cosmetic 號;off-diag 本就混合)。
%        d-fit 的 b 也用物理場(cc_phys = −c,c 來自頁1 負號版 B)→ b、v 一致;
%        d 對 b、v 的整體同號變換不變,故 d 數值與先前相同。
%  ===========================================================================
S_hall = 130;                                                % Hall 靈敏度 [V/T](EQ-730L:130 V/T,非舊版 130e-3)

% ---- sensor 幾何(下極=milled flat;上極=natural cone;n+ 出鋼)----  >>確認<< 沿用 B_bar 慣例
beta = atan2(3.0, 15.0);                                     % 錐半角 11.31°(ANSYS POLE_R 3.0 / CONE_LEN 15.0)
sensor_pl_lower = [4.572e-3; 0.41e-3; 0];                    % 下極 pole-local:4.572 沿軸、0.41 ⊥ milled flat
n_pl_lower      = [0; 1; 0];                                 % 下極法線 n+ = +up(⊥ 平面、出鋼)
sensor_pl_upper = [4.572e-3*cos(beta) - 0.41e-3*sin(beta); ...% 上極 pole-local:4.572 沿錐斜邊
                   4.572e-3*sin(beta) + 0.41e-3*cos(beta); ...%   + 0.41 ⊥ 錐面外
                   0];
n_pl_upper      = [-sin(beta); cos(beta); 0];               % 上極法線 n+ ⊥ 錐面(出鋼)

sensor_pos = zeros(3,6); sensor_n = zeros(3,6);            % 全域 sensor 位置 / 法線
disc_u = zeros(3,6); disc_v = zeros(3,6);                  % disc 平面兩個基底向量
for i = 1:6                                                % 逐極 P1..P6
    th = cnst.pole_angles(i)*pi/180;                      % 方位角
    if cnst.pole_is_lower(i)                              % --- 下極 ---
        pole_axis = [cos(th); sin(th); 0];               % 極軸 = 水平徑向
        up_hat = [0;0;1]; tip_z = -cnst.R_norm_z;        % up = +z;尖端 z
        s_pl = sensor_pl_lower; n_pl = n_pl_lower;        % 下極 offset / 法線
    else                                                 % --- 上極 ---
        inc = cnst.upper_incline;                        % 上極傾角(~36.6°,沿用 mt_constants)
        pole_axis = [cos(inc)*cos(th); cos(inc)*sin(th); sin(inc)];  % 傾斜極軸
        up_un = [0;0;1] - sin(inc)*pole_axis; up_hat = up_un/norm(up_un);  % up = ẑ 去掉沿軸分量
        tip_z = +cnst.R_norm_z; s_pl = sensor_pl_upper; n_pl = n_pl_upper; % 上極 offset / 法線
    end
    tip  = [cnst.R_norm_xy*cos(th); cnst.R_norm_xy*sin(th); tip_z];  % 極尖座標(WP 框)
    side = cross(pole_axis, up_hat); side = side/norm(side);         % 側向(⊥ axis & up)
    Rg   = [pole_axis, up_hat, side];                    % pole-local → global 旋轉矩陣
    sensor_pos(:,i) = tip + Rg*s_pl;                     % sensor 全域位置
    sensor_n(:,i)   = Rg*n_pl;                           % sensor 全域法線 n+
    u = up_hat - dot(up_hat,sensor_n(:,i))*sensor_n(:,i);% disc 基底 u = up_hat 投到 ⊥n 平面
    if norm(u)<1e-9, u = pole_axis - dot(pole_axis,sensor_n(:,i))*sensor_n(:,i); end  % 退化時改用 axis
    u = u/norm(u); disc_u(:,i) = u; disc_v(:,i) = cross(sensor_n(:,i),u);  % 正交基底 u,v
end

% ---- Ø0.3mm Hall disc 取樣點(同心環,721 點;面積平均近似)----
sensor_r = 0.15e-3; n_rings = 15; ring_dr = sensor_r/n_rings; % 半徑 0.15mm、15 環、環距 0.01mm
disc_local = [0 0];                                          % 中心點
for kk = 1:n_rings                                          % 每環點數 6k(正比半徑)
    rk = ring_dr*kk; nk = 6*kk; phi = (0:nk-1)*2*pi/nk;
    disc_local = [disc_local; rk*cos(phi(:)), rk*sin(phi(:))]; %#ok<AGROW>
end
Ndisc = size(disc_local,1);                                % = 721

% ---- 抽 v_jk:Vmat(i,j) = S_hall * disc 平均((−B)·n_i);列 i=sensor 極 P1..P6、欄 j=coil(sim)----
%  >>確認<< 用 coil1..coil6(v4 baseline,'all' dataset);不可用已刪的 coil1_pre_fine_mesh(RESULTS_MAP)
Vmat = zeros(6,6);                                         % sensor 電壓矩陣 [V]
for kc = 1:6                                               % 逐 coil(= sim j = kc)
    cn = sprintf('coil%d', kc);                           % coil 資料夾名
    da = import_ansys_data(fullfile(results_root, cn, 'standard'),'all',cn);  % 載 'all' 全域(sensor 在遠處錐面)
    zc = da.z - cnst.SPH_OFST;                            % ANSYS z → WP 框
    Fx = scatteredInterpolant(da.x,da.y,zc,da.bx,'linear','nearest');  % Bx 內插器
    Fy = scatteredInterpolant(da.x,da.y,zc,da.by,'linear','nearest');  % By
    Fz = scatteredInterpolant(da.x,da.y,zc,da.bz,'linear','nearest');  % Bz
    for i = 1:6                                           % 逐 sensor 極
        acc = 0;                                          % disc 上 (−B)·n+ 累加
        for p = 1:Ndisc                                  % 逐 disc 取樣點
            xp = sensor_pos(:,i) + disc_local(p,1)*disc_u(:,i) + disc_local(p,2)*disc_v(:,i);  % 取樣點全域座標
            Bv = [Fx(xp(1),xp(2),xp(3)); Fy(xp(1),xp(2),xp(3)); Fz(xp(1),xp(2),xp(3))];        % 該點 B
            acc = acc + Bv.' * sensor_n(:,i);            % 物理 signed B·n+(朝 n+ 正、朝 n− 負)
        end
        Vmat(i,kc) = S_hall * acc / Ndisc;               % 面積平均 × S_hall = sensor 電壓 [V]
    end
    fprintf('Page2: coil%d sensor 電壓抽取完成\n', kc);   % 進度
end

% ---- all-source 慣例:翻「下極激發」(P1/P3/P6)的欄,使每顆激發極都當 source ----
%  物理:翻線圈繞線方向 ⟺ 把該 coil 的 B 整個變號(磁靜場對電流線性)→ 後處理翻號
%  與「反向繞線重跑 FEM」bit-for-bit 等價,免重跑。翻後 self 對角全正、off-diag 幾乎全負,
%  與頁1 K̄_I(all-source)一致。下極 raw 是 sink,翻號後變 source。
exc_sign = ones(1,6);                                    % 各 sim(coil)激發極的 source 翻號
for j = 1:6                                              % 下極 P1/P3/P6 激發 → 翻號變 source
    if ismember(apdl_to_paper_idx(j), [1 3 6]), exc_sign(j) = -1; end
end
Vmat = Vmat .* exc_sign;                                 % all-source sensor 電壓(翻下極激發欄)

% ---- 解 d(依 d_final.pdf,含增益 g_H)----
%  模型:b_ij = g_H · S_i V_j d;閉式解 d = (1/g_H)(Σ_j V_j M V_j)^{-1}(Σ_j V_j cc_j)。
%  g_H = k_m/(ℓ̂²μ_0) = 1/(4πℓ̂²) → 1/g_H = 4πℓ̂²。M、cc 重用步驟9(在 ell_hat 下)。
gH = cnst.k_m / (ell_hat^2 * cnst.mu_0);                  % 增益 g_H = 1/(4πℓ̂²)  [1/m^2]
A = zeros(6,6); rhs = zeros(6,1);                         % 法矩陣與右端累加器
for j = 1:N_I                                             % 逐模擬
    Vj  = diag(Vmat(:,j));                                % V_j = diag(該 sim 的 6 個 all-source sensor 電壓)
    ccj = exc_sign(j) * (-c(:,j));                        % all-source cc_j(b 同步翻:−c=物理,×exc_sign=all-source)
    A   = A   + Vj * M * Vj;                              % Σ V_j M V_j
    rhs = rhs + Vj * ccj;                                 % Σ V_j cc_j
end
d = (A \ rhs) / gH;                                       % d_final = (1/g_H)(ΣVMV)^{-1}(ΣVcc) = 4πℓ̂²·(無增益 d)

% ---- sensor 模型殘差(R<=150 點上:b_ij vs S_i V_j d)----
num = 0; den = 0;                                         % 相對 RMSE 的分子/分母
for i = 1:Np                                              % 逐點
    Si = build_S(P(i,:), ell_hat, Pc);                   % 3x6(在 ell_hat)
    for j = 1:N_I                                         % 逐模擬
        bm = gH * Si * diag(Vmat(:,j)) * d;              % sensor 模型場 g_H·S_i V_j d(含增益,all-source)
        bf = exc_sign(j) * (-squeeze(B(i,:,j)).');       % all-source FEM 場 b_ij(−B=物理,×exc_sign=all-source)
        num = num + sum((bm-bf).^2); den = den + sum(bf.^2);
    end
end
nrmse_sensor = sqrt(num/den)*100;                        % 全域相對 RMSE [%]

fprintf('\n=============== PAGE 2: Hall-sensor 模型 (ell_hat=%.3f mm) ===============\n', ell_hat*1e3);
fprintf('  sensor 電壓 Vmat [V](列=sensor 極 P1..P6,欄=激發 coil1..6):\n');
for i=1:6, fprintf('   % .3e % .3e % .3e % .3e % .3e % .3e\n', Vmat(i,:)); end  % 印 6x6 電壓
fprintf('  d (6x1, P1..P6) [A/(V·...)]:\n'); fprintf('   % .4e\n', d);          % 印 d
fprintf('  sensor 模型在 R<=150um 的相對 RMSE = %.2f%%\n', nrmse_sensor);        % 殘差
fprintf('========================================================================\n');

% ---- 存 page-2 結果(供 gen_d_latex.m 出 LaTeX)----
out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit\fitting_d';
if ~exist(out_dir,'dir'); mkdir(out_dir); end
save(fullfile(out_dir,'calib_sensor_d.mat'), ...
     'd','gH','Vmat','ell_hat','nrmse_sensor','S_hall','sensor_pos','sensor_n','apdl_to_paper_idx','R_select');
fprintf('已存 fitting_d/calib_sensor_d.mat\n');

%% ===========================================================================
%  本地函式
%  ===========================================================================
function S = build_S(p_row, ell, Pc)
% BUILD_S  步驟6 的單點空間函數矩陣 S(pbar; ell)。
%   p_row : 1x3 工作點 p_i [m]
%   ell   : 純量特徵長度 ell_hat [m]
%   Pc    : 3x6 電荷位置矩陣(單位距離的極方向)
%   S     : 3x6,第 k 欄 = (pbar - Pc(:,k)) / ||pbar - Pc(:,k)||^3
    pbar = p_row(:) / ell;                                    % 3x1 正規化點 pbar_i = p_i/ell
    D    = pbar - Pc;                                         % 3x6 各電荷的差 (pbar - pbar_ck)(廣播相減)
    nrm  = vecnorm(D);                                        % 1x6 各欄歐氏範數 ||pbar - pbar_ck||
    S    = D ./ (nrm.^3);                                     % 3x6 每欄除以自己範數的立方 → s_ik
end

function J = costJ(ell, P, B, Pc)
% COSTJ  步驟7 的剖面化最小二乘成本 J(ell),對所有點與所有模擬加總。
%   ell : 純量 ell_hat [m]
%   P   : Np x 3 工作點
%   B   : Np x 3 x N_I 記錄的場 b_ij
%   Pc  : 3 x 6 電荷位置矩陣
%   J   : 純量成本(各模擬剖面化後殘差平方和之總和)
    Np  = size(P,1);                                          % 點數
    N_I = size(B,3);                                          % 模擬數
    M   = zeros(6,6);                                         % 累加器:sum_i S_i^T S_i(6x6,與 j 無關)
    c   = zeros(6,N_I);                                       % 累加器:sum_i S_i^T b_ij(每模擬 6x1)→ 6xN_I
    bb  = zeros(1,N_I);                                       % 累加器:sum_i b_ij^T b_ij(每模擬純量)
    for i = 1:Np                                              % 逐點 p_i
        Si = build_S(P(i,:), ell, Pc);                       % 該點的 3x6 空間矩陣
        M  = M + Si.' * Si;                                   % 把 S_i^T S_i(6x6)加進共用法矩陣
        for j = 1:N_I                                         % 逐模擬
            bij    = squeeze(B(i,:,j)).';                     % 3x1 該點記錄的場 b_ij
            c(:,j) = c(:,j) + Si.' * bij;                     % 把 S_i^T b_ij(6x1)加進模擬 j 的累加器
            bb(j)  = bb(j) + bij.' * bij;                     % 把 b_ij^T b_ij(純量)加進模擬 j 的累加器
        end
    end
    J = 0;                                                    % 初始化總成本
    for j = 1:N_I                                            % 加總各模擬剖面化後的殘差
        J = J + ( bb(j) - c(:,j).' * (M \ c(:,j)) );         % bb - c^T M^{-1} c(M\c 即解 M g = c)
    end
end
