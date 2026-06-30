%% main.m -- Hall-sensor-based hexapole model: 求每極常數 d（指定流程 driver）
% =========================================================================
%  這支是「Hall_sensor_base_fix_dir」包的主程式(driver)。本檔只串流程＋印結果＋存檔，
%  所有模型數學都在 code\function\ 裡。流程嚴格照使用者指定的四步：
%
%  ── 概念 ────────────────────────────────────────────────────────────────
%  把每顆磁極等效成一顆磁荷(放在極尖方向、距 WP 中心 ℓ̂)。Hall-sensor 模型主張：
%  「第 j 次激發、第 k 極的電荷強度 ∝ 該極 sensor 量到的電壓 v_jk × 一個每極常數 d_k」，
%  於是 WP 附近的場可寫成（no-gain，無 g_H）
%        b_ij = S_i · V_j · d                    (V_j = diag(Vmat(:,j)))
%  只有 6 個 d_k 是自由參數(V_j 由 FEM 固定)。本檔把 d 求出來。
%  ("fix_dir" = 電荷鎖在極軸方向 d̂ 上、不離軸。)
%
%  ── 流程（兩階段：先 fit ℓ̂、再求 shared d）────────────────────────────────
%    步驟1  載 6-coil FEM 場（R≤R_select 球內 air 節點）→ 工作點 P、場 B。
%    步驟2  抽電壓 V=S·B（真實 FEM 節點）：build_sensor_geometry 給 sensor 位置/法線；
%                 extract_Vmat 沿法線 n 開圓柱選真實節點、對 B·n 平均 → Vmat（含 all-source 翻號）。
%    步驟3  【先找 ℓ̂】minimize profiled cost J(ℓ̂)（消去 d、per-j；sensor_cost_lhat + fminbnd），
%                 不再從 fix_dir fit_KI_fixl 載 ℓ̂。
%    步驟4  【再找 d】用該 ℓ̂ 閉式解 shared 6×1 d（solve_d；反矩陣在 Σ_j 外）、cost J = Σ‖ε‖²。
%
%  ── 關鍵變數 ────────────────────────────────────────────────────────────
%    P(Np×3)   擬合區(R≤R_select 球)內工作點；B(Np×3×6) 對應場(負號版 −B^FEM)
%    ell_hat   特徵長度 ℓ̂（載入，給 build_S）
%    Vmat(6×6) sensor 電壓(all-source)；d(6×1) 每極常數(no-gain)；J 成本 [T²]
%
%  Current : I = 1 A = FEM 激發電流(per fit-current-matches-sim rule)。
%  Sign    : 物理 signed B·n+，all-source(翻下極激發 P1/P3/P6)(per charge-model-source-convention)。
% =========================================================================

clear; clc;

%% ---- config ----------------------------------------------------------------
R_select   = 150e-6;           % 取點半徑 [m]：(a) 決定載入哪個 R 的 ℓ̂、(b) 擬合/殘差只用此球內 air 節點
I_actual   = 1;                % 驅動電流 [A] = FEM 激發電流（1 A）
S_hall     = 130;              % Hall 靈敏度 [V/T]（EQ-730L：130 V/T）
N_I        = 6;                % FEM 模擬次數 = 6 個單線圈解
VARIANT    = 'sensor_spheres';       % 讀哪個 FEM 變體子夾：'standard'(baseline) | 'sensor_spheres'(sensor 加密)
ELL_LO     = 0.2e-3;           % ℓ̂ 一維搜尋下界 [m]（fit ℓ̂ 用）
ELL_HI     = 3.0e-3;           % ℓ̂ 一維搜尋上界 [m]

%% ---- paths -----------------------------------------------------------------
TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % mt_constants/import_ansys_data/filter_iron_nodes
addpath(fullfile(TREE,'code','function'));                                       % 本包模型輔助函式
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';  % FEM 結果根目錄
data_dir     = fullfile(TREE,'data');                                             % 規則#2：.mat 放本組 data/
if strcmp(VARIANT,'standard'); dsuf=''; else; dsuf=['_' VARIANT]; end              % baseline 保原檔名；變體加後綴
mat_out      = fullfile(data_dir, sprintf('calib_sensor_d%s.mat',dsuf));           % 本流程的解輸出（本組 data/）
if ~exist(fileparts(mat_out),'dir'); mkdir(fileparts(mat_out)); end

%  注意：ℓ̂ 不再從 fix_dir fit_KI_fixl 載入；改成本檔步驟3 自己 fit（需先有 P,B,Vmat）。

%% ---- 常數 + 電荷位置（極尖單位方向）---------------------------------------
cnst = mt_constants();                                 % 幾何常數（R_norm、極尖、SPH_OFST、k_m、mu_0 ...）
apdl_to_paper_idx = [1,3,6,5,2,4];                     % APDL coil j → paper pole 索引（判下極用）
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];  % 3×6 極尖位置（WP 框）
dhat = tip ./ vecnorm(tip);                            % 3×6 極尖單位方向
Pc   = dhat;                                           % 電荷(極)位置：放單位距離（build_S 用）

%% ---- 載 6-coil FEM 場、選 R 球內 air 節點、組殘差用的場 B -------------------
% 取 coil1 air 節點當工作點 p_i，限 R_select 球內（6 顆 coil 共用同一 mesh/順序）
d1   = import_ansys_data(fullfile(results_root, 'coil1', VARIANT),'wp','coil1');  % 載 'wp' 近場（variant 子夾）
air1 = filter_iron_nodes(d1.x,d1.y,d1.z,cnst,struct('visualize',false));             % 只留 air 節點（排除鐵內）
zwp1 = d1.z - cnst.SPH_OFST;                                                          % ANSYS z → WP 框 z
P_all = [d1.x(air1), d1.y(air1), zwp1(air1)];                                         % air 節點座標
insel = sum(P_all.^2,2) <= R_select^2;                                                % 取 R≤R_select 球內
P     = P_all(insel,:);  Np = size(P,1);                                              % 擬合區工作點
fprintf('選出 N_p = %d 個 air 節點（R ≤ %g µm 球內）\n', Np, R_select*1e6);

% 記錄 b_ij = −B^FEM（all-source 慣例的負號版；存 B(Np,3,N_I)）
B = zeros(Np,3,N_I);
for j = 1:N_I                                          % 逐線圈 j 載入
    if j == 1, dj = d1; airj = air1;                  % coil1 已載，重用
    else
        cn = sprintf('coil%d', j);
        dj = import_ansys_data(fullfile(results_root, cn, VARIANT),'wp',cn);
        airj = filter_iron_nodes(dj.x,dj.y,dj.z,cnst,struct('visualize',false));
    end
    Bj_all = -[dj.bx(airj), dj.by(airj), dj.bz(airj)];% 負號版場（與 PAGE1 慣例一致）
    B(:,:,j) = Bj_all(insel,:);                       % 只取 R 球內節點
end

%% ---- 步驟2：建 sensor 幾何 → 真實節點抽 sensor 電壓 Vmat（all-source）------
[sensor_pos, sensor_n] = build_sensor_geometry(cnst);                              % 6 顆 sensor 中心 + 法線 n+
[Vmat, exc_sign]       = extract_Vmat(results_root, cnst, apdl_to_paper_idx, ...    % 沿 n 圓柱選真實節點、V=S·B（含翻號）
                                      sensor_pos, sensor_n, S_hall, VARIANT);

%% ---- 步驟3：自己 fit ℓ̂（minimize profiled cost；不拿 fix_dir 的值）---------
%  J(ℓ̂)=Σ_j[‖b_j‖²−c_jᵀM_j⁻¹c_j]（消去 d、per-j profile；見 sensor_cost_lhat / cost_function.pdf）
costfun       = @(l) sensor_cost_lhat(l, P, B, Vmat, exc_sign, Pc, N_I);
opt           = optimset('TolX',1e-9,'Display','off');
[ell_hat, Jl] = fminbnd(costfun, ELL_LO, ELL_HI, opt);                              % 一維搜尋 ℓ̂
fprintf('步驟3：fit ℓ̂ = %.4f mm（min profiled cost = %.6e；搜尋 [%.2f, %.2f] mm）\n', ...
        ell_hat*1e3, Jl, ELL_LO*1e3, ELL_HI*1e3);
if abs(ell_hat-ELL_LO)<1e-9 || abs(ell_hat-ELL_HI)<1e-9
    warning('ℓ̂ 落在搜尋邊界 → 放寬 ELL_LO/ELL_HI 再確認（可能非內部極小）。');
end

%% ---- 步驟4：閉式解 shared d（6×1，無 g_H，照 V2 p2）+ cost J = Σ‖ε‖² --------
% solve_d 反矩陣在 Σ_j 外（對 i,j 加總）= shared d：d = (Σ_j V_j (Σ_i S_iᵀS_i) V_j)⁻¹ (Σ_j V_j (Σ_i S_iᵀ b_ij))
d = solve_d(P, B, Pc, ell_hat, Vmat, exc_sign, N_I);                       % 每極常數 d（no-gain，模型 b = S·V·d）
J = sensor_residual(P, B, Vmat, exc_sign, ell_hat, Pc, d, N_I);            % 殘差成本 cost J = Σ‖ε‖² [T²]

%% ---- 印結果 ----------------------------------------------------------------
fprintf('\n=============== Hall-sensor 模型結果（ℓ̂ = %.3f mm）===============\n', ell_hat*1e3);
fprintf('  sensor 電壓 Vmat [V]（列=sensor 極 P1..P6，欄=激發 coil1..6）：\n');
for i=1:6, fprintf('   % .3e % .3e % .3e % .3e % .3e % .3e\n', Vmat(i,:)); end
fprintf('  d（6×1，P1..P6，no-gain）：\n'); fprintf('   % .4e\n', d);
fprintf('  cost J = %.6e  [T^2]\n', J);
fprintf('===================================================================\n');

%% ---- 存解（交付物）--------------------------------------------------------
save(mat_out, 'd','Vmat','exc_sign','ell_hat','J', ...
              'S_hall','sensor_pos','sensor_n','apdl_to_paper_idx','R_select');
fprintf('已存 %s\n', mat_out);
