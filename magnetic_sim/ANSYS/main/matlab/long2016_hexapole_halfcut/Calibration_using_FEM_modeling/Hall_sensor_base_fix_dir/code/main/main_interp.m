%% main_interp.m -- Hall-sensor 模型求 d（standard 粗網格「內插版」Vmat）
% =========================================================================
%  與 main.m 同流程，差別只在 sensor 電壓 Vmat 的來源：
%    main.m         ：VARIANT 變體的「真實節點圓柱平均」(extract_Vmat)。
%    main_interp.m  ：standard 粗網格「真·FEM tet 重心內插」(extract_Vmat_interp)
%                     —— 在與加密版同數量/同位置的取樣點上內插 standard 場。
%  目的：看「未加密 + 內插」解出的 d / cost J 與加密版差多少。
%
%  重要：per-j cost 的 V_j 會抵銷（見 sensor_cost_lhat）→ ℓ̂ 與 Vmat 無關，
%        故 ℓ̂ 仍 = 0.856mm；只有 d 與 cost J 受 interp Vmat 影響。
%
%  Current : I = 1 A = FEM 激發電流。Sign : all-source（翻下極 P1/P3/P6）。
% =========================================================================

clear; clc;

%% ---- config ----
R_select   = 150e-6;           % 取點半徑 [m]（擬合/殘差用此球內 air 節點）
S_hall     = 130;              % Hall 靈敏度 [V/T]
N_I        = 6;                % 6 個單線圈解
VARIANT    = 'standard';       % 內插版的「場」一律用 standard 粗網格
N_UNIFORM  = 1000;             % 每 sensor 圓柱內均勻取樣點數（內插）
ELL_LO     = 0.2e-3;           % ℓ̂ 搜尋下界 [m]
ELL_HI     = 3.0e-3;           % ℓ̂ 搜尋上界 [m]

%% ---- paths ----
TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
addpath(fullfile(TREE,'code','function'));
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
MESH_CSV_DIR = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data\mesh\standard\csv';
data_dir     = fullfile(TREE,'data');                          % 規則#2：.mat 放本組 data/
mat_out      = fullfile(data_dir, 'calib_sensor_d_interp.mat');
if ~exist(fileparts(mat_out),'dir'); mkdir(fileparts(mat_out)); end

%% ---- 常數 + 電荷位置 ----
cnst = mt_constants();
apdl_to_paper_idx = [1,3,6,5,2,4];
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);
Pc   = dhat;

%% ---- 載 6-coil 場（standard），選 R 球內 air 節點 → P, B ----
d1   = import_ansys_data(fullfile(results_root, 'coil1', VARIANT),'wp','coil1');
air1 = filter_iron_nodes(d1.x,d1.y,d1.z,cnst,struct('visualize',false));
zwp1 = d1.z - cnst.SPH_OFST;
P_all = [d1.x(air1), d1.y(air1), zwp1(air1)];
insel = sum(P_all.^2,2) <= R_select^2;
P     = P_all(insel,:);  Np = size(P,1);
fprintf('選出 N_p = %d 個 air 節點（R ≤ %g µm 球內，standard）\n', Np, R_select*1e6);

B = zeros(Np,3,N_I);
for j = 1:N_I
    if j == 1, dj = d1; airj = air1;
    else
        cn = sprintf('coil%d', j);
        dj = import_ansys_data(fullfile(results_root, cn, VARIANT),'wp',cn);
        airj = filter_iron_nodes(dj.x,dj.y,dj.z,cnst,struct('visualize',false));
    end
    Bj_all = -[dj.bx(airj), dj.by(airj), dj.bz(airj)];
    B(:,:,j) = Bj_all(insel,:);
end

%% ---- 步驟2：sensor 電壓 Vmat（standard 內插版）----
[sensor_pos, sensor_n] = build_sensor_geometry(cnst);
[Vmat, exc_sign]       = extract_Vmat_interp(results_root, cnst, apdl_to_paper_idx, ...
                                             sensor_pos, sensor_n, S_hall, MESH_CSV_DIR, N_UNIFORM);

%% ---- 步驟3：fit ℓ̂（per-j profiled cost；V_j 會抵銷 → 同電荷擬合）----
costfun       = @(l) sensor_cost_lhat(l, P, B, Vmat, exc_sign, Pc, N_I);
opt           = optimset('TolX',1e-9,'Display','off');
[ell_hat, Jl] = fminbnd(costfun, ELL_LO, ELL_HI, opt);
fprintf('步驟3：fit ℓ̂ = %.4f mm（min profiled cost = %.6e）\n', ell_hat*1e3, Jl);

%% ---- 步驟4：解 shared d + cost J ----
d = solve_d(P, B, Pc, ell_hat, Vmat, exc_sign, N_I);
J = sensor_residual(P, B, Vmat, exc_sign, ell_hat, Pc, d, N_I);

%% ---- 印結果 ----
fprintf('\n=========== Hall-sensor 內插版結果（standard 粗網格，ℓ̂ = %.3f mm）===========\n', ell_hat*1e3);
fprintf('  sensor 電壓 Vmat [V]（列=sensor P1..P6，欄=激發 coil1..6）：\n');
for i=1:6, fprintf('   % .3e % .3e % .3e % .3e % .3e % .3e\n', Vmat(i,:)); end
fprintf('  d（6×1，P1..P6，no-gain）：\n'); fprintf('   % .4e\n', d);
fprintf('  cost J = %.6e  [T^2]\n', J);
fprintf('================================================================================\n');

%% ---- 存解 ----
save(mat_out, 'd','Vmat','exc_sign','ell_hat','J','S_hall','sensor_pos','sensor_n', ...
              'apdl_to_paper_idx','R_select','VARIANT','N_UNIFORM');
fprintf('已存 %s\n', mat_out);
