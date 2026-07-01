%% solve_DH_interp.m -- 滿 6×6 解耦 D_H（內插版 Vmat，1000 均勻取樣點）
% =========================================================================
%  同 solve_DH_full.m 的滿解耦模型，但 Vmat 改用「standard 粗網格 tet 重心內插」
%  在每 sensor 圓柱內均勻撒 1000 點得到（extract_Vmat_interp）。
%
%  模型場：b_ij = S_i · D_H · Vmat(:,j)，q_j = D_H·Vmat(:,j)。
%  閉式：  D_H = (M\C)/Vmat，M=ΣS_iᵀS_i、C_j=ΣS_iᵀb_ij。
%  恆等：  D_H·Vmat = Q_free（= M\C）→ cost_J 恆 = fix_dir 自由電荷下界（與點數無關）；
%          但 D_H 矩陣值 = Q_free·Vmat⁻¹ 隨 Vmat(點數) 變。
%  V 表示：每次模擬 j 的 V_j = diag(Vmat(:,j))（6×6 對角）→ 6 個 → Vtensor 6×6×6。
%
%  Current : I = 1 A = FEM 激發電流。Sign：all-source（翻下極 P1/P3/P6）。
% =========================================================================

clear; clc;

%% ---- config ----
R_select  = 150e-6;            % 取點/載 ℓ̂ 半徑 [m]
S_hall    = 130;              % Hall 靈敏度 [mV/mT]（配 B[mT]→V[mV]）
N_I       = 6;                % 6 個單線圈解
N_UNIFORM = 1000;             % 每 sensor 圓柱內均勻取樣點數（內插）

%% ---- paths ----
TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
addpath(fullfile(TREE,'code','function'));
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
MESH_CSV_DIR = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data\mesh\standard\csv';
data_dir     = fullfile(TREE,'data');                         % 規則#2：本組 .mat
fixdir_data  = fullfile(fileparts(TREE),'fix_dir','data');    % fit_fixl ℓ̂ 來源（fix_dir/data）
mat_out      = fullfile(data_dir, 'calib_DH_interp.mat');
if ~exist(fileparts(mat_out),'dir'); mkdir(fileparts(mat_out)); end

%% ---- 載 ℓ̂（fix_dir fit_KI_fixl @ R=150）----
RRR  = round(R_select*1e6);
ellf = fullfile(fixdir_data, sprintf('fit_fixl_R%03dum.mat', RRR));
assert(exist(ellf,'file')==2, 'ℓ̂ 來源不存在：%s', ellf);
SL = load(ellf,'ell','J'); ell_hat = SL.ell; J_fixdir = SL.J;   % ℓ̂ [µm]、cost [mT²]
fprintf('載入 ℓ̂ = %.2f µm；fix_dir 自由電荷 cost J = %.6e\n', ell_hat, J_fixdir);

%% ---- 常數 + 電荷位置 ----
cnst = mt_constants();
apdl_to_paper_idx = [1,3,6,5,2,4];
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);  Pc = dhat;

%% ---- 載 6-coil WP 場（standard），取 R 球內 air 節點 → P, B（負號版）----
d1   = import_ansys_data(fullfile(results_root,'coil1','standard'),'wp','coil1');
air1 = filter_iron_nodes(d1.x,d1.y,d1.z,cnst,struct('visualize',false));
zwp1 = d1.z - cnst.SPH_OFST;
P_all = [d1.x(air1), d1.y(air1), zwp1(air1)];
insel = sum(P_all.^2,2) <= R_select^2;                % 選球用公尺
P     = P_all(insel,:) * 1e6;  Np = size(P,1);        % m → µm（配 ell_hat[µm]）
B = zeros(Np,3,N_I);
for j = 1:N_I
    if j==1, dj=d1; airj=air1;
    else
        cn=sprintf('coil%d',j); dj=import_ansys_data(fullfile(results_root,cn,'standard'),'wp',cn);
        airj=filter_iron_nodes(dj.x,dj.y,dj.z,cnst,struct('visualize',false));
    end
    Bj=-1e3*[dj.bx(airj),dj.by(airj),dj.bz(airj)]; B(:,:,j)=Bj(insel,:);   % Tesla → ×1e3 原生 mT
end
fprintf('R≤%g µm 球內 air 節點 Np = %d\n', R_select*1e6, Np);

%% ---- sensor 電壓 Vmat（內插版，1000 點）----
[sensor_pos, sensor_n] = build_sensor_geometry(cnst);
[Vmat, exc_sign]       = extract_Vmat_interp(results_root, cnst, apdl_to_paper_idx, ...
                                             sensor_pos, sensor_n, S_hall, MESH_CSV_DIR, N_UNIFORM);

%% ---- M, C → Q_free → 滿 D_H = Q_free/Vmat ----
M = zeros(6,6); C = zeros(6,N_I);
for i = 1:Np
    Si = build_S(P(i,:), ell_hat, Pc);
    M  = M + Si.'*Si;
    for j = 1:N_I
        bij = exc_sign(j)*(-squeeze(B(i,:,j)).');
        C(:,j) = C(:,j) + Si.'*bij;
    end
end
Q_free = M \ C;
D_H    = Q_free / Vmat;                          % 滿 6×6 解耦

%% ---- cost_J = Σ‖S·D_H·Vmat(:,j) − b_ij‖² ----
cost_J = 0;
for j = 1:N_I
    qj = D_H * Vmat(:,j);
    for i = 1:Np
        bm = build_S(P(i,:),ell_hat,Pc)*qj;
        bf = exc_sign(j)*(-squeeze(B(i,:,j)).');
        cost_J = cost_J + sum((bm-bf).^2);
    end
end

%% ---- V_j = diag(Vmat(:,j))（6×6 對角）→ Vtensor 6×6×6 ----
Vtensor = zeros(6,6,N_I);
for j = 1:N_I, Vtensor(:,:,j) = diag(Vmat(:,j)); end

%% ---- 印結果 ----
offmask = ~eye(6); diagn = norm(diag(D_H)); offn = norm(D_H(offmask));
fprintf('\n=========== 內插版 滿 D_H 解耦（1000 點，ℓ̂=%.2f µm）===========\n', ell_hat);
fprintf('  D_H (6×6, 列=電荷極 P1..P6, 欄=sensor P1..P6；單位 mT/mV)：\n');
for i=1:6, fprintf('   % .3e % .3e % .3e % .3e % .3e % .3e\n', D_H(i,:)); end
fprintf('  ‖diag‖=%.3e ‖off‖=%.3e off/diag=%.2f\n', diagn, offn, offn/diagn);
fprintf('  cost_J (內插滿 D_H)       = %.6e  [mT^2]\n', cost_J);
fprintf('  cost_J (fix_dir 自由電荷) = %.6e  [mT^2]（下界）\n', J_fixdir);
fprintf('==================================================================\n');

%% ---- 存解 ----
save(mat_out, 'D_H','Vtensor','Vmat','cost_J','ell_hat','exc_sign', ...
              'S_hall','apdl_to_paper_idx','R_select','N_UNIFORM');
fprintf('已存 %s\n', mat_out);
