%% solve_DH_full.m -- Hall-sensor 滿 6×6 解耦矩陣 D_H（對照 main.m 的對角 d）
% =========================================================================
%  目的：解「滿矩陣」解耦校正 D_H（6×6），把交叉耦合的 sensor 電壓 V 反解成每極電荷：
%        q_j = D_H · Vmat(:,j)      （對照 main.m 的對角版 q_j = diag(d)·Vmat(:,j)）
%  模型場：b_ij = S_i · D_H · Vmat(:,j)（無 g_H、all-source；S_i = build_S 在 ℓ̂ 下）。
%
%  ── 閉式解（min_{D_H} Σ_ij‖S_i·D_H·Vmat(:,j) − b_ij‖²）────────────────────
%        D_H = (M \ C) / Vmat
%     其中 M = Σ_i S_iᵀS_i（6×6）、C = [c_1..c_6]、c_j = Σ_i S_iᵀ b_ij（all-source）。
%     `M\C` 即「每次激發的自由電荷最佳解」（= fix_dir fit_KI_fixl 在找的）→
%     q_j = D_H·Vmat(:,j) 還原自由電荷 → cost_J 應 ≈ fix_dir 的自由電荷下界(~0.004)。
%     物理：D_H 非對角 = cross-talk 反解權重（= 耦合矩陣 G 的逆 G⁻¹ 的非對角）。
%
%  ── 與 main.m 的關係 ──────────────────────────────────────────────────────
%  獨立對照腳本，不動 main（對角 d、LAB406 per-pole 模型）。重用同一組 function/ 與
%  同一個 ℓ̂（R=150，fix_dir fit_KI_fixl=0.856）。輸出 calib_DH_full.mat 到 fitting_d/。
%
%  Current : I = 1 A = FEM 激發電流。Sign：物理 signed B·n+、all-source（翻下極 P1/P3/P6）。
% =========================================================================

clear; clc;

%% ---- config ----------------------------------------------------------------
R_select = 150e-6;             % 取點/載 ℓ̂ 半徑 [m]
S_hall   = 130;                % Hall 靈敏度 [mV/mT]（EQ-730L；配 B[mT]→V[mV]）
N_I      = 6;                  % 6 個單線圈解

%% ---- paths -----------------------------------------------------------------
TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % mt_constants/import_ansys_data/filter_iron_nodes
addpath(fullfile(TREE,'code','function'));                                       % build_S / build_sensor_geometry / extract_Vmat（重用）
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
data_dir     = fullfile(TREE,'data');                         % 規則#2：本組 .mat
fixdir_data  = fullfile(fileparts(TREE),'fix_dir','data');    % fit_fixl ℓ̂ 來源（fix_dir/data）
mat_out      = fullfile(data_dir, 'calib_DH_full.mat');                          % 輸出（與 calib_sensor_d.mat 並列於 data/）
if ~exist(fileparts(mat_out),'dir'); mkdir(fileparts(mat_out)); end

%% ---- 載 ℓ̂（同 main：fix_dir fit_KI_fixl 在 R=150 的解）---------------------
RRR  = round(R_select*1e6);
ellf = fullfile(fixdir_data, sprintf('fit_fixl_R%03dum.mat', RRR));
assert(exist(ellf,'file')==2, 'ℓ̂ 來源不存在：%s（請先跑 fix_dir main）', ellf);
SL = load(ellf, 'ell','J');  ell_hat = SL.ell;  J_fixdir = SL.J;                 % ℓ̂ [µm]、fix_dir 自由電荷 cost（對照，mT²）
fprintf('載入 ℓ̂ = %.2f µm（fix_dir fit_KI_fixl）；fix_dir cost J = %.6e\n', ell_hat, J_fixdir);

%% ---- 常數 + 電荷位置 -------------------------------------------------------
cnst = mt_constants();
apdl_to_paper_idx = [1,3,6,5,2,4];
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);  Pc = dhat;

%% ---- 載 6-coil FEM、取 R 球內 air 節點、組負號版場 B ------------------------
d1   = import_ansys_data(fullfile(results_root, 'coil1', 'standard'),'wp','coil1');
air1 = filter_iron_nodes(d1.x,d1.y,d1.z,cnst,struct('visualize',false));
zwp1 = d1.z - cnst.SPH_OFST;
P_all = [d1.x(air1), d1.y(air1), zwp1(air1)];
insel = sum(P_all.^2,2) <= R_select^2;                % 選球用公尺（配 R_select[m]）
P     = P_all(insel,:) * 1e6;  Np = size(P,1);        % m → µm（配 ell_hat[µm]、build_S 比值）
B = zeros(Np,3,N_I);
for j = 1:N_I
    if j == 1, dj = d1; airj = air1;
    else
        cn = sprintf('coil%d', j);
        dj = import_ansys_data(fullfile(results_root, cn, 'standard'),'wp',cn);
        airj = filter_iron_nodes(dj.x,dj.y,dj.z,cnst,struct('visualize',false));
    end
    Bj = -1e3*[dj.bx(airj), dj.by(airj), dj.bz(airj)];   % 負號版場；ANSYS Tesla → ×1e3 原生 mT
    B(:,:,j) = Bj(insel,:);
end
fprintf('R≤%g µm 球內 air 節點 Np = %d\n', R_select*1e6, Np);

%% ---- sensor 幾何 → 真實節點抽 Vmat（all-source）---------------------------
[sensor_pos, sensor_n] = build_sensor_geometry(cnst);
[Vmat, exc_sign]       = extract_Vmat(results_root, cnst, apdl_to_paper_idx, sensor_pos, sensor_n, S_hall);

%% ---- 組 M、C（all-source 目標）→ 解滿矩陣 D_H = (M\C)/Vmat ------------------
M = zeros(6,6);  C = zeros(6,N_I);
for i = 1:Np
    Si = build_S(P(i,:), ell_hat, Pc);                 % 3×6 空間核
    M  = M + Si.'*Si;                                  % Σ S_iᵀS_i
    for j = 1:N_I
        bij    = exc_sign(j) * (-squeeze(B(i,:,j)).'); % all-source 物理場 b_ij
        C(:,j) = C(:,j) + Si.'*bij;                    % Σ S_iᵀ b_ij
    end
end
Q_free = M \ C;                                        % 每次激發的自由電荷最佳解（= fix_dir 找的）
D_H    = Q_free / Vmat;                                % 滿 6×6 解耦矩陣：q_j = D_H·Vmat(:,j)

%% ---- 直接評估 cost_J = Σ‖S·D_H·Vmat(:,j) − b_ij‖² --------------------------
cost_J = 0;
for j = 1:N_I
    qj = D_H * Vmat(:,j);                              % 該激發的電荷向量（= 還原的自由電荷）
    for i = 1:Np
        Si = build_S(P(i,:), ell_hat, Pc);
        bm = Si * qj;
        bf = exc_sign(j) * (-squeeze(B(i,:,j)).');
        cost_J = cost_J + sum((bm-bf).^2);
    end
end

%% ---- 印結果 + D_H 的非對角分析 --------------------------------------------
offmask = ~eye(6);
diagn   = norm(diag(D_H));
offn    = norm(D_H(offmask));
fprintf('\n=============== 滿 D_H 解耦結果（ℓ̂=%.2f µm, R=%d µm）===============\n', ell_hat, RRR);
fprintf('  D_H (6×6, 列=電荷極 P1..P6, 欄=sensor P1..P6；單位 mT/mV)：\n');
for i=1:6, fprintf('   % .3e % .3e % .3e % .3e % .3e % .3e\n', D_H(i,:)); end
fprintf('  ‖diag(D_H)‖ = %.3e ; ‖off-diag(D_H)‖ = %.3e ; off/diag = %.2f ; max|off| = %.3e\n', ...
        diagn, offn, offn/diagn, max(abs(D_H(offmask))));
fprintf('  --------------------------------------------------------------\n');
fprintf('  cost_J (滿 D_H 解耦)      = %.6e  [mT^2]\n', cost_J);
fprintf('  cost_J (fix_dir 自由電荷) = %.6e  [mT^2]（對照，= 下界）\n', J_fixdir);
fprintf('==================================================================\n');

%% ---- 存解 ------------------------------------------------------------------
save(mat_out, 'D_H','cost_J','ell_hat','Vmat','exc_sign','S_hall','apdl_to_paper_idx','R_select');
fprintf('已存 %s\n', mat_out);
