%% svd_SiHI_WP.m -- SVD of the current->field transfer T = S_i·Ĥ_I at the WP center (ACTUATOR frame)
% =========================================================================
%  目的：在「工作點中心」p = 0 評估 fix-ℓ 電荷模型的電流→場轉移矩陣
%           T = S_i · Ĥ_I   （3×6，單位 mT/A），並做 SVD  T = U·Σ·Wᵀ，回報 U、Σ、Wᵀ。
%
%  座標系（統一 actuator）：
%     文件 kernel 用 p̄_ci = Pc_base（理想 actuator 格 [+u −u +v −v +w −w]），故本腳本
%     **直接在 actuator 框建 S_i**（用 Pc_base；節點旋進 actuator：p_act = R_act·p）。
%     R_act = [d̂_P1, d̂_P3, d̂_P5]ᵀ（measure→actuator；R_act·dhat = Pc_base）。
%     ★ fix-ℓ 校正的 K̄/ℓ̂/^Bĝ_I（故 Ĥ_I）是極空間+純量、框無關 → Ĥ_I 不變；
%       只有「空間輸出」S_i/T/U 帶框。等價：T_actuator = R_act · T_measure。
%
%  由來：b_i = S_i·^Bĝ_I·K̄·F = S_i·Ĥ_I·F（Ĥ_I = ^Bĝ_I·K̄ = gB·Khat）。S_i·Ĥ_I 即作用在電流上的轉移。
%
%  WP 中心的解析簡化：
%     p = 0 ⇒ p̄ = R_act·p/ℓ̂ = 0 ⇒ S_i(:,k) = (0 − Pc_base_k)/‖0 − Pc_base_k‖³ = −Pc_base_k。
%     ∴ S_i(WP) = −Pc_base（3×6），且與 ℓ̂ 無關（電荷距中心恆 = ℓ̂，ℓ̂ 在 kernel 對消）。
%
%  SVD 幾何讀法（兩指標）：T 把單位電流球映成 3D 場橢球，半軸 = σ1..σ3、方向 = U(:,1..3)。
%     形狀大小（增益）：σ1（最大）、(σ1σ2σ3)^{1/3}（幾何平均）；越大 = 致動越強。
%     均勻度（球度）  ：iso = σ3/σ1 ∈(0,1]（=1 等向）、κ = σ1/σ3；越接近 1 = 越接近正球。
%     Σ、Wᵀ 框無關；U 為 actuator 框（= R_act·U_measure）。
%
%  輸入：fix_dir/data/fit_fixl_R150um_gap200um_mueq.mat（gB=^Bĝ_I、Khat=K̄；ell 中心點不用）。
%  輸出：console（actuator 框 U/Σ/Wᵀ/T/σ/κ/指標）+ fix_dir/data/svd_SiHI_WP_gap200um_mueq.mat。
% =========================================================================

clear; clc;

%% ---- config ----------------------------------------------------------------
VARIANT = 'gap200um_mueq';       % 校正變體（gap200 2 段式 μ_eff）
p_wp    = [0; 0; 0];             % 工作點 = WP 中心 [m]

%% ---- paths -----------------------------------------------------------------
TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\fix_dir'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % mt_constants
data_dir = fullfile(TREE,'data');

%% ---- 載入校正參數：Ĥ_I = ^Bĝ_I·K̄ = gB·Khat（框無關）----------------------
S = load(fullfile(data_dir, sprintf('fit_fixl_R150um_%s.mat', VARIANT)), 'ell','gB','Khat');
ell_m  = S.ell * 1e-6;           % µm → m（中心點不影響結果，僅供記錄）
Hhat_I = S.gB * S.Khat;          % Ĥ_I（6×6，mT/A）

%% ---- actuator 框：Pc_base（文件理想格）+ R_act（measure→actuator）----------
cnst = mt_constants();
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];   % 3×6 [m]（measure）
dhat = tip ./ vecnorm(tip);                                      % 3×6 measure 單位方向
Pc_base = [ 1 -1  0  0  0  0;                                    % actuator 理想格 [+u -u +v -v +w -w]
            0  0  1 -1  0  0;
            0  0  0  0  1 -1];
R_act = [dhat(:,1), dhat(:,3), dhat(:,5)].';                     % measure→actuator（û=P1, v̂=P3, ŵ=P5）
assert(max(abs(R_act*dhat - Pc_base), [], 'all') < 1e-9, 'R_act*dhat 必須 = Pc_base');

%% ---- 組 S_i（actuator 框；WP 中心自然得 −Pc_base、ℓ̂-free）------------------
pbar = (R_act * p_wp) / ell_m;   % actuator 框無因次工作點（此處 = 0）
Dk   = pbar - Pc_base;           % 3×6：p̄ − Pc_base_k
Si   = Dk ./ (vecnorm(Dk).^3);   % 3×6 庫倫核（p=0 → Si = −Pc_base）
chk_center = max(abs(Si + Pc_base), [], 'all');   % 驗中心解析：應 ≈ 0

%% ---- T = S_i·Ĥ_I（actuator），SVD ------------------------------------------
T = Si * Hhat_I;                 % 3×6，mT/A（actuator 框：列 = b_u/b_v/b_w）
T_meas = (-dhat) * Hhat_I;       % measure 框 WP 中心 T（供交叉驗證）
chk_rot = max(abs(T - R_act*T_meas), [], 'all');      % 應 ≈ 0（T_act = R_act·T_meas）
[U, Sig, W] = svd(T);            % T = U·Σ·Wᵀ
WT   = W.';                      % Wᵀ（6×6）
sv   = diag(Sig);                % 3 個奇異值 [mT/A]
recon = norm(U*Sig*W.' - T, 'fro') / norm(T,'fro');   % 重建檢查

%% ---- 兩指標：形狀大小（增益）+ 均勻度（球度）------------------------------
gain_max = sv(1);                        % 最大增益（橢球最長半軸）
gain_gm  = (sv(1)*sv(2)*sv(3))^(1/3);    % 幾何平均增益（橢球「體積」半徑）
iso      = sv(3) / sv(1);                % 均勻度 ∈(0,1]（=1 完美等向/正球）
kappa    = sv(1) / sv(3);                % 條件數 σ1/σ3（=1 完美等向）

%% ---- 參考：P1 激發在 WP 中心的場（= T 第 1 欄，actuator 框）----------------
F    = [1;0;0;0;0;0];            % P1 激發電流向量（paper 序，1 A）
b_P1 = T * F;                    % 3×1 [mT]（= T(:,1)）

%% ---- console 回報（actuator 框）-------------------------------------------
fprintf('\n===== SVD of T = S_i·Ĥ_I at WP center (fix-ℓ, ACTUATOR frame, variant=%s) =====\n', VARIANT);
fprintf('工作點 p_wp = [% .0f % .0f % .0f]；ℓ̂ = %.2f µm（中心點不影響 S_i）\n', p_wp, S.ell);
fprintf('框：actuator（Pc_base；列 = b_u/b_v/b_w = P1/P3/P5 軸）。Σ、Wᵀ 框無關；U = R_act·U_meas。\n');
fprintf('驗：max|S_i+Pc_base| = %.2e、max|T − R_act·T_meas| = %.2e、‖UΣWᵀ−T‖/‖T‖ = %.2e（皆應≈0）\n\n', ...
        chk_center, chk_rot, recon);

fprintf('T = S_i·Ĥ_I  [mT/A]（列 = b_u/b_v/b_w；欄 = 電流 P1..P6）：\n');
for r = 1:3, fprintf('  % .4e % .4e % .4e % .4e % .4e % .4e\n', T(r,:)); end

fprintf('\nU（3×3，磁場端主軸；actuator 框；列 = b_u/b_v/b_w，欄 = 主方向 1..3）：\n');
for r = 1:3, fprintf('  % .6f % .6f % .6f\n', U(r,:)); end

fprintf('\nΣ（3×6，奇異值 [mT/A]）：\n');
for r = 1:3, fprintf('  % .4e % .4e % .4e % .4e % .4e % .4e\n', Sig(r,:)); end

fprintf('\nWᵀ（6×6，電流端主軸；列 = 主組合 1..6，欄 = 電流 P1..P6）：\n');
for r = 1:6, fprintf('  % .6f % .6f % .6f % .6f % .6f % .6f\n', WT(r,:)); end

fprintf('\nσ = [% .4e % .4e % .4e] mT/A\n', sv);
fprintf('指標｜形狀大小(增益)：σ1=%.4f、幾何平均=(σ1σ2σ3)^{1/3}=%.4f mT/A\n', gain_max, gain_gm);
fprintf('指標｜均勻度(球度)  ：iso=σ3/σ1=%.4f（→1 越圓）、κ=σ1/σ3=%.4f\n', iso, kappa);
fprintf('參考：P1 激發 (F=e1) 在 WP 中心的場 b_P1 = [% .4e % .4e % .4e] mT（actuator）\n', b_P1);
fprintf('=====================================================================\n');

%% ---- 存 .mat ---------------------------------------------------------------
out_mat = fullfile(data_dir, sprintf('svd_SiHI_WP_%s.mat', VARIANT));
save(out_mat, 'T','U','Sig','W','WT','sv','kappa','iso','gain_max','gain_gm', ...
              'recon','chk_center','chk_rot','p_wp','Pc_base','R_act','dhat', ...
              'Si','T_meas','Hhat_I','ell_m','F','b_P1','VARIANT');
fprintf('已存 %s\n', out_mat);
