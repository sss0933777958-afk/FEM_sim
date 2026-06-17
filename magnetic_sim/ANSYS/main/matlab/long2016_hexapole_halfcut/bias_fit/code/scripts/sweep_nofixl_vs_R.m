%% sweep_nofixl_vs_R.m
%  ============================================================================
%  no_fix_l(18-param bias)模型的「取樣半徑 R 趨勢掃描」。
%  對每個 R(50:5:500 um)做一次 bias fit(actuator 框,ell_hat + e_hat joint),
%  存下 ell_hat(R)、gB_hat(R)、K̄_I(R)、e_hat(R),供趨勢圖與每半徑 LaTeX 使用。
%
%  與既有 sweep_alln_vs_R.m(舊全-K̂ 模型)的差異:
%    * 模型 = no_fix_l bias(文件 no_fix_l.pdf):actuator 框,Pc = Pc_base + E(e),
%      g_j profile 掉後 gauge 出 gB·K̄_I(與 calib_fem_bias.m 同骨架)。
%    * 每半徑用「全部 in-ball 節點」(無 cap、確定性),趨勢平滑且更超定(不過擬合)。
%    * 不算 NRMSE(驗證才用)。
%    * fit 用 lsqnonlin(profiled 殘差向量),避開 fminunc 在 R2025b 的 optim:fminusub
%      訊息目錄 bug(lsqnonlin 在既有 sweep 已驗證無此問題)。
%  ============================================================================
clear; clc;

addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');         % ansys_path / matlab_path
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\hexapole-long2016\analysis');   % mt_constants / import_ansys_data / filter_iron_nodes
model = 'long2016_hexapole_halfcut';

% ---- 頂層參數 --------------------------------------------------------------
cnst      = mt_constants();
R_um      = 50:5:500;                       % 取樣半徑掃描格 [um](91 點)
Nmax      = Inf;                            % 全取樣(無 cap;用全部 in-ball 節點)
ell0      = 0.5e-3;                          % bias fit 的 ell 初始猜 [m]
I_actual  = 1;                              % 模型電流 = FEM 激發(1A;per fit-current-matches-sim)
apdl_to_paper_idx = [1,3,6,5,2,4];          % APDL coil j → paper pole 索引
coil_sign = [1 -1 1 -1 -1 1];               % all-source 顯示翻上極(P2/P4/P5),存供 LaTeX 用
dataset   = 'all';                          % standard mesh 的 'all' dataset

% ---- actuator 旋轉 R(measure→actuator;原點在 WP 中心)----------------------
tip   = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];  % 3x6 極尖位置(measure 框)
dhat  = tip ./ vecnorm(tip);                                    % 3x6 單位方向
R_act = [dhat(:,1), dhat(:,3), dhat(:,5)].';                    % û=P1, v̂=P3, ŵ=P5(per dissertation +x_a/+y_a/+z_a)
Pc_base = [ 1 -1  0  0  0  0;                                   % 文件理想格點 [+û -û +v̂ -v̂ +ŵ -ŵ]
            0  0  1 -1  0  0;
            0  0  0  0  1 -1];
assert(abs(det(R_act)-1) < 1e-9, 'R_act 必須是 proper rotation');
assert(max(abs(R_act*dhat - Pc_base), [], 'all') < 1e-9, 'R_act*dhat 必須等於 Pc_base');

% ---- 電流矩陣 F(permutation, rank 6)---------------------------------------
N_I = 6;
F   = zeros(6, N_I);
for j = 1:N_I, F(apdl_to_paper_idx(j), j) = 1; end

% ---- 載 6 顆 coil 一次,旋進 actuator 框 ------------------------------------
fprintf('loading 6 coils (standard ''%s'') ...\n', dataset);
d1   = import_ansys_data(ansys_path(model,'coil1','standard'), dataset, 'coil1');
air1 = filter_iron_nodes(d1.x,d1.y,d1.z,cnst,struct('visualize',false));
zwp1 = d1.z - cnst.SPH_OFST;
P_meas = [d1.x(air1), d1.y(air1), zwp1(air1)];                 % 共享 air 點集(measure 框)
Nair   = size(P_meas,1);
Pa     = (R_act * P_meas.').';                                 % 旋進 actuator 框(一次)
r2     = sum(P_meas.^2, 2);                                    % 各點半徑^2(範數旋轉不變)

Ba = zeros(Nair, 3, N_I);                                      % 各 coil 的 all-source B(actuator 框)
for k = 1:N_I
    if k == 1, dk = d1; airk = air1;
    else
        cn = sprintf('coil%d', k);
        dk = import_ansys_data(ansys_path(model,cn,'standard'), dataset, cn);
        airk = filter_iron_nodes(dk.x,dk.y,dk.z,cnst,struct('visualize',false));
    end
    Bk = -[dk.bx(airk), dk.by(airk), dk.bz(airk)];            % all-source(取負)
    Ba(:,:,k) = (R_act * Bk.').';                             % 旋進 actuator 框
end
fprintf('  loaded; Nair = %d\n', Nair);

% ---- lsqnonlin 設定 --------------------------------------------------------
opts = optimoptions('lsqnonlin','Display','off', ...
                    'MaxFunctionEvaluations',1e5,'MaxIterations',4e3, ...
                    'FunctionTolerance',1e-20,'StepTolerance',1e-12);

% ---- 掃描 ------------------------------------------------------------------
nR     = numel(R_um);
ell_R  = nan(1,nR);  gB_R = nan(1,nR);  npts = nan(1,nR);
Ksave  = zeros(6,6,nR);   Esave = zeros(18,nR);             % K̄_I(R) 與 [e_hat;e6z](R)
fprintf('\n  R[um]  npts/coil   ell[mm]  gB[1e-3]\n');
for ri = 1:nR
    R = R_um(ri)*1e-6;
    idx = find(r2 < R^2);                                   % in-ball 節點(共享,全取樣)
    if numel(idx) > Nmax, idx = idx(1:Nmax); end            % Nmax=Inf → 不裁(全用)
    npts(ri) = numel(idx);

    P = Pa(idx, :);                                         % actuator 框近場點(Np x 3)
    Bstack = zeros(3*numel(idx), N_I);                      % 堆疊 b_ij(3Np x N_I)
    for j = 1:N_I, Bstack(:,j) = reshape(Ba(idx,:,j).', [], 1); end

    x0 = [ell0*1e3; zeros(17,1)];                          % 初值:ell(mm)、e=0
    xf = lsqnonlin(@(x) bias_resid(x, P, Bstack, Pc_base), x0, [], [], opts);
    ell = xf(1)*1e-3;  e_hat = xf(2:18);
    e6z = e_hat(1) - e_hat(4) + e_hat(8) - e_hat(11) + e_hat(15);
    Pc  = make_Pc(e_hat, Pc_base);
    [KbarI, gB] = gauge_KI(ell, Pc, P, Bstack, F);

    ell_R(ri) = ell;  gB_R(ri) = gB;
    Ksave(:,:,ri) = KbarI;  Esave(:,ri) = [e_hat; e6z];
    fprintf('  %4d   %7d     %.4f   %6.3f\n', R_um(ri), npts(ri), ell*1e3, gB*1e3);
end

% ---- 存檔(不含 NRMSE)------------------------------------------------------
data_dir = matlab_path(model, 'charge_fit', 'fitting_trend');
save(fullfile(data_dir,'sweep_nofixl_vs_R.mat'), ...
     'R_um','ell_R','gB_R','Ksave','Esave','npts','Nmax','I_actual', ...
     'apdl_to_paper_idx','coil_sign','R_act','Pc_base');
fprintf('\nwrote %s\n', fullfile(data_dir,'sweep_nofixl_vs_R.mat'));

%% ===========================================================================
%  本地函式(模型骨架同 calib_fem_bias.m)
%  ===========================================================================
function r = bias_resid(x, P, Bstack, Pc_base)
% BIAS_RESID  lsqnonlin 殘差:g_j 用最小二乘 profile 掉後,各模擬的殘差向量堆疊。
    A = build_A(x(1)*1e-3, make_Pc(x(2:18), Pc_base), P);   % 3Np x 6
    M = A.' * A;                                            % 6 x 6
    r = [];
    for j = 1:size(Bstack,2)
        gj = M \ (A.' * Bstack(:,j));                       % 該模擬最佳電荷向量
        r  = [r; A*gj - Bstack(:,j)];                      %#ok<AGROW> 殘差(平方和 = J)
    end
end

function A = build_A(ell, Pc, P)
% BUILD_A  堆疊空間函數矩陣 A(3Np x 6):第 k 欄為各點 (pbar-Pc_k)/||.||^3 沿點堆疊。
    Np   = size(P, 1);
    pbar = P / ell;
    A    = zeros(3*Np, 6);
    for k = 1:6
        D  = pbar - Pc(:,k).';
        r3 = sum(D.^2, 2).^1.5;
        Sk = D ./ r3;
        A(:,k) = reshape(Sk.', 3*Np, 1);
    end
end

function Pc = make_Pc(e17, Pc_base)
% MAKE_PC  由 17 自由偏置組裝 Pc = Pc_base + E;e6z = e1x-e2x+e3y-e4y+e5z(約束)。
    E = zeros(3, 6);
    E(:,1) = e17(1:3);   E(:,2) = e17(4:6);
    E(:,3) = e17(7:9);   E(:,4) = e17(10:12);
    E(:,5) = e17(13:15);
    E(1,6) = e17(16);    E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);
    Pc = Pc_base + E;
end

function [KbarI, gB] = gauge_KI(ell, Pc, P, Bstack, F)
% GAUGE_KI  Step 6-8:g_j、G、H = gB*K̄_I,用 gauge k̄_I(1,1)=5/6 拆出。
    A   = build_A(ell, Pc, P);
    M   = A.' * A;
    C   = A.' * Bstack;
    G   = M \ C;
    H   = G * F.' / (F * F.');
    g11 = G(1,1);                                           % 文件 step8 用 g11(= h11,此標定)
    KbarI = (5/(6*g11)) * H;
    gB    = (6/5) * H(1,1);
end
