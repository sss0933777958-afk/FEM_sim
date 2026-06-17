%% calib_fem_bias.m
%  ============================================================================
%  依「no_fix_l.pdf」實作【Calibration using FEM modeling, 18-parameter (bias) model】。
%
%  本檔同時跑文件的兩個模型(每個步驟文件都並列 a/b 兩版):
%    (A) One-parameter model   : 只 fit 特徵長度 ell_hat,6 顆電荷鎖在 ±ell_hat 對稱格點。
%    (B) 18-parameter model    : ell_hat + 1x17 偏置向量 e_hat(放掉電荷格點的剛性),
%                                第 6 顆電荷的 z 分量 e6z 由其餘分量約束 → 17 自由 + ell = 18 參數。
%
%  與既有 calib_fem.m 的關鍵差異:
%    * calib_fem.m 留在「量測(measure)座標」,用 Pc=dhat(實際極尖方向)避開旋轉。
%    * 文件的 bias 是逐分量定義在「致動器(actuator)座標」(無 bias 時 ell 固定的方向 =
%      三條互相正交的對極軸、原點在工作空間中心),所以本檔先把節點座標與 B 旋進 actuator 框,
%      理想格點就成為文件的 Pc = [+u -u +v -v +w -w](單位軸),再套 bias 與 e6z 約束。
%    * 因點電荷核對旋轉協變,one-parameter 段的 ell_hat 與 calib_fem.m 完全相同(旋轉不變)→
%      可當「座標/旋轉正確」的 sanity(預期 ell_hat ~ 0.856 mm)。
%
%  資料來源(使用者拍板):standard mesh 的 'all' dataset(README 記載全網格 ~494,873 節點);
%  Step 1 仍在 work space 用近場球 R_select 選點(點電荷是近場模型,遠場無法被 6 顆電荷表示)。
%
%  只讀 magnetic_sim/ANSYS/main/ANSYS_data、只寫 magnetic_sim/ANSYS/main/MATLAB_data/charge_fit;不搬不改任何資料夾結構。
%  ============================================================================
clear; clc;                                                  % 清空工作區與命令視窗

% ---- paths(用 helper,不硬寫磁碟根路徑)----------------------------------
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');        % ansys_path / matlab_path
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\hexapole-long2016\analysis');  % mt_constants / import_ansys_data / filter_iron_nodes
model = 'long2016_hexapole_halfcut';                         % 物理模型(model-first 結構的第一層)

% ---- 頂層參數 --------------------------------------------------------------
cnst       = mt_constants();        % 幾何常數(R_norm、極尖位置、SPH_OFST 等)
R_select   = 150e-6;                % <== 單一旋鈕:Step 1 近場 work space 取點半徑 [m]
                                    %     預設 150e-6(文件驗證值,ell_hat~0.856;此球內 'all'='wp' 節點集相同)。
                                    %     ⚠ 點電荷是近場模型:勿放到涵蓋全 494k(遠場主導會讓 ell_hat 失真)。
ell_design = cnst.R_norm;           % Step 5 ell_hat 初始猜 [m] = 尖端到 WP 距離(500 um)
N_I        = 6;                     % FEM 模擬次數 = 6 個單線圈解
apdl_to_paper_idx = [1,3,6,5,2,4];  % APDL coil j → paper pole 索引(coil1→P1, coil2→P3, ...)
dataset    = 'all';                 % <== 使用者拍板:standard mesh 的 'all'(全網格)dataset

%% ===========================================================================
%  Step 2a - 致動器座標系 [u v w] = I(文件:單位正交基底);原點在工作空間中心。
%  由本設計三條互相正交的「對極軸」實現:u=P1 尖端方向、v=P3、w=P5(已驗證互相正交)。
%  R(3x3,rows=actuator 基底,以 measure 框表示)把 measure→actuator:R*dhat = 理想格點。
%  ===========================================================================
tip   = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];  % 3x6 極尖位置(WP/measure 座標)[m]
dhat  = tip ./ vecnorm(tip);                                    % 3x6 極尖單位方向(measure 框)
uhat  = dhat(:,1);  vhat = dhat(:,3);  what = dhat(:,5);        % actuator 軸:P1=+x_a / P3=+y_a / P5=+z_a
                                                               %   (per notation-glossary / dissertation Fig.2.1:
                                                               %    對極 (P1,P2)/(P3,P4)/(P5,P6),挑一顆/對 = +x_a/+y_a/+z_a)
R     = [uhat, vhat, what].';                                   % 3x3 measure→actuator 旋轉

Pc_base = [ 1 -1  0  0  0  0;                                   % 文件理想格點 [+u -u +v -v +w -w](paper 序 P1..P6)
            0  0  1 -1  0  0;
            0  0  0  0  1 -1];
assert(abs(det(R)-1) < 1e-9, 'R 必須是 proper rotation(det=+1)');                    % 右手系自洽
assert(max(abs(R*dhat - Pc_base), [], 'all') < 1e-9, 'R*dhat 必須等於理想格點 Pc_base');% 旋轉/指派自洽
fprintf('Step2: actuator 旋轉 R 建好(det=%.6f);R*dhat == [+u -u +v -v +w -w] 已驗證\n', det(R));

%% ===========================================================================
%  Step 1 - 在 work space 選 N_p 個點 p_i。
%  載 coil1 'all' → 過濾鐵芯取 air → WP 為原點 → |p|<=R_select 近場球 → 旋進 actuator 框。
%  6 顆 coil 共用同一 mesh(節點集/順序一致),故 air 遮罩與 insel 在 6 顆之間通用。
%  ===========================================================================
d1   = import_ansys_data(ansys_path(model,'coil1','standard'), dataset, 'coil1');     % coil1 全網格(x,y,z,bx,by,bz)
air1 = filter_iron_nodes(d1.x,d1.y,d1.z,cnst,struct('visualize',false));             % air 遮罩(排除鐵芯)
zwp1 = d1.z - cnst.SPH_OFST;                                                         % ANSYS z → WP 為原點
P_all = [d1.x(air1), d1.y(air1), zwp1(air1)];                                        % air 節點座標(measure/WP 框)[m]

insel  = sum(P_all.^2, 2) <= R_select^2;                                             % 近場球選點(measure 框做,filter 認 measure 幾何)
P_meas = P_all(insel, :);                                                            % (Np x 3) 選出的 work-space 點(measure)
P      = (R * P_meas.').';                                                           % (Np x 3) 旋進 actuator 框
Np     = size(P, 1);                                                                 % N_p
fprintf('Step1: standard ''%s'' dataset,R_select = %g um → N_p = %d 個近場點\n', dataset, R_select*1e6, Np);

%% ===========================================================================
%  Step 2 - F = [I_1 ... I_{N_I}](6xN_I)須 rank 6。單線圈:第 j 次激發 APDL coil j(1A),
%  對應 paper pole apdl_to_paper_idx(j) → I_j = e_{pj(j)},F 為置換矩陣。 I = 1A(= FEM 激發,
%  per fit-current-matches-sim;不可塞操作電流 0.6A)。
%  ===========================================================================
F = zeros(6, N_I);                                           % 6x6 電流矩陣(pole 空間 P1..P6)
for j = 1:N_I, F(apdl_to_paper_idx(j), j) = 1; end           % 第 j 次:paper pole pj(j) 通 1A
assert(rank(F) == 6, 'F 必須 rank 6');                       % rank 6 要求

%% ===========================================================================
%  Step 3 - 記錄 b_ij = b^FEM(p_i, I_j)。符號:Bn = -B^FEM(整體變號,all-source 慣例,
%  與 calib_fem 一致;ell_hat 對整體變號不變)。旋進 actuator 框後堆成 Bstack(3Np x N_I)以利向量化。
%  ===========================================================================
B = zeros(Np, 3, N_I);                                       % b_ij 容器(點, xyz, 模擬;actuator 框)
for j = 1:N_I                                                % 逐一掃 6 個單線圈模擬
    if j == 1                                                % coil1 已載 → 重用
        dj = d1;  airj = air1;
    else
        cn = sprintf('coil%d', j);                          % 'coil2'..'coil6'
        dj = import_ansys_data(ansys_path(model,cn,'standard'), dataset, cn);          % 載 coil j 全網格
        airj = filter_iron_nodes(dj.x,dj.y,dj.z,cnst,struct('visualize',false));      % air 遮罩(幾何相同→相同)
    end
    Bj_all = -[dj.bx(airj), dj.by(airj), dj.bz(airj)];      % (Nair x 3) coil j 場 [T];取負(all-source)
    Bj_sel = Bj_all(insel, :);                              % 只留近場選點(measure 框)
    B(:,:,j) = (R * Bj_sel.').';                            % 旋進 actuator 框 → b_ij
end
Bstack = zeros(3*Np, N_I);                                  % 堆疊:每欄 = 該模擬的 b_ij(點-major、xyz 交錯)
for j = 1:N_I, Bstack(:,j) = reshape(B(:,:,j).', 3*Np, 1); end
fprintf('Step3: 已記錄 b_ij(%d 模擬 x %d 點),旋入 actuator 框\n', N_I, Np);

%% ===========================================================================
%  Step 4/5 - 成本函數 + 最小化。
%  成本 J 用最小二乘把每次模擬的 6x1 電荷向量 g_j 剖面化(profile out),
%  只剩非線性參數(ell,或 ell+e)。(local 函式 costJ / build_A / make_Pc 在檔尾。)
%  ===========================================================================

% --- (A) One-parameter model: Minimize J(ell) ⟹ ell_hat(初始猜 ell_design),e = 0 ---
obj1   = @(x) costJ(x, Pc_base, P, Bstack);                 % 成本對 ell 的一維函式(Pc 固定理想格點)
fopt1  = optimset('TolX',1e-9,'Display','iter');           % 容差 1 nm
[ell_hat_1param, J1] = fminbnd(obj1, 0.2e-3, 2.0e-3, fopt1);% 一維最小化
fprintf('Step5a: one-param ell_hat = %.4f mm (J=%.6e) [sanity vs calib_fem ~0.856 mm]\n', ell_hat_1param*1e3, J1);

% --- (B) 18-parameter model: Minimize J(ell, e) ⟹ ell_hat, e_hat(初始猜 ell_design, e = 0) ---
%  決策變數 x = [ell_mm; e(17)];ell 用 mm 尺度化改善條件數。e 為無因次(加在單位格點上)。
objx = @(x) costJ(x(1)*1e-3, make_Pc(x(2:18), Pc_base), P, Bstack);
x0   = [ell_design*1e3; zeros(17,1)];                       % 初始猜:ell_design(mm)、e = 0
% R2025b 的 fminunc 在收斂後組「結束訊息」時會踩到 'optim:fminusub' 訊息目錄載入 bug
% (最佳化本身正常收斂,只壞在最後組字串)。用 OutputFcn 每代存當前 x,再 try/catch 把
% 收斂解救回來,繞過此 toolbox bug(fminbnd 無此問題)。
global FMINUNC_LAST_X; FMINUNC_LAST_X = x0;                 %#ok<GVMIS> 每代由 capture_x 更新
optx = optimoptions('fminunc','Display','iter','Algorithm','quasi-newton', ...
                    'MaxFunctionEvaluations',5e4,'MaxIterations',3e3, ...
                    'OptimalityTolerance',1e-12,'StepTolerance',1e-12, ...
                    'FiniteDifferenceType','central', ...
                    'OutputFcn',@capture_x);                % 中央差分梯度(18 維)+ 每代存 x
try
    [xfit, J18] = fminunc(objx, x0, optx);                  % 18 維最小化
catch ME                                                    % 繞過 fminunc 結束訊息 bug
    if contains(ME.message,'fminusub') || contains(ME.message,'message catalog')
        xfit = FMINUNC_LAST_X(:);                           % 收斂解(最後一代 OutputFcn 存下)
        J18  = objx(xfit);                                  % 重算該點成本
        fprintf('[note] 繞過 fminunc 結束訊息 bug(optim:fminusub),已救回收斂解 J=%.6e\n', J18);
    else
        rethrow(ME);                                        % 其他錯誤照常拋出
    end
end
ell_hat = xfit(1)*1e-3;                                     % ell_hat [m]
e_hat   = xfit(2:18);                                       % 17 自由偏置
e6z     = e_hat(1) - e_hat(4) + e_hat(8) - e_hat(11) + e_hat(15);  % 約束:e6z = e1x-e2x+e3y-e4y+e5z
Pc_18   = make_Pc(e_hat, Pc_base);                          % 18-param 電荷格點(actuator 框)
E_hat   = Pc_18 - Pc_base;                                  % 偏置矩陣 E(3x6,含算出的 e6z)
fprintf('Step5b: 18-param ell_hat = %.4f mm (J=%.6e),e6z(約束) = %+.4e\n', ell_hat*1e3, J18, e6z);

%% ===========================================================================
%  Step 6-8 - g_j、G、g_B*K̄_I(對兩個模型各做一次)。
%    g_j = (Σ S^T S)^{-1}(Σ S^T b_ij);G = [g_1..g_{N_I}];
%    H = G F^T (F F^T)^{-1} = gB_hat*K̄_I;gB_hat = (6/5) h11;K̄_I = 5/(6 g11) H。
%  其中 g11 = G(1,1)(模擬1 電荷向量 g_1 的第一分量)、h11 = H(1,1)。本標定 coil1 激發 P1
%  且 F 置換不動第1欄 → H(:,1)=G(:,1) → g11 = h11,故文件兩式自洽(K̄_I = H/gB_hat),
%  K̄_I(1,1)=5/6。文件寫 g11 正確(非筆誤)。
%  ===========================================================================
[KbarI_1,  gB_1,  G1, H1]  = gauge_KI(ell_hat_1param, Pc_base, P, Bstack, F);  % one-param
[KbarI_18, gB_18, G18, H18] = gauge_KI(ell_hat,        Pc_18,   P, Bstack, F);  % 18-param

% all-source 顯示:翻上極 P2/P4/P5 欄(純符號重貼,不改場;對齊 fit_KI_R150)
coil_sign     = [1 -1 1 -1 -1 1];                          % paper P1..P6(下極+1、上極-1)
KbarI_1_disp  = KbarI_1  .* coil_sign;                     % 顯示用(對角全正)
KbarI_18_disp = KbarI_18 .* coil_sign;

%% ===========================================================================
%  殘差量化 - canonical NRMSE(worst-of-6,R_select 球內;重建 B_model=S*g_j vs all-source B_FEM)。
%  ===========================================================================
[nrmse_1,  ev_1]  = nrmse_worst(ell_hat_1param, Pc_base, P, B);   % one-param
[nrmse_18, ev_18] = nrmse_worst(ell_hat,        Pc_18,   P, B);   % 18-param

%% ===========================================================================
%  報告
%  ===========================================================================
fprintf('\n================= no_fix_l.pdf 校正結果(actuator 框, N_p=%d, R_select=%g um)=================\n', Np, R_select*1e6);
fprintf('  (A) one-parameter : ell_hat = %.4f mm,  gB_hat = %.4e,  J = %.4e,  NRMSE(worst) = %.3f%%\n', ...
        ell_hat_1param*1e3, gB_1, J1, nrmse_1);
fprintf('  (B) 18-parameter  : ell_hat = %.4f mm,  gB_hat = %.4e,  J = %.4e,  NRMSE(worst) = %.3f%%\n', ...
        ell_hat*1e3, gB_18, J18, nrmse_18);
fprintf('  J 改善 = %.4e → %.4e (%.1f%%);NRMSE 改善 = %.3f%% → %.3f%%\n', ...
        J1, J18, (J1-J18)/J1*100, nrmse_1, nrmse_18);

fprintf('\n  e_hat(17 自由偏置, actuator 框, 無因次):\n');
elab = {'e1x','e1y','e1z','e2x','e2y','e2z','e3x','e3y','e3z','e4x','e4y','e4z','e5x','e5y','e5z','e6x','e6y'};
for q = 1:17, fprintf('    %-4s = %+.4e\n', elab{q}, e_hat(q)); end
fprintf('    e6z (約束) = %+.4e\n', e6z);

fprintf('\n  K̄_I (one-param, all-source 顯示) =\n');
for i=1:6, fprintf('    % .4f % .4f % .4f % .4f % .4f % .4f\n', KbarI_1_disp(i,:)); end
fprintf('  K̄_I (18-param, all-source 顯示) =\n');
for i=1:6, fprintf('    % .4f % .4f % .4f % .4f % .4f % .4f\n', KbarI_18_disp(i,:)); end
fprintf('==================================================================================================\n');

%% ===========================================================================
%  存檔(只寫 MATLAB_data/charge_fit;不動任何結構)
%  ===========================================================================
out_dir = matlab_path(model, 'charge_fit');                 % .../kuo/MATLAB_data/<model>/charge_fit
save(fullfile(out_dir, 'calib_bias.mat'), ...
     'R', 'Pc_base', 'R_select', 'Np', 'dataset', 'apdl_to_paper_idx', 'F', 'coil_sign', ...
     'ell_hat_1param', 'KbarI_1', 'gB_1', 'J1', 'nrmse_1', 'ev_1', ...
     'ell_hat', 'e_hat', 'e6z', 'E_hat', 'Pc_18', 'KbarI_18', 'gB_18', 'J18', 'nrmse_18', 'ev_18');
fprintf('已存 %s\n', fullfile(out_dir, 'calib_bias.mat'));

%% ===========================================================================
%  本地函式
%  ===========================================================================
function stop = capture_x(x, ~, ~)
% CAPTURE_X  fminunc OutputFcn:每代把當前 x 存進 global,讓 R2025b 結束訊息 bug 後仍能救回收斂解。
    global FMINUNC_LAST_X
    FMINUNC_LAST_X = x;
    stop = false;
end

function A = build_A(ell, Pc, P)
% BUILD_A  堆疊空間函數矩陣 A(3Np x 6):第 k 欄為各點 (pbar-Pc_k)/||pbar-Pc_k||^3 的
%   x,y,z 沿點堆疊(點-major、xyz 交錯),與 Bstack 的堆疊一致。pbar = p/ell。
    Np   = size(P, 1);                                       % 點數
    pbar = P / ell;                                          % Np x 3 正規化點
    A    = zeros(3*Np, 6);                                   % 預配置
    for k = 1:6                                              % 逐電荷(欄)
        D  = pbar - Pc(:,k).';                               % Np x 3:pbar - Pc_k(廣播)
        r3 = sum(D.^2, 2).^1.5;                              % Np x 1:||.||^3
        Sk = D ./ r3;                                        % Np x 3:該電荷的核
        A(:,k) = reshape(Sk.', 3*Np, 1);                    % 沿點堆疊成 3Np x 1
    end
end

function Pc = make_Pc(e17, Pc_base)
% MAKE_PC  由 17 自由偏置組裝 18-param 電荷格點 Pc = Pc_base + E(actuator 框)。
%   e17 = [e1x e1y e1z e2x e2y e2z e3x e3y e3z e4x e4y e4z e5x e5y e5z e6x e6y]
%   第 6 顆電荷的 z 由約束:e6z = e1x - e2x + e3y - e4y + e5z。
    E = zeros(3, 6);                                         % 偏置矩陣
    E(:,1) = e17(1:3);   E(:,2) = e17(4:6);                  % e1, e2
    E(:,3) = e17(7:9);   E(:,4) = e17(10:12);                % e3, e4
    E(:,5) = e17(13:15);                                     % e5
    E(1,6) = e17(16);    E(2,6) = e17(17);                   % e6x, e6y
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);   % e6z(約束)
    Pc = Pc_base + E;                                        % 文件 page2 的 Pc_c(e)
end

function J = costJ(ell, Pc, P, Bstack)
% COSTJ  剖面化最小二乘成本:每次模擬 j 把 g_j = M^{-1} c_j 解掉後的殘差平方和之總和。
%   M = A^T A = Σ_i S_i^T S_i;  c_j = A^T b_j = Σ_i S_i^T b_ij;  bb_j = Σ_i b_ij^T b_ij。
    A  = build_A(ell, Pc, P);                               % 3Np x 6
    M  = A.' * A;                                            % 6 x 6
    C  = A.' * Bstack;                                       % 6 x N_I(各欄 c_j)
    bb = sum(Bstack.^2, 1);                                  % 1 x N_I(各 bb_j)
    J  = 0;                                                  % 累加成本
    for j = 1:size(Bstack, 2)
        J = J + (bb(j) - C(:,j).' * (M \ C(:,j)));          % bb - c^T M^{-1} c
    end
end

function [KbarI, gB, G, H] = gauge_KI(ell, Pc, P, Bstack, F)
% GAUGE_KI  Step 6-8:在給定 (ell,Pc) 下求 g_j、G、H = gB*K̄_I,並用 gauge k̄_I(1,1)=5/6 拆出。
    A   = build_A(ell, Pc, P);                               % 3Np x 6
    M   = A.' * A;                                           % 6 x 6
    C   = A.' * Bstack;                                      % 6 x N_I(c_j)
    G   = M \ C;                                             % 6 x N_I(各欄 g_j = M^{-1} c_j)
    H   = G * F.' / (F * F.');                               % 6 x 6 = gB*K̄_I
    g11 = G(1,1);                                            % G 的 (1,1) = g_1 第一分量(文件 step8 用 g11)
    h11 = H(1,1);                                            % H 的 (1,1)
    assert(abs(g11-h11) < 1e-12*max(1,abs(h11)), 'g11≠h11:F 第1欄被置換移動,gauge 兩式不自洽');
    gB  = (6/5) * h11;                                       % gB_hat = (6/5) h11(文件)
    KbarI = (5/(6*g11)) * H;                                 % K̄_I = 5/(6 g11) H(文件;此處 g11=h11 → K̄_I(1,1)=5/6)
end

function [nrmse, ev] = nrmse_worst(ell, Pc, P, B)
% NRMSE_WORST  canonical 相對 RMSE(每 coil:sqrt(mean||Bm-Bf||^2)/max||Bf|| x100),取 6 顆最差。
%   Bm = A*g_j(剖面化最佳電荷向量重建場);Bf = 該模擬的 all-source b_ij(= 餵進來的 B)。
    Np  = size(P,1); N_I = size(B,3);                        % 點/模擬數
    A   = build_A(ell, Pc, P);                               % 3Np x 6
    M   = A.' * A;                                           % 6 x 6
    ev  = zeros(1, N_I);                                     % 各 coil 的 NRMSE
    for j = 1:N_I
        Bf      = B(:,:,j);                                  % Np x 3
        bstackj = reshape(Bf.', 3*Np, 1);                    % 3Np x 1
        gj      = M \ (A.' * bstackj);                       % 6 x 1 最佳電荷向量
        Bm      = reshape(A * gj, 3, Np).';                  % Np x 3 重建場
        rmse    = sqrt(mean(sum((Bm - Bf).^2, 2)));          % 場向量 RMSE
        denom   = max(sqrt(sum(Bf.^2, 2)));                  % 該 coil 的 |B| max
        ev(j)   = rmse / denom * 100;                        % 相對 RMSE [%]
    end
    nrmse = max(ev);                                         % worst-of-6
end
