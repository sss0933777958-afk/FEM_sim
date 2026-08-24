function compare_gain_iso_fix_vs_bias()
% COMPARE_GAIN_ISO_FIX_VS_BIAS  bias 相對 fix 的逐點 C/kappa RMS 相對誤差（R≤150µm 球）。
%   在 R≤150µm 球內每個真實 FEM 節點 p，各自算 fix（無 bias）與 bias 模型的
%     C(p) = ∏σ_k = σ₁σ₂σ₃、kappa(p) = σ_min/σ_max = σ₃/σ₁，^BG_I(p)=S(p)·Ĥ_I（singular_value.pdf）。
%   回報 RMS 相對誤差（參考=沒 bias / fix）：
%     relRMS_X = sqrt( Σ_i (X_bias,i − X_fix,i)² / Σ_i X_fix,i² ),  X∈{C, kappa}.
%   ★ 兩模型都在 **actuator 框**同一組節點上算：R_act·dhat = Pc_base（load_coils_actuator 已 assert）
%     → fix 電荷在 actuator 框 = Pc_base（在軸）、bias = Pc=make_Pc(ê)（離軸）；只差 (ℓ̂, Ĥ, 電荷位置)。
%   純數值、印 console；不存檔、不畫圖。

    % [MODIFIED 2026-08-08] 脫離 backup（規則 no-backup-data）→ live config。
    CALROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CALROOT,'function'), fullfile(CALROOT,'utils'), fullfile(CALROOT,'common_path'));
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');  % ansys_path
    here  = fileparts(mfilename('fullpath'));
    nofix = fileparts(fileparts(here));                        % .../current_base
    calroot = fileparts(nofix);                               % .../Calibration_using_FEM_modeling
    addpath(fullfile(nofix,'code','main_function'));          % load_coils_actuator/select_ball/make_Pc

    cnst = model_config('long2016_hexapole_halfcut','tip40um');
    apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];

    %% ---- R≤150µm 球內真實節點（actuator 框，8774）----
    D = load_coils_actuator('long2016_hexapole_halfcut', cnst, apdl_to_paper_idx, 'all', 'gap_200um');
    [P, ~, npts] = select_ball(D, 150e-6);                    % P: Np×3 [m] actuator 框
    Pc_base = D.Pc_base;                                      % 3×6 在軸單位晶格（actuator）

    %% ---- 兩模型參數 ----
    Sf = load(fullfile(calroot,'current_base','data','fit_fixl_R150um_gap200um_mueq.mat'), 'ell','gB','Khat');
    Sb = load(fullfile(nofix,        'data','fit_bias_R150um_gap200um_mueq.mat'), 'ell','gB','Khat','e_hat');
    ellf = Sf.ell*1e-6;   Hf = Sf.gB*Sf.Khat;   Pcf = Pc_base;                       % fix：在軸電荷
    ellb = Sb.ell*1e-6;   Hb = Sb.gB*Sb.Khat;   Pcb = make_Pc(Sb.e_hat, Pc_base);    % bias：離軸電荷

    %% ---- 逐節點 C / kappa ----
    Cvf = zeros(npts,1); kapf = zeros(npts,1);
    Cvb = zeros(npts,1); kapb = zeros(npts,1);
    for i = 1:npts
        p = P(i,:).';
        svf = svd(((p/ellf - Pcf) ./ (vecnorm(p/ellf - Pcf).^3)) * Hf);
        svb = svd(((p/ellb - Pcb) ./ (vecnorm(p/ellb - Pcb).^3)) * Hb);
        Cvf(i)=prod(svf); kapf(i)=svf(3)/svf(1);
        Cvb(i)=prod(svb); kapb(i)=svb(3)/svb(1);
    end

    RMSE_C     = sqrt( mean((Cvb-Cvf).^2) );                 % 絕對 RMSE = sqrt(mean(Δ²))
    RMSE_kappa = sqrt( mean((kapb-kapf).^2) );
    relRMS_C     = sqrt( sum((Cvb-Cvf).^2) / sum(Cvf.^2) );  % RMS 相對誤差（÷ Σ fix²）
    relRMS_kappa = sqrt( sum((kapb-kapf).^2) / sum(kapf.^2) );

    fprintf('\n=== bias vs no-bias（fix）逐點 C/kappa，R≤150µm 球內 %d 節點 ===\n', npts);
    fprintf('sanity  mean C     : fix=%.4g  bias=%.4g (mT/A)^3\n', mean(Cvf), mean(Cvb));
    fprintf('sanity  mean kappa : fix=%.4f  bias=%.4f\n',           mean(kapf), mean(kapb));
    fprintf('C     逐點差：max|Δ|=%.4g (mT/A)^3\n', max(abs(Cvb-Cvf)));
    fprintf('kappa 逐點差：max|Δ|=%.4f\n',           max(abs(kapb-kapf)));
    fprintf('----\n');
    fprintf('RMSE           C = %.4g (mT/A)^3 | kappa = %.4f\n', RMSE_C, RMSE_kappa);
    fprintf('RMS 相對誤差   C = %.4f%%        | kappa = %.4f%%\n', 100*relRMS_C, 100*relRMS_kappa);
end

% ---- local make_Pc (self-contained; charge lattice Pc = Pc_base + E(e17)) ----
function Pc = make_Pc(e17, Pc_base)
    E = zeros(3,6);
    E(:,1)=e17(1:3); E(:,2)=e17(4:6); E(:,3)=e17(7:9); E(:,4)=e17(10:12); E(:,5)=e17(13:15);
    E(1,6)=e17(16);  E(2,6)=e17(17);  E(3,6)=e17(1)-e17(4)+e17(8)-e17(11)+e17(15);
    Pc = Pc_base + E;
end
