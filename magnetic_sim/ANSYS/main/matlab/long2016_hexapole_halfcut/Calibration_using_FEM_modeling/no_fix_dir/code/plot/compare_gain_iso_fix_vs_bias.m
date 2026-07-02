function compare_gain_iso_fix_vs_bias()
% COMPARE_GAIN_ISO_FIX_VS_BIAS  bias 相對 fix 的逐點 gain/iso RMS 相對誤差（R≤150µm 球）。
%   在 R≤150µm 球內每個真實 FEM 節點 p，各自算 fix（無 bias）與 bias 模型的
%     gain(p) = ‖T‖_F = √Σσ²、iso(p) = σ_max/σ_min，T=S(p)·Ĥ_I。
%   回報 RMS 相對誤差（參考=沒 bias / fix）：
%     relRMS_X = sqrt( Σ_i (X_bias,i − X_fix,i)² / Σ_i X_fix,i² ),  X∈{gain, iso}.
%   ★ 兩模型都在 **actuator 框**同一組節點上算：R_act·dhat = Pc_base（load_coils_actuator 已 assert）
%     → fix 電荷在 actuator 框 = Pc_base（在軸）、bias = Pc=make_Pc(ê)（離軸）；只差 (ℓ̂, Ĥ, 電荷位置)。
%   純數值、印 console；不存檔、不畫圖。

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');   % mt_constants/import/filter
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');  % ansys_path
    here  = fileparts(mfilename('fullpath'));
    nofix = fileparts(fileparts(here));                        % .../no_fix_dir
    calroot = fileparts(nofix);                               % .../Calibration_using_FEM_modeling
    addpath(fullfile(nofix,'code','function'));               % load_coils_actuator/select_ball/make_Pc

    cnst = mt_constants();
    apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];

    %% ---- R≤150µm 球內真實節點（actuator 框，8774）----
    D = load_coils_actuator('long2016_hexapole_halfcut', cnst, apdl_to_paper_idx, 'all', 'gap200um_mueq');
    [P, ~, npts] = select_ball(D, 150e-6);                    % P: Np×3 [m] actuator 框
    Pc_base = D.Pc_base;                                      % 3×6 在軸單位晶格（actuator）

    %% ---- 兩模型參數 ----
    Sf = load(fullfile(calroot,'fix_dir','data','fit_fixl_R150um_gap200um_mueq.mat'), 'ell','gB','Khat');
    Sb = load(fullfile(nofix,        'data','fit_bias_R150um_gap200um_mueq.mat'), 'ell','gB','Khat','e_hat');
    ellf = Sf.ell*1e-6;   Hf = Sf.gB*Sf.Khat;   Cf = Pc_base;                 % fix：在軸電荷
    ellb = Sb.ell*1e-6;   Hb = Sb.gB*Sb.Khat;   Cb = make_Pc(Sb.e_hat, Pc_base);   % bias：離軸電荷

    %% ---- 逐節點 gain/iso ----
    gf = zeros(npts,1); isf = zeros(npts,1);
    gb = zeros(npts,1); isb = zeros(npts,1);
    for i = 1:npts
        p = P(i,:).';
        svf = svd(((p/ellf - Cf) ./ (vecnorm(p/ellf - Cf).^3)) * Hf);
        svb = svd(((p/ellb - Cb) ./ (vecnorm(p/ellb - Cb).^3)) * Hb);
        gf(i)=norm(svf); isf(i)=svf(1)/svf(3);
        gb(i)=norm(svb); isb(i)=svb(1)/svb(3);
    end

    RMSE_gain = sqrt( mean((gb-gf).^2) );                 % 絕對 RMSE = sqrt(mean(Δ²))
    RMSE_iso  = sqrt( mean((isb-isf).^2) );
    relRMS_gain = sqrt( sum((gb-gf).^2) / sum(gf.^2) );   % RMS 相對誤差（÷ Σ fix²）
    relRMS_iso  = sqrt( sum((isb-isf).^2) / sum(isf.^2) );

    fprintf('\n=== bias vs no-bias（fix）逐點 gain/iso，R≤150µm 球內 %d 節點 ===\n', npts);
    fprintf('sanity  mean gain: fix=%.4f  bias=%.4f mT/A（對照 σ_tot 18.605/18.612）\n', mean(gf), mean(gb));
    fprintf('sanity  mean iso : fix=%.4f  bias=%.4f      （對照 iso_tot 1.141/1.144）\n', mean(isf), mean(isb));
    fprintf('gain 逐點差：max|Δ|=%.4f mT/A\n', max(abs(gb-gf)));
    fprintf('iso  逐點差：max|Δ|=%.4f\n',      max(abs(isb-isf)));
    fprintf('----\n');
    fprintf('RMSE           gain = %.4f mT/A | iso = %.4f\n', RMSE_gain, RMSE_iso);
    fprintf('RMS 相對誤差   gain = %.4f%%    | iso = %.4f%%\n', 100*relRMS_gain, 100*relRMS_iso);
end
