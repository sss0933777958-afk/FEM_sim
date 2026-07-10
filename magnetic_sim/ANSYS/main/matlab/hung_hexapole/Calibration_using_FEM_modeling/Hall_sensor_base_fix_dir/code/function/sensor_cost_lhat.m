function J = sensor_cost_lhat(ell, P, B, Vmat, exc_sign, Pc, N_I)
% SENSOR_COST_LHAT  Hall-sensor 模型「消去 d 後」的 ℓ̂ profiled cost（照 cost_function.pdf）。
% -------------------------------------------------------------------------
% 用途（步驟1：找 ℓ̂）：
%   把每極常數 d 在每次激發 j 內以最小二乘閉式消掉，得到只剩 ℓ̂ 的成本：
%
%     J(ℓ̂) = Σ_{j=1}^{N_I} [ Σ_{i=1}^{N_p} b_ijᵀb_ij
%                            − (Σ_i V_jᵀS_iᵀb_ij)ᵀ (Σ_i V_jᵀS_iᵀS_iV_j)⁻¹ (Σ_i V_jᵀS_iᵀb_ij) ]
%
%   對 ℓ̂ 一維最小化此 J 即得特徵長度（main.m 用 fminbnd 包它）。
%   注意：此 cost 的反矩陣「包在 Σ_j 內」(per-excitation profile)，是文件 cost_function.pdf
%   的原樣；最終 d（步驟2 solve_d）才是 6×1 shared（反矩陣在 Σ_j 外）。兩步刻意分開。
%
% 與 solve_d 的符號一致：
%   b_ij = exc_sign(j)·(−B(i,:,j))（−B 轉物理場、×exc_sign 翻下極 all-source）；
%   V_j  = diag(Vmat(:,j))（Vmat 已含 all-source 翻號）。
%   S_i  = build_S(P_i, ℓ̂, Pc) 只透過 ℓ̂ 進來（p̄=p/ℓ̂）→ J 是 ℓ̂ 的純量函數。
%
% 輸入：
%   ell      純量 ℓ̂ [m]（被 fminbnd 掃的變數）
%   P        Np×3 工作點 [m]（WP 框，R≤R_select 球內）
%   B        Np×3×N_I 記錄的場（負號版 −B^FEM，與 solve_d 同）
%   Vmat     6×6 all-source sensor 電壓
%   exc_sign 1×6 all-source 翻號（下極激發 = −1）
%   Pc       3×6 電荷位置（極尖單位方向）
%   N_I      模擬數（6）
% 輸出：
%   J        純量 cost(ℓ̂) [T²]
%
%   (calls build_S)
% -------------------------------------------------------------------------
    Np    = size(P,1);
    M     = zeros(6,6);          % 內層 Σ_i S_iᵀS_i（與 j 無關，只算一次）
    c     = zeros(6,N_I);        % 內層 c_j = Σ_i S_iᵀb_ij（每激發一欄）
    bnorm = zeros(1,N_I);        % Σ_i b_ijᵀb_ij（資料項，每激發一個）
    for i = 1:Np                                     % 逐工作點：建一次 S_i，分配進 M / 各 c_j
        Si    = build_S(P(i,:), ell, Pc);           % 3×6 空間核（在 ℓ̂ 下）
        M     = M + Si.'*Si;                         % 累加 S_iᵀS_i
        for j = 1:N_I
            bij      = exc_sign(j) * (-squeeze(B(i,:,j)).');   % all-source 物理場 b_ij（3×1）
            c(:,j)   = c(:,j)   + Si.'*bij;          % 累加 S_iᵀb_ij
            bnorm(j) = bnorm(j) + bij.'*bij;         % 累加 b_ijᵀb_ij
        end
    end
    J = 0;
    for j = 1:N_I                                    % 外層 Σ_j：每激發各 profile 一個 d_j 後的殘差
        Vj  = diag(Vmat(:,j));
        Mjp = Vj * M * Vj;                           % Σ_i V_jᵀS_iᵀS_iV_j（6×6）
        cjp = Vj * c(:,j);                           % Σ_i V_jᵀS_iᵀb_ij（6×1）
        J   = J + bnorm(j) - cjp.' * (Mjp \ cjp);    % ‖b_j‖² − c_jᵀM_j⁻¹c_j
    end
end
