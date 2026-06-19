function [d, gH] = solve_d(Vmat, exc_sign, M, c, ell_hat, cnst, N_I)
% SOLVE_D  解 Hall-sensor 模型的每極常數 d（依 d_final.pdf，含增益 g_H）。
% -------------------------------------------------------------------------
% 模型（PAGE2）：
%   把每極電荷強度用「sensor 電壓 × 每極常數 d」表示：
%       b_ij = g_H · S_i · V_j · d         （V_j = diag(Vmat(:,j))）
%   即第 j 次模擬、第 k 極的等效電荷強度 = g_H · v_jk · d_k。
%   只有 6 個 d_k 是自由參數（V_j 由 FEM 固定、g_H 由 ℓ̂ 決定），
%   對所有點 i、模擬 j 做最小二乘。
%
% 閉式解（對 d 的法方程）：
%   最小化 Σ_ij‖g_H S_i V_j d − b_ij‖² ⇒
%       d = (1/g_H) (Σ_j V_j M V_j)⁻¹ (Σ_j V_j cc_j)
%   其中 M = Σ_i S_iᵀS_i、c_j = Σ_i S_iᵀ b_ij（PAGE1 在 ℓ̂ 已算好、重用）。
%   增益 g_H = k_m/(ℓ̂²μ0) = 1/(4πℓ̂²)（由 Q̂=Φ̂/μ0、b=(k_m/ℓ̂²)S Q̂ 推得）。
%
% 符號：PAGE1 的 c 用的是「負號版場」(b=−B^FEM)；這裡 cc_j = exc_sign(j)·(−c_j)
%       把它轉回物理場、再套 all-source 翻號，與 Vmat 的翻號一致。
%
% 輸入：
%   Vmat     6×6 all-source sensor 電壓（extract_Vmat）
%   exc_sign 1×6 各模擬的 all-source 翻號
%   M        6×6 法矩陣 Σ S_iᵀS_i（PAGE1 在 ℓ̂ 算好）
%   c        6×N_I 右端 Σ S_iᵀ b_ij（負號版場）
%   ell_hat  特徵長度 ℓ̂ [m]
%   cnst     常數（用 k_m、mu_0 算 g_H）
%   N_I      模擬數（6）
% 輸出：
%   d   6×1 每極常數（= d_final，已含增益）
%   gH  純量增益 g_H [1/m²]
%
%   (ported verbatim from calib_fem.m PAGE 2)
% -------------------------------------------------------------------------
    gH = cnst.k_m / (ell_hat^2 * cnst.mu_0);                  % 增益 g_H = k_m/(ℓ̂²μ0) = 1/(4πℓ̂²) [1/m²]
    A = zeros(6,6); rhs = zeros(6,1);                         % 法矩陣 A = ΣV_j M V_j、右端 rhs = ΣV_j cc_j 的累加器
    for j = 1:N_I                                             % 逐模擬 j 累加
        Vj  = diag(Vmat(:,j));                                % V_j = diag(該模擬的 6 個 all-source sensor 電壓)
        ccj = exc_sign(j) * (-c(:,j));                        % all-source 物理右端：−c_j 轉回物理場，再 ×exc_sign 翻號
        A   = A   + Vj * M * Vj;                              % 累加 V_j M V_j（6×6）
        rhs = rhs + Vj * ccj;                                 % 累加 V_j cc_j（6×1）
    end
    d = (A \ rhs) / gH;                                       % d_final = (1/g_H)(ΣVMV)⁻¹(ΣVcc) = 4πℓ̂²·(無增益 d)
end
