function d = solve_d(P, B, Pc, ell_hat, Vmat, exc_sign, N_I)
% SOLVE_D  解 Hall-sensor 模型每極常數 d（閉式解，無增益 g_H）。
% -------------------------------------------------------------------------
% 模型：b_ij = S_i · V_j · d           （V_j = diag(Vmat(:,j))；無 g_H）
%   對所有工作點 i、所有模擬 j 做最小二乘 min_d Σ_ij‖S_i V_j d − b_ij‖²，
%   其閉式解（本檔直接實作此式，內外層雙重加總集中在這裡）：
%
%     d = ( Σ_{j=1}^{N_I} V_j ( Σ_{i=1}^{N_p} S_iᵀ S_i ) V_j )⁻¹
%         ( Σ_{j=1}^{N_I} V_j ( Σ_{i=1}^{N_p} S_iᵀ b_ij ) )
%
%   內層（對 i）：M = Σ_i S_iᵀ S_i（6×6，與 j 無關）、c_j = Σ_i S_iᵀ b_ij。
%   外層（對 j）：A = Σ_j V_j M V_j、rhs = Σ_j V_j c_j。最後 d = A⁻¹ rhs。
%
% 符號：b_ij 取 all-source 物理場 = exc_sign(j)·(−B(i,:,j))
%       （−B = 物理場；×exc_sign 翻下極激發 P1/P3/P6，與 Vmat 的 all-source 一致）。
%
% 輸入：
%   P        Np×3 工作點 [m]
%   B        Np×3×N_I 記錄的場（負號版 −B^FEM）
%   Pc       3×6 電荷位置（極尖單位方向 d̂）
%   ell_hat  特徵長度 ℓ̂ [m]（build_S 用）
%   Vmat     6×6 all-source sensor 電壓（extract_Vmat）
%   exc_sign 1×6 all-source 翻號（下極激發 = −1）
%   N_I      模擬數（6）
% 輸出：
%   d        6×1 每極常數（no-gain，直接照公式；無 g_H）
%
%   (calls build_S)
% -------------------------------------------------------------------------
    Np  = size(P,1);                                         % 工作點數 N_p
    A   = zeros(6,6);                                        % 外層累加器：A = Σ_j V_j (Σ_i S_iᵀS_i) V_j（法矩陣）
    rhs = zeros(6,1);                                        % 外層累加器：rhs = Σ_j V_j (Σ_i S_iᵀb_ij)（右端）
    for j = 1:N_I                                            % === 外層 Σ_{j=1}^{N_I}：第 j 次模擬 ===
        Vj = diag(Vmat(:,j));                               % V_j = diag(該模擬的 6 個 all-source sensor 電壓)
        Mj = zeros(6,6);                                     % 內層累加器：Mj = Σ_i S_iᵀ S_i（每 j 相同，公式裡寫在 j 內）
        cj = zeros(6,1);                                     % 內層累加器：cj = Σ_i S_iᵀ b_ij（該模擬的右端）
        for i = 1:Np                                         % === 內層 Σ_{i=1}^{N_p}：該模擬累加所有工作點 ===
            Si  = build_S(P(i,:), ell_hat, Pc);             % 該點 3×6 空間核 S_i（在 ℓ̂ 下）
            bij = exc_sign(j) * (-squeeze(B(i,:,j)).');     % all-source 物理場 b_ij（−B 轉物理、×exc_sign 翻下極）
            Mj  = Mj + Si.' * Si;                           % 累加 S_iᵀ S_i 進 Mj
            cj  = cj + Si.' * bij;                          % 累加 S_iᵀ b_ij 進 cj
        end
        A   = A   + Vj * Mj * Vj;                           % 累加 V_j (Σ_i S_iᵀS_i) V_j 進 A
        rhs = rhs + Vj * cj;                                % 累加 V_j (Σ_i S_iᵀb_ij) 進 rhs
    end
    d = A \ rhs;                                            % 閉式解 d = A⁻¹ rhs（無 g_H，直接照公式；不顯式求逆）
end
