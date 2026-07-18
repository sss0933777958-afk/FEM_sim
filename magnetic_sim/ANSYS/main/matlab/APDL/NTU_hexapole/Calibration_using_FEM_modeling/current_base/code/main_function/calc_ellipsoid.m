function m = calc_ellipsoid(P, Hhat_I, ell_m, Pc)
% CALC_ELLIPSOID  集總模型在控制範圍（R≤150µm 球）的性能量：flux ellipsoid 體積 C 與均勻度 kappa。
% -------------------------------------------------------------------------
%   對球內每個真實 FEM 節點 p，算致動增益矩陣 ^BG_I(p)=S(p)·Ĥ_I（3×6，mT/A）的 SVD
%   （σ₁≥σ₂≥σ₃），依 singular_value.pdf 定義：
%     C(p) = ∏_{k=1}^{3} σ_k = σ₁σ₂σ₃   （flux-generating ellipsoid 體積；(mT/A)³）
%     kappa(p) = σ_min/σ_max = σ₃/σ₁     （condition number / isotropy；≤1，→1 等向）
%   代表值 = 對所有節點取算術平均（= 體積平均，不是加總）。
%
%   S(p) 核：p̄=p/ℓ̂，S(:,k)=(p̄−Pc_k)/‖p̄−Pc_k‖³。
%   **frame-agnostic**：P 與 Pc 只要在同一框即可（no_fix bias = actuator 框）。
%   C/kappa 是奇異值純量 → 框無關。  (Was calc_range_metrics; renamed 2026-07-16.)
%
% 輸入：
%   P       Np×3  球內節點座標 [m]（與 Pc 同框）
%   Hhat_I  6×6   Ĥ_I = ^Bĝ_I·K̄ = gB·Khat [mT/A]
%   ell_m   純量  ℓ̂ [m]
%   Pc      3×6   等效磁荷正規化位置（與 P 同框）
% 輸出 struct m：.C_mean/.kappa_mean（mean）、.C_min/.kappa_worst（worst）、.C_std/.kappa_std、.Np
% -------------------------------------------------------------------------
    Np = size(P,1);
    C   = zeros(Np,1);   kap = zeros(Np,1);
    for i = 1:Np
        pbar = P(i,:).' / ell_m;                 % 3×1 無因次
        D    = pbar - Pc;                         % 3×6
        S    = D ./ (vecnorm(D).^3);             % 3×6 庫倫核
        sv   = svd(S * Hhat_I);                   % σ₁≥σ₂≥σ₃ [mT/A]
        C(i)   = prod(sv);                        % C = ∏σ_k = σ₁σ₂σ₃  [(mT/A)³]
        kap(i) = sv(3) / sv(1);                   % kappa = σ_min/σ_max = σ₃/σ₁
    end
    m = struct('C_mean',mean(C), 'kappa_mean',mean(kap), ...
               'C_min',min(C), 'kappa_worst',min(kap), ...
               'C_std',std(C), 'kappa_std',std(kap), 'Np',Np);
end
