function m = calc_range_metrics(P, Hhat_I, ell_m, dhat)
% CALC_RANGE_METRICS  集總模型在控制範圍（R≤150µm 球）的性能量：體積 C 與均勻度 kappa。
% -------------------------------------------------------------------------
%   對球內每個真實 FEM 節點 p，算致動增益矩陣 ^BG_I(p)=S(p)·Ĥ_I（3×6，mT/A）的 SVD
%   （σ₁≥σ₂≥σ₃），依 singular_value.pdf 定義：
%     C(p) = ∏_{k=1}^{3} σ_k = σ₁σ₂σ₃   （flux-generating ellipsoid 體積；(mT/A)³）
%     kappa(p) = σ_min/σ_max = σ₃/σ₁     （condition number / isotropy；≤1，→1 等向）
%   代表值 = 對所有節點取算術平均（= 體積平均，不是加總）。
%
%   S(p) 核（同 charge_residual）：p̄=p/ℓ̂，S(:,k)=(p̄−d̂_k)/‖p̄−d̂_k‖³。
%   C/kappa 是奇異值純量、框無關 → P/dhat 用同一框即可。
%
% 輸入：
%   P       Np×3  球內節點座標 [m]（select_ball 已篩好）
%   Hhat_I  6×6   Ĥ_I = ^Bĝ_I·K̄ = gB·Khat [mT/A]
%   ell_m   純量  ℓ̂ [m]
%   dhat    3×6   電荷單位方向（與 P 同框）
% 輸出 struct m：
%   .C_mean        總代表體積 = mean_p ∏σ_k      [(mT/A)³]
%   .kappa_mean    總代表均勻度 = mean_p σ₃/σ₁   [--]
%   .C_min        最弱體積 = min_p ∏σ_k          [(mT/A)³]（worst-case 控制能力）
%   .kappa_worst  最不均 = min_p σ₃/σ₁           [--]（κ 越小越不均）
%   .C_std/.kappa_std  空間標準差
%   .Np           節點數
% -------------------------------------------------------------------------
    Np = size(P,1);
    C   = zeros(Np,1);   kap = zeros(Np,1);
    for i = 1:Np
        pbar = P(i,:).' / ell_m;                 % 3×1 無因次
        D    = pbar - dhat;                       % 3×6
        S    = D ./ (vecnorm(D).^3);             % 3×6 庫倫核
        sv   = svd(S * Hhat_I);                   % σ₁≥σ₂≥σ₃ [mT/A]
        C(i)   = prod(sv);                        % C = ∏σ_k = σ₁σ₂σ₃  [(mT/A)³]
        kap(i) = sv(3) / sv(1);                   % kappa = σ_min/σ_max = σ₃/σ₁
    end
    m = struct('C_mean',mean(C), 'kappa_mean',mean(kap), ...
               'C_min',min(C), 'kappa_worst',min(kap), ...
               'C_std',std(C), 'kappa_std',std(kap), 'Np',Np);
end
