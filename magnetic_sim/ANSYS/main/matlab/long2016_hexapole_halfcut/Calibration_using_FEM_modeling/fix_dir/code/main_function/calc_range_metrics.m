function m = calc_range_metrics(P, Hhat_I, ell_m, dhat)
% CALC_RANGE_METRICS  集總模型在控制範圍（R≤150µm 球）的代表增益 A 與均勻度 iso。
% -------------------------------------------------------------------------
%   對球內每個真實 FEM 節點 p，算電流→場轉移 T(p)=S_i(p)·Ĥ_I（3×6，mT/A）的 SVD：
%     P(p)   = ‖T‖_F = √(σ₁²+σ₂²+σ₃²)   （gain，符號 P）
%     iso(p) = σ_max/σ_min = σ₁/σ₃        （條件數；≥1，=1 等向）
%   代表值 = 對所有節點取算術平均（= 體積平均，不是加總；加總會 ∝ 節點數）。
%
%   S_i(p) 核（同 charge_residual）：p̄=p/ℓ̂，S(:,k)=(p̄−d̂_k)/‖p̄−d̂_k‖³。
%   gain/iso 是奇異值純量、框無關 → P/dhat 用 measure 框即可。
%
% 輸入：
%   P       Np×3  球內節點座標 [m]（WP/measure 框；select_ball 已篩好）
%   Hhat_I  6×6   Ĥ_I = ^Bĝ_I·K̄ = gB·Khat [mT/A]
%   ell_m   純量  ℓ̂ [m]
%   dhat    3×6   極尖單位方向（measure 框）
% 輸出 struct m：
%   .sigma_tot  總代表增益 = mean_p ‖T‖_F   [mT/A]（gain 總量符號 σ_tot）
%   .iso_tot    總代表均勻度 = mean_p σ₁/σ₃  [--]（iso_tot）
%   .sigma_min  最弱增益 = min_p ‖T‖_F        [mT/A]（worst-case 控制能力）
%   .iso_worst  最不均 = max_p σ₁/σ₃          [--]
%   .sigma_std/.iso_std  空間標準差
%   .Np         節點數
% -------------------------------------------------------------------------
    Np = size(P,1);
    g  = zeros(Np,1);  is = zeros(Np,1);
    for i = 1:Np
        pbar = P(i,:).' / ell_m;                 % 3×1 無因次
        D    = pbar - dhat;                       % 3×6
        S    = D ./ (vecnorm(D).^3);             % 3×6 庫倫核
        sv   = svd(S * Hhat_I);                   % σ₁≥σ₂≥σ₃ [mT/A]
        g(i)  = norm(sv);                         % P = ‖T‖_F = √Σσ²
        is(i) = sv(1) / sv(3);                    % σ_max/σ_min
    end
    m = struct('sigma_tot',mean(g), 'iso_tot',mean(is), ...
               'sigma_min',min(g), 'iso_worst',max(is), ...
               'sigma_std',std(g), 'iso_std',std(is), 'Np',Np);
end
