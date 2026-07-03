function m = calc_range_metrics(P, Hhat_I, ell_m, Pc)
% CALC_RANGE_METRICS  集總模型在控制範圍（R≤150µm 球）的總代表增益 σ_tot 與均勻度 iso_tot。
% -------------------------------------------------------------------------
%   對球內每個真實 FEM 節點 p，算電流→場轉移 T(p)=S_i(p)·Ĥ_I（3×6，mT/A）的 SVD：
%     gain(p) = ‖T‖_F = √(σ₁²+σ₂²+σ₃²)   [mT/A]
%     iso(p)  = σ_max/σ_min = σ₁/σ₃        [--]（條件數；≥1，=1 等向）
%   代表值 = 對所有節點取算術平均（= 體積平均，不是加總）。
%
%   S_i(p) 核：p̄=p/ℓ̂，S(:,k)=(p̄−Pc_k)/‖p̄−Pc_k‖³。
%   **frame-agnostic**：P 與 Pc 只要在同一框即可（no_fix bias = actuator 框：P=select_ball 節點、
%   Pc=make_Pc 電荷；fix = measure 框：Pc=dhat）。gain/iso 是奇異值純量 → 框無關。
%
% 輸入：
%   P       Np×3  球內節點座標 [m]（與 Pc 同框）
%   Hhat_I  6×6   Ĥ_I = ^Bĝ_I·K̄ = gB·Khat [mT/A]
%   ell_m   純量  ℓ̂ [m]
%   Pc      3×6   等效磁荷正規化位置（與 P 同框；bias=make_Pc、fix=dhat）
% 輸出 struct m：.sigma_tot/.iso_tot（mean）、.sigma_min/.iso_worst、.sigma_std/.iso_std、.Np
% -------------------------------------------------------------------------
    Np = size(P,1);
    g  = zeros(Np,1);  is = zeros(Np,1);
    for i = 1:Np
        pbar = P(i,:).' / ell_m;                 % 3×1 無因次
        D    = pbar - Pc;                         % 3×6
        S    = D ./ (vecnorm(D).^3);             % 3×6 庫倫核
        sv   = svd(S * Hhat_I);                   % σ₁≥σ₂≥σ₃ [mT/A]
        g(i)  = norm(sv);                         % gain = ‖T‖_F = √Σσ²
        is(i) = sv(1) / sv(3);                    % σ_max/σ_min
    end
    m = struct('sigma_tot',mean(g), 'iso_tot',mean(is), ...
               'sigma_min',min(g), 'iso_worst',max(is), ...
               'sigma_std',std(g), 'iso_std',std(is), 'Np',Np);
end
