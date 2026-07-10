function J = sensor_residual(P, B, Vmat, exc_sign, ell_hat, Pc, d, N_I)
% SENSOR_RESIDUAL  算 sensor 模型的殘差 ε 與成本 J（= solve_d 最小化的目標；無 g_H）。
% -------------------------------------------------------------------------
% 用途：
%   解出 d 後，算 sensor 模型場 b_model = S_i·V_j·d（無 g_H，與 no-gain d 一致）與
%   all-source FEM 物理場 b_FEM 的殘差成本：
%       殘差 ε_ij = b_model − b_FEM
%       成本 J    = Σ_ij ‖ε_ij‖²            （最小二乘目標值，單位 T²）
%   （依使用者要求：只回 cost J，不算相對 RMSE。）
%
% 符號對齊：
%   b_model = S_i diag(V_j) d                （Vmat 已是 all-source 翻號版；無 g_H）
%   b_FEM   = exc_sign(j)·(−B(i,:,j))        （−B = 物理場；×exc_sign = all-source）
%
% 輸入：
%   P        Np×3 工作點 [m]
%   B        Np×3×N_I 記錄的場（負號版 −B^FEM）
%   Vmat     6×6 all-source sensor 電壓
%   exc_sign 1×6 all-source 翻號
%   ell_hat  特徵長度 ℓ̂ [m]
%   Pc       3×6 電荷位置（極尖單位方向）
%   d        6×1 每極常數（solve_d 解出，no-gain）
%   N_I      模擬數（6）
% 輸出：
%   J        成本（殘差平方和 Σ‖ε‖²）[T²]
%
%   (calls build_S)
% -------------------------------------------------------------------------
    Np = size(P,1);                                          % 工作點數
    J  = 0;                                                  % 成本累加器 J = Σ_j (Σ_i ‖ε_ij‖²)
    for j = 1:N_I                                            % === 外層：第 j 次模擬 ===
        Vj = diag(Vmat(:,j));                               % V_j = diag(該模擬的 6 個 all-source sensor 電壓)
        for i = 1:Np                                        % === 內層：該模擬累加所有工作點 i ===
            Si = build_S(P(i,:), ell_hat, Pc);             % 該點 3×6 空間核（在 ℓ̂ 下）
            bm = Si * Vj * d;                              % sensor 模型場 b_model = S_i V_j d（無 g_H、all-source）
            bf = exc_sign(j) * (-squeeze(B(i,:,j)).');     % all-source FEM 物理場 b_FEM（−B=物理，×exc_sign=all-source）
            J  = J + sum((bm-bf).^2);                      % 累加殘差平方 ‖ε_ij‖²
        end
    end
end
