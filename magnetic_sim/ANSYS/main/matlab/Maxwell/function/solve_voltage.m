function [D_bar, gV_hat, G, rm] = solve_voltage(l_hat, e, Pc_base, P, Bstack, V, K22_set)
%SOLVE_VOLTAGE  電壓側解參數：優化器解完 → 重建 S/Pc_bar → G → D̄ / ĝ_V / 𝒞κ。
%   [D_bar, gV_hat, G, rm] = SOLVE_VOLTAGE(l_hat, e, Pc_base, P, Bstack, V, K22_set)
%     l_hat : 有效長度 [m]；e：17×1 偏移；Pc_base：3×6 理想電荷格
%     P     : Np×3 [m]；Bstack：3Np×6 all-source [mT]；V：6×6 sensor 電壓矩陣 [mV]（all-source）
%     K22_set（可選）：K̄(2,2) 物理約束值（[]=自由擬合）。見下方 [ADDED] 段。
%   Pc=make_Pc(e) → A=build_S(l,Pc) → G=(AᵀA)\(AᵀBstack)（6×N_I profiled 電荷）。
%   H_V=(G·Vᵀ)/(V·Vᵀ)；gV_hat=(6/5)H_V(1,1)[mT/mV]；D_bar=(5/(6·H_V(1,1)))·H_V（gauge D̄(1,1)=5/6）；
%   rm=𝒞/κ（逐節點 svd(S·(gV_hat·D_bar))）。
    if nargin < 7, K22_set = []; end
    Pc = make_Pc(e, Pc_base);
    A  = build_S(l_hat, Pc, P);
    G  = (A.'*A) \ (A.'*Bstack);                     % 6×N_I

    % [ADDED 2026-08-10] K̄(2,2) 物理約束 —— **與 solve_current 逐字相同的 G 修改**。
    %   規則 calibration-transfer-matrix-output.md 附則二寫死：G 一改，voltage 必須吃同一顆 G，
    %   否則同一顆物理電荷會有兩個版本、K̄_I 與 D̄ 不自洽。
    %   （兩側的 G 本來就只吃 WP 場、與 sensor 無關，逐位相同；故這裡重算的 G 與 current 端一致。）
    if ~isempty(K22_set)
        rm_K22_free = (5/6) * G(2,2) / G(1,1);       % 約束前的自由擬合值（供追溯）
        G(2,2)      = (6/5) * K22_set * G(1,1);
    end

    H_V    = (G * V.') / (V * V.');                  % Ĥ_V（un-gauged）
    gV_hat = (6/5) * H_V(1,1);                       % ĝ_V [mT/mV]
    D_bar  = (5/(6*H_V(1,1))) * H_V;                 % D̄，gauge D̄(1,1)=5/6
    rm     = control_metrics(P, gV_hat*D_bar, l_hat, Pc);   % 𝒞/κ（物理 Ĥ_V=ĝ·D̄）
    resid    = A*G - Bstack;                          % [ADDED] 擬合殘差 ε [mT]（與 solve_current 同定義）
    rm.RMSPE = sqrt(sum(resid(:).^2) / sum(Bstack(:).^2)) * 100;   % RMSPE [%] = sqrt(Σε²/Σb²)·100
    % [ADDED 2026-08-15 使用者拍板] PDF 改印 NMAE（L1 版）；RMSPE 仍算、只存 .mat 供追溯。
    rm.NMAE  = sum(abs(resid(:))) / sum(abs(Bstack(:))) * 100;     % NMAE [%] = Σ|ε|/Σ|b|·100
    if ~isempty(K22_set)                              % [ADDED] 約束值與自由值一併回報，供 .mat 追溯
        rm.K22_set = K22_set;   rm.K22_free = rm_K22_free;
    end
end

% ---- 控制範圍指標：逐節點 svd(S·Ĥ) → 𝒞=∏σ、κ=σ₃/σ₁ ----
function rm = control_metrics(P, Hhat, l_hat, Pc)
    Np  = size(P, 1);
    C   = zeros(Np, 1);
    kap = zeros(Np, 1);
    for i = 1:Np
        pbar = P(i,:) / l_hat;
        d  = pbar - Pc.';
        r3 = sum(d.^2, 2).^1.5;
        S  = (d ./ r3).';
        sv = svd(S * Hhat);
        C(i)   = prod(sv);
        kap(i) = sv(3) / sv(1);
    end
    rm.C_mean = mean(C);   rm.kappa_mean  = mean(kap);
    rm.C_min  = min(C);    rm.kappa_worst = min(kap);
    rm.C_std  = std(C);    rm.kappa_std   = std(kap);
    rm.Np     = Np;
end

% ---- 電荷格 Pc_bar = Pc_base + E(e)（含 e6z 約束）----
function Pc = make_Pc(e17, Pc_base)
    E = zeros(3, 6);
    E(:,1) = e17(1:3);     E(:,2) = e17(4:6);
    E(:,3) = e17(7:9);     E(:,4) = e17(10:12);
    E(:,5) = e17(13:15);
    E(1,6) = e17(16);      E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);
    Pc = Pc_base + E;
end

% ---- 每點無因次 Coulomb kernel（3Np×6）----
function A = build_S(l_hat, Pc, P)
    Np   = size(P, 1);
    pbar = P / l_hat;
    A    = zeros(3*Np, 6);
    for k = 1:6
        d  = pbar - Pc(:,k).';
        r3 = sum(d.^2, 2).^1.5;
        A(:,k) = reshape((d ./ r3).', 3*Np, 1);
    end
end
