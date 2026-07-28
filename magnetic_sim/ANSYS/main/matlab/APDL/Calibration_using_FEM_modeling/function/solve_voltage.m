function [D_bar, gV_hat, G, rm] = solve_voltage(l_hat, e, Pc_base, P, Bstack, V)
%SOLVE_VOLTAGE  電壓側解參數：優化器解完 → 重建 S/Pc_bar → G → D̄ / ĝ_V / 𝒞κ。
%   [D_bar, gV_hat, G, rm] = SOLVE_VOLTAGE(l_hat, e, Pc_base, P, Bstack, V)
%     l_hat : 有效長度 [m]；e：17×1 偏移；Pc_base：3×6 理想電荷格
%     P     : Np×3 [m]；Bstack：3Np×6 all-source [mT]；V：6×6 sensor 電壓矩陣 [mV]（all-source）
%   Pc=make_Pc(e) → A=build_S(l,Pc) → G=(AᵀA)\(AᵀBstack)（6×N_I profiled 電荷）。
%   H_V=(G·Vᵀ)/(V·Vᵀ)；gV_hat=(6/5)H_V(1,1)[mT/mV]；D_bar=(5/(6·H_V(1,1)))·H_V（gauge D̄(1,1)=5/6）；
%   rm=𝒞/κ（逐節點 svd(S·(gV_hat·D_bar))）。
    Pc = make_Pc(e, Pc_base);
    A  = build_S(l_hat, Pc, P);
    G  = (A.'*A) \ (A.'*Bstack);                     % 6×N_I

    H_V    = (G * V.') / (V * V.');                  % Ĥ_V（un-gauged）
    gV_hat = (6/5) * H_V(1,1);                       % ĝ_V [mT/mV]
    D_bar  = (5/(6*H_V(1,1))) * H_V;                 % D̄，gauge D̄(1,1)=5/6
    rm     = control_metrics(P, gV_hat*D_bar, l_hat, Pc);   % 𝒞/κ（物理 Ĥ_V=ĝ·D̄）
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
