function [KI_bar, gI_hat, G, rm] = solve_current(l_hat, e, Pc_base, P, Bstack, F, K22_set)
%SOLVE_CURRENT  電流側解參數：優化器解完 → 重建 S/Pc_bar → G → K̄_I / ĝ_I / 𝒞κ。
%   [KI_bar, gI_hat, G, rm] = SOLVE_CURRENT(l_hat, e, Pc_base, P, Bstack, F, K22_set)
%     l_hat : 有效長度 [m]；e：17×1 偏移；Pc_base：3×6 理想電荷格
%     P     : Np×3 [m]；Bstack：3Np×6 all-source [mT]；F：6×N_I coil→pole 連接矩陣
%     K22_set（可選）：K̄(2,2) 物理約束值（[]=自由擬合）。見下方 [ADDED] 段。
%   Pc=make_Pc(e) → S=build_S(l,Pc) → G=(SᵀS)\(SᵀBstack)（6×N_I profiled 電荷）。
%   H_I=G·Fᵀ(FFᵀ)⁻¹；gI_hat=(6/5)H_I(1,1)[mT/A]；KI_bar=(5/(6·G(1,1)))·H_I（gauge K̄(1,1)=5/6）；
%   rm=𝒞/κ（逐節點 svd(S·(gI_hat·KI_bar))）。
    if nargin < 7, K22_set = []; end
    Pc = make_Pc(e, Pc_base);
    S  = build_S(l_hat, Pc, P);
    G  = (S.'*S) \ (S.'*Bstack);                     % 6×N_I

    % [ADDED 2026-08-10] K̄(2,2) 物理約束（規則 calibration-transfer-matrix-output.md 附則二）：
    %   「上磁極全錐應強於半切下磁極」，但該方向對資料近乎退化 → 在資料無法分辨的區間內設定它。
    %   **只改 G(2,2)** → (6/5)·K22_set·G(1,1)；ℓ̂ / e / 其餘 35 個 G 元素不動。
    %   因 gauge 除的是 G(1,1)（未動），K̄ 其餘元素與 ĝ_I 完全不變，只有 (2,2) 改變。
    %   ⚠ voltage 端必須吃**同一顆 G**（solve_voltage 帶同一個 K22_set）。
    %   ⚠ K̄(2,2)=(5/6)·G(2,2)/G(1,1) 僅在 F=identity 時成立（Maxwell map = identity）→ 此處強制檢查。
    if ~isempty(K22_set)
        assert(isequal(F, eye(size(F,1))), ...
               'solve_current: K22_set 只支援 F=identity 的 map（本分支 Maxwell 為 identity）');
        rm_K22_free = (5/6) * G(2,2) / G(1,1);       % 約束前的自由擬合值（供追溯）
        G(2,2)      = (6/5) * K22_set * G(1,1);
    end

    H_I    = G * F.' / (F * F.');                    % Ĥ_I（un-gauged）
    gI_hat = (6/5) * H_I(1,1);                       % ĝ_I [mT/A]
    KI_bar = (5/(6*G(1,1))) * H_I;                   % K̄_I，gauge K̄(1,1)=5/6
    rm     = control_metrics(P, gI_hat*KI_bar, l_hat, Pc);   % 𝒞/κ（物理 Ĥ_I=ĝ·K̄）
    resid    = S*G - Bstack;                          % 擬合殘差 ε [mT]（模型 S·G vs FEM b）
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
        S  = (d ./ r3).';                     % 3×6
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
function S = build_S(l_hat, Pc, P)
    Np   = size(P, 1);
    pbar = P / l_hat;
    S    = zeros(3*Np, 6);
    for k = 1:6
        d  = pbar - Pc(:,k).';
        r3 = sum(d.^2, 2).^1.5;
        S(:,k) = reshape((d ./ r3).', 3*Np, 1);
    end
end
