function [KI_bar, gI_hat, G, rm] = solve_current(l_hat, e, Pc_base, P, Bstack, F, K22_set, Peval, Beval)
%SOLVE_CURRENT  電流側解參數：優化器解完 → 重建 S/Pc_bar → G → K̄_I / ĝ_I / 𝒞κ。
%   [KI_bar, gI_hat, G, rm] = SOLVE_CURRENT(l_hat, e, Pc_base, P, Bstack, F, K22_set, Peval, Beval)
%     l_hat : 有效長度 [m]；e：17×1 偏移；Pc_base：3×6 理想電荷格
%     P     : Np×3 [m]；Bstack：3Np×6 all-source [mT]；F：6×N_I coil→pole 連接矩陣
%     K22_set（可選）：K̄(2,2) 物理約束值（[]=自由擬合）。見下方 [ADDED] 段。
%     Peval / Beval（可選）：**評估**用的點雲與場。見下方 [ADDED 2026-08-26]。
%   Pc=make_Pc(e) → S=build_S(l,Pc) → G=(SᵀS)\(SᵀBstack)（6×N_I profiled 電荷）。
%   H_I=G·Fᵀ(FFᵀ)⁻¹；gI_hat=(6/5)|H_I(1,1)|[mT/A]；KI_bar=(5/(6·G(1,1)))·H_I（gauge K̄(1,1)=5/6）；
%   rm=𝒞/κ（逐節點 svd(S·(gI_hat·KI_bar))）。
%
%   [ADDED 2026-08-26 使用者拍板] 誤差指標改成 **out-of-sample**：
%     給了 (Peval, Beval) 時，NMAE / RMSPE / NMAE_L1 在**那組點**上計算 —— 用的仍是
%     由 (P, Bstack) 擬合出來的同一顆 G，只是把 kernel 重建在 Peval 上。main.m 傳的是
%     「R 內全部 .fld 格點」，所以 PDF 上的 NMAE = 「用 N 個取樣點校正出來的模型，
%     拿去擬合全格點的場」。
%     不給就退回舊行為（in-sample，於 P 上計算），既有呼叫端逐字不變。
%     in-sample 版一律另存 rm.NMAE_in / rm.RMSPE_in / rm.NMAE_L1_in 供追溯。
%     ⚠ 為什麼要改：in-sample 對點數有系統性樂觀 —— 35 點的 in-sample NMAE 是 2.60%，
%       比全格點擬合的 2.81% 還「好」，純粹因為點少比較好配。同一顆模型拿去配全格點
%       其實是 2.81%，跟全格點模型一樣。
%     ⚠ 𝒞 / κ **仍在 P（擬合點）上計算**，未一併改；要改需另議（會動到所有既有紀錄）。
    if nargin < 7, K22_set = []; end
    if nargin < 8, Peval   = []; end
    if nargin < 9, Beval   = []; end
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
    % [MODIFIED 2026-08-17 使用者拍板] ĝ_I 取 |H_I(1,1)|，增益一律回報**大小**。
    %   起因：ĝ_I 由**單一元素** H_I(1,1)（=P1 自激發電荷）定義，是所有輸出量裡最脆弱的。
    %   R<=40um 取樣時六顆電荷近簡併，G 的下極對角整層翻負（P1/P3/P6 −33/−26/−34、
    %   上極 +45/+48/+45）→ ĝ_I 被 (1,1) 的負號直接繼承，畫出 −40 mT/A 的假谷。
    %   ⚠ 對 R>=60um 的既有結果**完全無影響**（G(1,1) 恆正，abs 是恆等）。
    %   ⚠ 副作用（僅在 G(1,1)<0 的退化區）：K̄_I 的 gauge 仍除**帶號**的 G(1,1)，
    %     故此時 ĝ_I·K̄_I = −H_I，規則 calibration-transfer-matrix-output.md 的
    %     恆等式 ᴮĤ_I = ĝ_I·K̄_I = H_I 不再成立。該區間本來就判定為不可用，不另處理。
    %   ⚠ 𝒞 / κ 不受影響：control_metrics 取 svd，整體變號不改奇異值。
    gI_hat = (6/5) * abs(H_I(1,1));                  % ĝ_I [mT/A]
    KI_bar = (5/(6*G(1,1))) * H_I;                   % K̄_I，gauge K̄(1,1)=5/6
    rm     = control_metrics(P, gI_hat*KI_bar, l_hat, Pc);   % 𝒞/κ（物理 Ĥ_I=ĝ·K̄）
    % [MODIFIED 2026-08-17 使用者拍板] NMAE 為**向量範數**版：
    %   NMAE = [ Σ_j Σ_i ‖b_ij − S_i·ᴮĝ_I·K̄_I·F_j‖ / N_p ] / b̄ · 100，
    %   b̄ = Σ_j Σ_i ‖b_ij‖ / N_p（同結構 → N_p 對消，比值無因次）。
    %   ⚠ 與更舊版（逐**分量**取絕對值 Σ|ε|/Σ|b|）不同：這裡先把每個「點×激發」的 3 維殘差
    %     取歐氏長度 ‖ε_ij‖ 再相加。同一份殘差下新值必 ≥ 舊值（‖·‖₂ ≥ 各分量平均）。
    %   ⚠ S_i·ĝ·K̄·F_j ≡ S_i·G_j（因 Ĥ_I·F = G），故實作用 S*G；僅在退化區 h11<0
    %     （ĝ 已取 abs）兩者差一個整體負號，該區間本就判定不可用。
    % [MODIFIED 2026-08-26] 同一組公式現在算兩次：in-sample（P）與 evaluation（Peval）。
    [rm.RMSPE_in, rm.NMAE_in, rm.NMAE_L1_in] = fitmetrics(S, G, Bstack);
    if isempty(Peval)
        rm.RMSPE = rm.RMSPE_in;   rm.NMAE = rm.NMAE_in;   rm.NMAE_L1 = rm.NMAE_L1_in;
        rm.NMAE_on = 'fit points';   rm.Np_eval = size(P,1);
    else
        Se = build_S(l_hat, Pc, Peval);               % 同一顆 G，kernel 重建在評估點上
        [rm.RMSPE, rm.NMAE, rm.NMAE_L1] = fitmetrics(Se, G, Beval);
        rm.NMAE_on = 'eval points';  rm.Np_eval = size(Peval,1);
    end
    if ~isempty(K22_set)                              % [ADDED] 約束值與自由值一併回報，供 .mat 追溯
        rm.K22_set = K22_set;   rm.K22_free = rm_K22_free;
    end
end

% ---- 擬合誤差指標：給定 kernel / 電荷 / 場，算 RMSPE、NMAE、舊 L1 版 ----
%   reshape(·,3,[]) 依 Bstack 的堆疊順序 [bx;by;bz] 逐點成行，欄 = (點, 激發)。
function [rmspe, nmae, nmae_l1] = fitmetrics(S, G, B)
    resid   = S*G - B;                                       % ε [mT]
    rmspe   = sqrt(sum(resid(:).^2) / sum(B(:).^2)) * 100;
    e_ij    = sqrt(sum(reshape(resid, 3, []).^2, 1));         % 每點每激發的 ‖ε‖
    b_ij    = sqrt(sum(reshape(B,     3, []).^2, 1));         % 同上的 ‖b‖
    nmae    = sum(e_ij) / sum(b_ij) * 100;
    nmae_l1 = sum(abs(resid(:))) / sum(abs(B(:))) * 100;      % 舊 L1 分量版，留存追溯
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
