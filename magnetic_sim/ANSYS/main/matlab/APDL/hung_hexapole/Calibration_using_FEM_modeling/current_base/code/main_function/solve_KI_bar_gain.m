function [K_bar, ghat_I_B, G] = solve_KI_bar_gain(ell, Pc, P, Bstack, F)
%SOLVE_KI_BAR_GAIN  Back out K-bar_I and the current-base gain ^Bg_I from the fit.
%   [K_bar, ghat_I_B, G] = SOLVE_KI_BAR_GAIN(ell, Pc, P, Bstack, F)
%   g_j = (A'A)\(A'b_j) -> G = [g_1 ... g_NI] (= D^v) -> H_I = G F'(FF')^-1 (= Hhat_I).
%   Gauge with K_bar(1,1)=5/6:  ^Bg_I = (6/5) h11,  K_bar = (5/(6 g11)) H_I.
%   (paper steps 6-8.  Was gauge_KI; renamed 2026-07-16.  Also returns G [mT].)
    A   = build_S_matrix(ell, Pc, P);
    M   = A.' * A;
    C   = A.' * Bstack;
    G   = M \ C;                       % G = D^v (profiled per-excitation charges)
    H_I = G * F.' / (F * F.');         % Hhat_I = G F'(FF')^-1 (current-side, un-gauged)
    g11 = G(1,1);
    K_bar    = (5/(6*g11)) * H_I;      % K-bar (gauge K_bar(1,1)=5/6)
    ghat_I_B = (6/5) * H_I(1,1);       % ^Bg_I (current-side gain, mT/A)
end
