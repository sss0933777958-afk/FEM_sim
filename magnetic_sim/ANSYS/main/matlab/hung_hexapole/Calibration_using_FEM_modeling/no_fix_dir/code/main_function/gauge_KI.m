function [K_bar, ghat_I_B] = gauge_KI(ell, Pc, P, Bstack, F)
%GAUGE_KI  Paper steps 6-8: back out ^Bg_I and K-bar from the profiled charge solution.
%   [K_bar, ghat_I_B] = GAUGE_KI(ell, Pc, P, Bstack, F)
%   g_j = (A'A)\(A'b_j) -> G = [g_1 ... g_NI] (= D^v) -> H_I = G F'(FF')^-1 (= Hhat_I).
%   Gauge with K_bar(1,1)=5/6:  ^Bg_I = (6/5) h11,  K_bar = (5/(6 g11)) H_I.
%   (paper steps 6-8.  g11 = G(1,1) = H_I(1,1) = h11 in this permutation calibration.)
    A   = build_A(ell, Pc, P);
    M   = A.' * A;
    C   = A.' * Bstack;
    G   = M \ C;                       % G = D^v (profiled per-excitation charges)
    H_I = G * F.' / (F * F.');         % Hhat_I = G F'(FF')^-1 (current-side, un-gauged)
    g11 = G(1,1);
    K_bar    = (5/(6*g11)) * H_I;      % K-bar (gauge K_bar(1,1)=5/6)
    ghat_I_B = (6/5) * H_I(1,1);       % ^Bg_I (current-side gain, T/A)
end
