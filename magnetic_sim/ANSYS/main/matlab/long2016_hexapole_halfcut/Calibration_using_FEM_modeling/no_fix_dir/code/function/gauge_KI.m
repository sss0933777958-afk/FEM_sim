function [KbarI, gB] = gauge_KI(ell, Pc, P, Bstack, F)
%GAUGE_KI  Steps 6-8: back out gB and K-bar_I from the profiled charge solution.
%   [KbarI, gB] = GAUGE_KI(ell, Pc, P, Bstack, F)
%   g_j  = (A'A)\(A'b_j)  ->  G = [g_1 ... g_NI]  ->  H = G F'(FF')^-1.
%   Gauge with k-bar_I(1,1)=5/6:  gB = (6/5) h11,  K-bar_I = (5/(6 g11)) H.
%   (no_fix_l.pdf steps 6-8.  g11 = G(1,1) = h11 in this permutation calibration.)
    A   = build_A(ell, Pc, P);
    M   = A.' * A;
    C   = A.' * Bstack;
    G   = M \ C;
    H   = G * F.' / (F * F.');
    g11 = G(1,1);
    KbarI = (5/(6*g11)) * H;
    gB    = (6/5) * H(1,1);
end
