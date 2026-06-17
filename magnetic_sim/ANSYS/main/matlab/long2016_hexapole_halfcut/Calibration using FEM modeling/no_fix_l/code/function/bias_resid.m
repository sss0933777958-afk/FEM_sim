function r = bias_resid(x, P, Bstack, Pc_base)
%BIAS_RESID  lsqnonlin residual for the 18-param bias model (g_j profiled out).
%   r = BIAS_RESID(x, P, Bstack, Pc_base)
%     x       : [ell_mm; e17(17x1)]
%     P       : Np x 3 sample points (actuator frame)
%     Bstack  : 3Np x N_I stacked FEM field (column j = simulation j)
%     Pc_base : 3 x 6 ideal charge lattice
%   For each simulation the best charge vector g_j is solved by least squares
%   (profiled), and its residual A*g_j - b_j is stacked.  sum(r.^2) = cost J.
    A = build_A(x(1)*1e-3, make_Pc(x(2:18), Pc_base), P);   % 3Np x 6
    M = A.' * A;                                            % 6 x 6
    r = [];
    for j = 1:size(Bstack,2)
        gj = M \ (A.' * Bstack(:,j));
        r  = [r; A*gj - Bstack(:,j)]; %#ok<AGROW>
    end
end
