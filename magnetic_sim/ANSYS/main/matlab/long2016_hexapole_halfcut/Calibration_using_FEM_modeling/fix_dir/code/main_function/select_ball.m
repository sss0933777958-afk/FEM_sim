function [coil, nmin] = select_ball(C, R)
%SELECT_BALL  Select FEM nodes inside the ball ||p|| <= R (WP frame) per coil.
%   [coil, nmin] = SELECT_BALL(C, R)
%     C : struct array from load_coils
%     R : ball radius [m]
%   Returns:
%     coil(k).p    Mx3   selected node positions [m]
%     coil(k).bfem 3M x1 stacked [Bx; By; Bz] (sign-negated FEM field)
%     coil(k).pj   paper pole index
%     nmin         minimum selected node count across the 6 coils
    coil = struct('p',{},'bfem',{},'pj',{});
    nmin = inf;
    for k = 1:numel(C)
        ii   = find(sum(C(k).P.^2, 2) < R^2);
        nmin = min(nmin, numel(ii));
        Bn   = C(k).Bn(ii,:);
        coil(k).p    = C(k).P(ii,:);
        coil(k).bfem = [Bn(:,1); Bn(:,2); Bn(:,3)];
        coil(k).pj   = C(k).pj;
    end
end
