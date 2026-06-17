function [P, Bstack, npts] = select_ball(D, R)
%SELECT_BALL  Select FEM nodes inside the ball ||p|| <= R and stack the field.
%   [P, Bstack, npts] = SELECT_BALL(D, R)
%     D : struct from load_coils_actuator (.Pa, .r2, .Ba)
%     R : ball radius [m]
%   Returns:
%     P      Np x 3        in-ball sample points (actuator frame)
%     Bstack 3Np x N_I     stacked [Bx;By;Bz] per simulation (column j)
%     npts   scalar        number of in-ball nodes (shared across the 6 coils)
    idx  = find(D.r2 < R^2);
    npts = numel(idx);
    P    = D.Pa(idx, :);
    N_I  = size(D.Ba, 3);
    Bstack = zeros(3*npts, N_I);
    for j = 1:N_I
        Bstack(:,j) = reshape(D.Ba(idx,:,j).', [], 1);
    end
end
