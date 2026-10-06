function fit = sp_fit(P, b, l0, use_bias)
%SP_FIT  Single-pole point-charge fit in the tip frame.
%   fit = SP_FIT(P, b, l0, use_bias)
%     P        Np x 3 sample positions [m], tip frame (origin = tip apex, pole axis +x)
%     b        Np x 3 field [mT]
%     l0       initial l_hat [m]
%     use_bias false: charge on the axis; true: transverse offset e_y, e_z also fitted
%   Model (normalised by l_hat):  b_i = G * S_i,  S_i = d_i/|d_i|^3,  d_i = p_i/l_hat - (1, e_y, e_z)
%     -> charge at l_hat*(1, e_y, e_z); l_hat is the distance along x from the tip.
%     e_x is pinned to 0 (only l_hat*(1+e_x) is identifiable). G [mT] is solved in closed form.
%   Returns .l_hat [m], .e = [e_x; e_y; e_z], .G [mT], .J (residual sum of squares [mT^2]),
%           .predict(Pq) -> Nq x 3 [mT], .exitflag

    Pu = P * 1e6;                                             % fit in um for scaling
    x0 = l0 * 1e6;  lb = 1e-3;  ub = 5e4;                     % l_hat > 0: charge inside the iron
    if use_bias, x0 = [x0; 0; 0];  lb = [lb; -Inf; -Inf];  ub = [ub; Inf; Inf]; end
    opts = optimoptions('lsqnonlin', 'Display','off', 'MaxFunctionEvaluations',1e5, ...
                        'MaxIterations',4e3, 'FunctionTolerance',1e-20, 'StepTolerance',1e-12);
    [xf, ~, ~, ef] = lsqnonlin(@(x) resid(x, Pu, b), x0, lb, ub, opts);
    [r, G] = resid(xf, Pu, b);

    fit.l_hat = xf(1) * 1e-6;
    fit.e     = [0; 0; 0];
    if use_bias, fit.e(2:3) = xf(2:3); end
    fit.G     = G;
    fit.J     = sum(r.^2);
    fit.exitflag = ef;
    c = [1, fit.e(2), fit.e(3)];
    fit.predict = @(Pq) G * shape(Pq / fit.l_hat, c);
end

function [r, G] = resid(x, Pu, b)
    c = [1, 0, 0];
    if numel(x) == 3, c(2:3) = x(2:3); end
    S = shape(Pu / x(1), c);
    G = sum(S(:) .* b(:)) / sum(S(:).^2);
    r = G * S(:) - b(:);
end

function S = shape(Pn, c)
    d = Pn - c;
    S = d ./ vecnorm(d, 2, 2).^3;
end
