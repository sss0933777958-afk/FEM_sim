function f = base_model(L, H_hat, u, c)
%BASE_MODEL  Force model   f_ijk = c * (H_hat*u_j)' * L_k * (H_hat*u_j).
%
%   f = BASE_MODEL(L, H_hat, u)
%   f = BASE_MODEL(L, H_hat, u, c)
%     L     : Np x 6 x 6 x 3   spatial-distribution tensor, from build_L
%     H_hat : 6 x 6            transfer matrix acting on the excitation
%     u     : 6 x N_I          excitation matrix, column j is u_j
%     c     : scalar coefficient, default 1
%     f     : 3 x Np x N_I     force [pN];  f(k,i,j) = component k, point i,
%                              excitation j
%
%   [RENAMED 2026-09-02] from current_base_model. All four arguments are
%   generic symbols: the function does not know, and does not need to know,
%   which base it serves. One expression covers both bases and both forms in
%   which the model is written.
%
%   EXPANDED FORM -- used to GENERATE the measurement f_m. The scalar prefactor
%   stays outside and H_hat holds the gauged, dimensionless matrix Mbar:
%
%       current   H_hat = K_I bar,  u = F  [A],   c = F g_I  [pN/A^2]
%       voltage   H_hat = D bar,    u = V  [mV],  c = F g_V  [pN/mV^2]
%
%   MERGED FORM -- used to EVALUATE the model inside the fit. The prefactor is
%   folded into the matrix, so the coefficient is 1:
%
%       current   H_hat = F H_I hat,  u = F,  c = 1
%       voltage   H_hat = F H_V hat,  u = V,  c = 1
%
%   The two agree because  F H hat = sqrt(F g) * Mbar : the square root sits on
%   both sides of L and reproduces c exactly. Routing the generated measurement
%   and the fitted model through this one function guarantees they cannot drift
%   apart in definition.
%
%   The excitation u comes from the flux calibration, not from a convention set
%   here: Fmap for a current calibration, V (the sensor voltages produced by the
%   project's voltage read-out) for a voltage one.
%
%   IDENTIFIABILITY.  f is invariant under (H_hat*u_j) -> -(H_hat*u_j), so with
%   a one-at-a-time excitation set the sign of each column of H_hat is not fixed
%   by f alone. That ambiguity is resolved later, in solve_force, using the
%   physical requirement diag(Mbar) > 0.
%
%   See also BUILD_L, SOLVE_FORCE, LOAD_FLUX_CALIB, FITTING_FORCE.

    if nargin < 4 || isempty(c), c = 1; end
    validateattributes(H_hat, {'numeric'}, {'size',[6 6],'real','finite'}, mfilename, 'H_hat');
    validateattributes(c, {'numeric'}, {'scalar','real','finite'}, mfilename, 'c');
    validateattributes(u, {'numeric'}, {'2d','nrows',6,'real','finite'}, mfilename, 'u');
    assert(size(L,2) == 6 && size(L,3) == 6 && size(L,4) == 3, ...
           'base_model:LSize', 'L must be Np x 6 x 6 x 3');

    Np = size(L, 1);
    NI = size(u, 2);

    W  = H_hat * u;                 % 6 x N_I, column j is w_j = H_hat*u_j
    Lr = reshape(L, Np, 36, 3);     % collapse the (n,m) charge pair into one index

    f = zeros(3, Np, NI);
    for j = 1:NI
        % outer product w*w' flattened in the same column-major (n,m) order as Lr
        q = reshape(W(:,j) * W(:,j).', 36, 1);
        for k = 1:3
            f(k,:,j) = c * (Lr(:,:,k) * q).';
        end
    end
end
