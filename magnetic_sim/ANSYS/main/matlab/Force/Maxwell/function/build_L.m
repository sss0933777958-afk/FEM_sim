function L = build_L(P, l_hat, e, Pc_base)
%BUILD_L  Position-dependent 6x6x3 dimensionless spatial-distribution tensor.
%
%   L = BUILD_L(P, l_hat, e, Pc_base)
%     P       : Np x 3  sample positions [um], actuator frame (same as the flux fit)
%     l_hat   : effective length [um]  (P and l_hat only ever appear as P/l_hat,
%               so any consistent length unit works; this package uses um)
%     e       : 17x1    dimensionless charge-grid offset
%     Pc_base : 3 x 6   ideal charge grid
%     L       : Np x 6 x 6 x 3
%
%   With  pbar = P/l_hat  and  rn = pbar - Pc(:,n),  rm = pbar - Pc(:,m):
%
%     L(:,n,m,:) = [ (1 - 3(rn.rm)/|rn|^2) rn + (1 - 3(rn.rm)/|rm|^2) rm ]
%                  / ( |rn|^3 |rm|^3 )
%
%   This is the gradient, with respect to the dimensionless position, of the
%   pairwise charge-interaction kernel (rn.rm)/(|rn|^3 |rm|^3) that appears when
%   B.B is expanded over the six point charges. It is symmetric in (n,m), so
%   only the 21 upper-triangular pairs are evaluated and mirrored.
%
%   The force follows as  f_k = s * g' * L_k * g  with g the charge vector.
%   See BASE_MODEL.
%
%   See also CURRENT_BASE_MODEL, LOAD_FLUX_CALIB.

    validateattributes(P, {'numeric'}, {'2d','ncols',3,'real','finite'}, mfilename, 'P');
    validateattributes(l_hat, {'numeric'}, {'scalar','real','positive','finite'}, mfilename, 'l_hat');
    validateattributes(Pc_base, {'numeric'}, {'size',[3 6],'real','finite'}, mfilename, 'Pc_base');
    e = e(:);
    assert(numel(e) == 17, 'build_L:eSize', 'e must be 17x1 (got %d)', numel(e));

    Pc   = make_Pc(e, Pc_base);          % 3x6 dimensionless charge grid
    pbar = P / l_hat;                    % Np x 3, dimensionless
    Np   = size(pbar, 1);

    % offset vector and norm for every charge
    R  = zeros(Np, 3, 6);
    Rn = zeros(Np, 6);
    for n = 1:6
        d        = pbar - Pc(:,n).';
        R(:,:,n) = d;
        Rn(:,n)  = sqrt(sum(d.^2, 2));
    end
    assert(all(Rn(:) > 0), 'build_L:onCharge', ...
           'a sample point coincides with a charge position (kernel is singular there)');

    L = zeros(Np, 6, 6, 3);
    for n = 1:6
        rn = R(:,:,n);   rn2 = Rn(:,n).^2;   rn3 = Rn(:,n).^3;
        for m = n:6
            rm = R(:,:,m);   rm2 = Rn(:,m).^2;   rm3 = Rn(:,m).^3;
            dt = sum(rn .* rm, 2);                                  % Np x 1
            v  = ((1 - 3*dt./rn2) .* rn + (1 - 3*dt./rm2) .* rm) ./ (rn3 .* rm3);
            L(:,n,m,:) = reshape(v, Np, 1, 1, 3);
            if m ~= n
                L(:,m,n,:) = L(:,n,m,:);                            % symmetric in (n,m)
            end
        end
    end
end

% ---- charge grid Pc_bar = Pc_base + E(e), including the e6z constraint -------
%  Verbatim copy of make_Pc in Flux/Maxwell/function/fitting.m. It MUST stay
%  identical: if the e6z constraint differs, e_hat can never reproduce e_m and
%  the self-consistency check would fail for the wrong reason.
function Pc = make_Pc(e17, Pc_base)
    E = zeros(3, 6);
    E(:,1) = e17(1:3);     E(:,2) = e17(4:6);
    E(:,3) = e17(7:9);     E(:,4) = e17(10:12);
    E(:,5) = e17(13:15);
    E(1,6) = e17(16);      E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);   % e6z constraint (paper step 2c)
    Pc = Pc_base + E;
end
