% i_cons.m -- Ihat_cons: the two-step (constrained) solution of Long's model at every fitting point and direction.
% Long's model (R_norm = 550 um, theoretical K_I, e = 0), unit force r_hat, min ||Ihat||^2 subject to fhat = r_hat and
% Long's region constraint c_x = I_1 + I_2 >= CMIN, c_y = I_3 + I_4 >= CMIN, c_z = I_5 + I_6 <= 0 (CMIN = sqrt(10)/24).
% Positions: the 64 points nn_grid_oct_h15_a45 (0:15:45 um), directions every 6 deg.
% The constrained problem is non-convex; four existing nn_label.m runs (same model and constraint, different starts)
% are combined, per position and direction the feasible solution of least norm:
%   cmin  40 random starts             sc  started from the origin constrained labels
%   rc    position relay from the origin   lc  started from Long's LUT
% Feasible = constraints hold to 1e-6 and Ihat ~= 0 (the run stores 0 where it found nothing).
% Output utils/data/long2016_hexapole_halfcut/i_cons.mat:
%   Ihat_cons 6 x N_phi x N_theta x N_p [A], P_act N_p x 3 [um], ph, th [deg], SRC (run index per entry, 1..4 in RUN order),
%   NOK (number of runs that agree with the least norm to 1e-6), RUN, CMIN
% PHI = true: the same merge for Long's region on the pair sums of Phi_hat, our single-parameter model (VAR flux, 6 dim),
% the 300 blue points nn_grid_sph_r45_mid, three runs (cmin, sc, rc; no LUT) -> phi_cons.mat (Phihat_cons, ...), no eta
% PCO = true: report/406/11.24/11.24.pptx slides 2-3 at the origin (pbar = 0), every 6 deg direction, closed form:
%   r_hat = 6 C dPhi_hat, C = diag(c_x, c_y, c_z), c = pair sums, c_x + c_y + c_z = 0 (exact at pbar = 0 for sum(Phi_hat) = 0);
%   {c_x, c_y, c_z}(phi, theta) = argmin J = 1/2 (c_x^2 + c_y^2 + c_z^2) + 1/72 ||C^-1 r_hat||^2  s.t.  c_x + c_y + c_z = 0;
%   [Phi_1; Phi_3; Phi_5], [Phi_2; Phi_4; Phi_6] = 1/2 c_0 +- 1/12 C^-1 r_hat.
%   J is convex in each of the six sign sectors of (c_x, c_y, c_z) -> one fminsearch per sector, least J kept.
%   Sign: (c, Phi) and (-c, -Phi) give the same force and J; per direction the sign of the nn_sign.m-aligned global
%   optimum nn_label_o_flux_single_sgn.mat is taken (positive inner product), so the surfaces are not cut by sign flips.
%   Check: f_k = Phi_hat' L_k(0) Phi_hat of the single calibration (point charges, e = 0) against r_hat.
%   Output utils/data/long2016_hexapole_halfcut/cons_o.mat: IOPT 6 x N_phi x N_theta (Phi_hat, P1..P6), c0 3 x N_phi x N_theta,
%   J, ph (-90:6:90 elevation), th (0:6:360 azimuth), P_act = [0 0 0], S25.z (6 x 1742 distinct directions), nP (1742 x 1)
if ~exist('PCO','var'), PCO = false; end
% CPOS = '<grid file tag>' (e.g. '_axs5_r150' -> nn_grid_axs5_r150.mat): position-dependent constant constraints,
%   C:\Users\Kuo\Downloads\c_position (1).pdf Sec. 3, 5, 6 (symbols as there).  At every position pbar_i = p_i / l_hat
%   (single calibration, e = 0) ONE c(pbar_i) = [c_x c_y c_z] for all 6 deg directions (phi_j azimuth 0:6:354,
%   theta_m from +z 0:6:180, the 1742 distinct directions, poles once):
%     fhat_k = 1/4 C' L_k C + 1/2 C' L_k dPhihat + 1/4 dPhihat' L_k dPhihat = rhat_k(phi_j, theta_m),  k = x, y, z
%     (C = [c_x c_x c_y c_y c_z c_z]', dPhihat = [dPhihat_x -dPhihat_x dPhihat_y -dPhihat_y dPhihat_z -dPhihat_z]')
%     J_i = 1/2 sum_j sum_m { ||dPhihat(phi_j, theta_m; pbar_i)||^2 + c_x^2 + c_y^2 + c_z^2 },  c_x + c_y + c_z = 0
%   Given c, dPhihat of each direction solves the three force equations (Newton from the linear-term solution,
%   lsqnonlin for the directions Newton leaves); c*(pbar_i) = argmin J_i by fmincon, algorithm 'sqp', linear equality
%   c_x + c_y + c_z = 0, analytic gradient; one start per sign sector (two positive pair sums; -c, -dPhihat give the
%   same forces), least J_i kept.  Points independent (c not a position function here).
%   Origin (pbar = 0): closed form.  dPhihat = C^-1 rhat / 6 gives J = 1/2 [ 1/36 sum_k S_k / c_k^2 + N_dir ||c||^2 ],
%   S_k = sum_jm rhat_k^2; stationarity N_dir c_k - S_k / (36 c_k^3) = lambda (monotone in c_k on each sign), so the two
%   same-sign sums of a sector with equal S are equal: c = a [1 1 -2] (negative axis q, S_p of the other two equal),
%   a^4 = (8 S_p + S_q) / (864 N_dir).  On this grid S_x = S_y (phi uniform), so the sector q = z is closed form; the
%   sectors q = x, y (S unequal) are solved by sqp and compared, least J kept.
%   CPI: indices of the points to solve ([] = all); an existing c_pos<CPOS>.mat is updated (unsolved points NaN).
%   CREL = true: relay outward from the origin along each half-axis: the parent of a point is the point one shell
%   closer on the same half-axis (the origin for the first shell, PAR = index); points solved in order of |p|; the
%   only sqp start is the parent's c*, and the inner Newton starts from the parent's dPhihat (the linear-term start and
%   lsqnonlin for the directions that do not converge from it).  The origin itself stays the closed form.
%   Output utils/data/long2016_hexapole_halfcut/c_pos<CPOS>.mat: c_star 3 x N_p, J N_p x 1, dPhihat 3 x N_dir x N_p,
%   pbar N_p x 3, P_act N_p x 3 [um], phi, theta N_dir x 1 [deg], l_hat [um], NF N_p x 1 (directions without a solution)
if ~exist('CPOS','var'), CPOS = ''; end
if ~exist('CPI','var'), CPI = []; end
if ~exist('CREL','var'), CREL = false; end
clearvars -except PHI PCO CPOS CPI CREL;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
DAT  = fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut');
if ~isempty(CPOS)
    addpath(fullfile(CALROOT, 'function'));
    cal = load_flux_calib('long2016_hexapole_halfcut', 'tip40um', 'calib_current_maxwell_axshN43_R150_single', 'current');
    G   = load(fullfile(DAT, ['nn_grid' CPOS '.mat']), 'P_act');  P_act = G.P_act;  N_p = size(P_act, 1);
    l_hat = cal.l_m;  pbar = P_act / l_hat;
    [TH, PH] = ndgrid(0:6:180, 0:6:354);  U = (TH > 0 & TH < 180) | PH == 0;     % poles once
    theta = TH(U);  phi = PH(U);  N_dir = numel(phi);
    rhat = [sind(theta).*cosd(phi), sind(theta).*sind(phi), cosd(theta)].';  rhat(abs(rhat) < 1e-12) = 0;
    T  = kron(eye(3), [1; -1]);  Uc = kron(eye(3), [1; 1]);        % dPhihat = T dPhihat_xyz, C = Uc c
    Lb = build_L(P_act, l_hat, cal.e_m, cal.Pc_base);                % N_p x 6 x 6 x 3
    % fmincon (sqp) on c in R^3 with Aeq = [1 1 1], beq = 0; J(c) and its gradient (implicit differentiation of the
    % force equations) from cpos_J_; one start per sign sector (two positive sums; -c is the same solution), least J kept
    C0 = 0.3 * [1 1 -2; 1 -2 1; -2 1 1].' / 2;                     % starts, columns (|c| as at the origin)
    opt = optimoptions('fmincon', 'Algorithm', 'sqp', 'SpecifyObjectiveGradient', true, 'Display', 'off', ...
                       'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-12, 'ConstraintTolerance', 1e-12, ...
                       'MaxIterations', 500, 'MaxFunctionEvaluations', 3000);
    OUTF = fullfile(DAT, ['c_pos' CPOS '.mat']);
    if exist(OUTF, 'file')                                          % update an existing file point by point
        load(OUTF, 'c_star', 'J', 'dPhihat', 'NF', 'SEC', 'EXITF', 'JS', 'CF');
    else
        c_star = nan(3, N_p);  J = nan(N_p, 1);  dPhihat = nan(3, N_dir, N_p);  NF = nan(N_p, 1);
        SEC = nan(N_p, 1);  EXITF = nan(N_p, 1);  JS = nan(N_p, size(C0, 2));  CF = false(N_p, 1);
    end
    if isempty(CPI), CPI = 1:N_p; end
    S_k = sum(rhat.^2, 2);                                          % S_x, S_y, S_z
    rp = vecnorm(P_act, 2, 2);  PAR = zeros(N_p, 1);                % relay parent: one shell closer, same half-axis
    for i = 1:N_p
        if rp(i) == 0, continue, end
        u_ = P_act(i,:) / rp(i);  k_ = find(abs(vecnorm(P_act - u_ .* rp, 2, 2)) < 1e-9 & rp < rp(i) - 1e-9);
        [~, m_] = max(rp(k_));  if ~isempty(k_), PAR(i) = k_(m_); else, PAR(i) = find(rp == 0, 1); end
    end
    [~, o_] = sort(rp(CPI));  CPI = CPI(o_);                        % inner points first
    for i = CPI(:).'
        L = cell(1, 3);
        for k = 1:3, L_k = squeeze(Lb(i,:,:,k));  L{k} = (L_k + L_k.') / 2; end
        Jb = Inf;  CF(i) = false;  Dp = [];
        if CREL && rp(i) > 0                                        % relay: the parent's c* and dPhihat
            assert(~isnan(J(PAR(i))), 'i_cons:relay', 'parent %d of point %d not solved', PAR(i), i);
            Dp = dPhihat(:, :, PAR(i));
            [c_, Jv, ef] = fmincon(@(c) cpos_J_(c, L, T, Uc, rhat, Dp), c_star(:, PAR(i)), [], [], [1 1 1], 0, [], [], [], opt);
            JS(i, :) = NaN;  SEC(i) = NaN;  EXITF(i) = ef;  Jb = Jv;
        end
        for s = 1:size(C0, 2) * ~(CREL && rp(i) > 0)
            q = find(C0(:, s) < 0);  pp = setdiff(1:3, q);         % negative axis q, the two positive ones pp
            if all(pbar(i, :) == 0) && abs(S_k(pp(1)) - S_k(pp(2))) < 1e-9 * N_dir
                a = ((8*S_k(pp(1)) + S_k(q)) / (864*N_dir))^(1/4);  % closed form at the origin
                c_s = zeros(3, 1);  c_s(pp) = a;  c_s(q) = -2*a;  Jv = cpos_J_(c_s, L, T, Uc, rhat);  ef = NaN;  cf = true;
            else
                [c_s, Jv, ef] = fmincon(@(c) cpos_J_(c, L, T, Uc, rhat), C0(:, s), [], [], [1 1 1], 0, [], [], [], opt);
                cf = false;
            end
            JS(i, s) = Jv;
            if Jv < Jb, Jb = Jv;  c_ = c_s;  SEC(i) = s;  EXITF(i) = ef;  CF(i) = cf; end
        end
        [Jv, ~, D, nf] = cpos_J_(c_, L, T, Uc, rhat, Dp);
        c_star(:, i) = c_;  J(i) = Jv;  dPhihat(:, :, i) = D;  NF(i) = nf;
        fprintf('p_%2d = [%6.1f %6.1f %6.1f] um: c* = [%8.5f %8.5f %8.5f], J = %.4f (sectors %s), closed form %d, unsolved %d\n', ...
                i, P_act(i,:), c_, Jv, mat2str(JS(i,:), 6), CF(i), nf);
    end
    save(OUTF, 'c_star', 'J', 'dPhihat', 'pbar', 'P_act', 'phi', 'theta', 'l_hat', 'NF', 'SEC', 'EXITF', 'JS', 'C0', 'CF', 'S_k', 'PAR');
    fprintf('-> %s (solved points: %d of %d)\n', OUTF, nnz(~isnan(J)), N_p);
    return
end
if PCO
    ph = -90:6:90;  th = 0:6:360;  P_act = [0 0 0];
    IOPT = zeros(6, numel(ph), numel(th));  c0 = zeros(3, numel(ph), numel(th));  J = zeros(numel(ph), numel(th));
    SG = [1 1 -1; 1 -1 1; -1 1 1; -1 -1 1; -1 1 -1; 1 -1 -1];      % six sign sectors of (c_x, c_y, c_z), sum = 0
    opt = optimset('TolX', 1e-12, 'TolFun', 1e-14, 'MaxIter', 4000, 'MaxFunEvals', 8000, 'Display', 'off');
    for a = 1:numel(ph)
        for b = 1:numel(th)
            r_hat = [cosd(ph(a))*cosd(th(b)); cosd(ph(a))*sind(th(b)); sind(ph(a))];
            r_hat(abs(r_hat) < 1e-12) = 0;
            Jbest = Inf;
            for s = 1:size(SG, 1)
                sg = SG(s, :).';
                Jf = @(u) jfun_(u, sg, r_hat);
                [u, Jv] = fminsearch(Jf, [0 0], opt);
                if Jv < Jbest - 1e-12, Jbest = Jv;  cb = cvec_(u, sg); end
            end
            if nnz(cb > 0) < 2, cb = -cb; end                     % sign convention: two positive pair sums
            C = diag(cb);  q = C \ r_hat;
            Po = cb/2 + q/12;  Pe = cb/2 - q/12;                  % [Phi_1; Phi_3; Phi_5], [Phi_2; Phi_4; Phi_6]
            IOPT(:, a, b) = [Po(1); Pe(1); Po(2); Pe(2); Po(3); Pe(3)];
            c0(:, a, b) = cb;  J(a, b) = Jbest;
        end
    end
    % sign per direction matched to the nn_sign.m-aligned global optimum (same force for Phi and -Phi)
    Gs = load(fullfile(DAT, 'nn_label_o_flux_single_sgn.mat'), 'IOPT', 'ph', 'th');
    assert(isequal(Gs.ph, ph) && isequal(Gs.th, th), 'i_cons:pco', 'aligned labels on another grid');
    flp = sum(IOPT .* Gs.IOPT(:, :, :, 1), 1) < 0;
    IOPT = IOPT .* (1 - 2*flp);  c0 = c0 .* (1 - 2*flp);
    fprintf('cons_o: %d of %d grid entries flipped to the aligned optimum\n', nnz(flp), numel(flp));
    % force check with the single calibration at the origin (point charges, e = 0)
    addpath(fullfile(CALROOT, 'function'));
    cal = load_flux_calib('long2016_hexapole_halfcut', 'tip40um', 'calib_current_maxwell_axshN43_R150_single', 'current');
    Lb  = build_L(P_act, cal.l_m, cal.e_m, cal.Pc_base);
    [PH, TH] = ndgrid(ph, th);  U = TH < 360 & (abs(PH) < 90 | TH == 0);
    z  = reshape(IOPT, 6, []);  z = z(:, U(:));
    R  = [cosd(PH(U)).*cosd(TH(U)), cosd(PH(U)).*sind(TH(U)), sind(PH(U))].';
    f  = zeros(3, size(z, 2));
    for k = 1:3, L_k = squeeze(Lb(1,:,:,k));  f(k,:) = sum(z .* (((L_k + L_k.')/2) * z), 1); end
    e  = 100 * vecnorm(f - R) ./ vecnorm(R);
    nP = vecnorm(z).';
    S25 = struct('z', z, 'phi', TH(U), 'theta', 90 - PH(U), 'e', e(:));      % fourier2d_new angles: phi azimuth, theta from +z
    save(fullfile(DAT, 'cons_o.mat'), 'IOPT', 'c0', 'J', 'ph', 'th', 'P_act', 'S25', 'nP');
    fprintf('cons_o: %d directions | ||Phi_hat|| mean %.4f, std %.4f, min %.4f, median %.4f, max %.4f | sum(Phi) max %.1e | force error max %.2e %%\n', ...
            numel(nP), mean(nP), std(nP), min(nP), median(nP), max(nP), max(abs(sum(z, 1))), max(e));
    return
end
RUN  = {'cmin', 'rc', 'lc', 'sc'};  FMT = 'nn_label_oct_h15_a45_%s_lg_flux5_single.mat';
if PHI, RUN = {'cmin', 'rc', 'sc'};  FMT = 'nn_label_sph_r45_mid_%s_flux_single.mat'; end
for i = 1:numel(RUN)
    S = load(fullfile(DAT, sprintf(FMT, RUN{i})), 'IOPT', 'P_act', 'ph', 'th', 'CMIN', 'LONG', 'CREG');
    assert(S.CREG && (S.LONG || PHI), 'i_cons:run', '%s is not a constrained run', RUN{i});
    if i == 1
        P_act = S.P_act;  ph = S.ph;  th = S.th;  CMIN = S.CMIN;  sz = size(S.IOPT);
        X = zeros(6, prod(sz(2:end)), numel(RUN));  nrm = inf(numel(RUN), prod(sz(2:end)));
    end
    assert(isequal(S.P_act, P_act) && isequal(S.ph, ph) && isequal(S.th, th) && S.CMIN == CMIN, 'i_cons:grid', '%s: grid differs', RUN{i});
    c  = reshape(S.IOPT, 6, []);  n = vecnorm(c);
    ok = c(1,:) + c(2,:) >= CMIN - 1e-6 & c(3,:) + c(4,:) >= CMIN - 1e-6 & c(5,:) + c(6,:) <= 1e-6 & n > 0;
    X(:,:,i) = c;  nrm(i, ok) = n(ok);
end
[nmin, SRC] = min(nrm, [], 1);
assert(all(isfinite(nmin)), 'i_cons:none', '%d entries without a feasible solution', nnz(~isfinite(nmin)));
NOK = sum(abs(nrm - nmin) <= 1e-6 * nmin, 1);
Ihat_cons = zeros(6, numel(nmin));
for i = 1:numel(RUN), Ihat_cons(:, SRC == i) = X(:, SRC == i, i); end
Ihat_cons = reshape(Ihat_cons, sz);  SRC = reshape(SRC, sz(2:end));  NOK = reshape(NOK, sz(2:end));
if PHI
    Phihat_cons = Ihat_cons;  [PH, TH] = ndgrid(ph, th);  U = TH < 360 & (abs(PH) < 90 | TH == 0);
    nP = reshape(squeeze(vecnorm(Phihat_cons, 2, 1)), [], size(Phihat_cons, 4));  nP = nP(U(:), :);
    save(fullfile(DAT, 'phi_cons.mat'), 'Phihat_cons', 'P_act', 'ph', 'th', 'SRC', 'NOK', 'RUN', 'CMIN');
    fprintf('Phihat_cons: %d positions x %d directions | from %s: %s | only one run at the least norm: %d | ||Phi|| mean %.4f, median %.4f, max %.4f -> phi_cons.mat\n', ...
            size(nP, 2), size(nP, 1), strjoin(RUN, '/'), mat2str(histcounts(SRC(repmat(U, 1, 1, size(nP, 2))), 0.5:numel(RUN)+0.5)), ...
            nnz(NOK(repmat(U, 1, 1, size(nP, 2))) == 1), mean(nP(:)), median(nP(:)), max(nP(:)));
    return
end
save(fullfile(DAT, 'i_cons.mat'), 'Ihat_cons', 'P_act', 'ph', 'th', 'SRC', 'NOK', 'RUN', 'CMIN');
% eta = ||Ihat_cons|| / ||Ihat_opt|| (>= 1), Ihat_opt = the proved global minimum without constraint (same model, positions,
% directions; nn_label_oct_h15_a45_g_lg_flux5_single.mat), on the distinct directions (theta = 360 dropped, poles once)
G = load(fullfile(DAT, 'nn_label_oct_h15_a45_g_lg_flux5_single.mat'), 'IOPT', 'P_act', 'ph', 'th');
assert(isequal(G.P_act, P_act) && isequal(G.ph, ph) && isequal(G.th, th), 'i_cons:glob', 'global labels on another grid');
[PH, TH] = ndgrid(ph, th);  U = TH < 360 & (abs(PH) < 90 | TH == 0);
eta = squeeze(vecnorm(Ihat_cons, 2, 1) ./ vecnorm(G.IOPT, 2, 1));       % N_phi x N_theta x N_p
eta = reshape(eta, [], size(eta, 3));  eta = eta(U(:), :);              % N_j x N_p
save(fullfile(DAT, 'i_cons.mat'), 'eta', '-append');
% eta_best = ||Ihat_best|| / ||Ihat_opt||, Ihat_best = the physics-loss fit (2D Fourier, M = N = 4, fit_lg_f2d_I.mat) evaluated
% by ls_fit.m TEST at the 64 points every 3 deg (test_lg_f2d_I_all.mat, S25.z), taken on the same 1742 directions
T  = load(fullfile(DAT, 'test_lg_f2d_I_all.mat'), 'S25');  T = T.S25;
assert(isequal(T.P, P_act), 'i_cons:best', 'test positions differ');
[~, k_] = ismember([PH(U) TH(U)], [T.phi T.theta], 'rows');  assert(all(k_ > 0), 'i_cons:dirs', 'grid direction missing in the test file');
nO = reshape(squeeze(vecnorm(G.IOPT, 2, 1)), [], numel(P_act)/3);  nO = nO(U(:), :);   % N_j x N_p
% Ihat_best is NOT rescaled to 1 pN (eta < 1 where its force is slightly short or off r_hat); ||f|| is kept for reference,
% f_k = Ihat' K_I' L_k K_I Ihat (Long's model: l = 550 um, e = 0, K_I = I_6 - 11'/6)
addpath(fullfile(CALROOT, 'function'));
cal = load_flux_calib('long2016_hexapole_halfcut', 'tip40um', 'calib_current_maxwell_axshN43_R150_single', 'current');
cal.l_m = 550;  cal.e_m(:) = 0;
U5 = null(ones(1, 6));  K_I = U5 * U5.';  Lb = build_L(P_act, cal.l_m, cal.e_m, cal.Pc_base);   % N_p x 6 x 6 x 3
Ib = T.z(:, k_, :);  fn = zeros(size(Ib, 2), size(Ib, 3));               % ||f||, N_j x N_p
for p = 1:size(Ib, 3)
    f_ = zeros(size(Ib, 2), 3);
    for c = 1:3, L_k = squeeze(Lb(p,:,:,c));  Q_k = K_I.' * ((L_k + L_k.')/2) * K_I;  f_(:,c) = sum(Ib(:,:,p) .* (Q_k * Ib(:,:,p)), 1).'; end
    fn(:, p) = vecnorm(f_, 2, 2);
end
eta_best = squeeze(vecnorm(Ib, 2, 1)) ./ nO;                             % N_j x N_p, not scaled (decided 2026-11-10)
save(fullfile(DAT, 'i_cons.mat'), 'eta_best', 'fn', '-append');
fprintf('eta_best: %d values | mean %.4f, std %.4f, min %.4f, median %.4f, max %.4f | eta = 1 (1e-6): %.2f %%, < 1: %.2f %%\n', ...
        numel(eta_best), mean(eta_best(:)), std(eta_best(:)), min(eta_best(:)), median(eta_best(:)), max(eta_best(:)), ...
        100*mean(abs(eta_best(:) - 1) < 1e-6), 100*mean(eta_best(:) < 1 - 1e-6));
fprintf('eta: %d values | mean %.4f, std %.4f, min %.4f, median %.4f, max %.4f | eta < 1 - 1e-6: %d, eta = 1 (1e-6): %.2f %%\n', ...
        numel(eta), mean(eta(:)), std(eta(:)), min(eta(:)), median(eta(:)), max(eta(:)), nnz(eta < 1 - 1e-6), 100*mean(abs(eta(:) - 1) < 1e-6));
fprintf('Ihat_cons: %d positions x %d x %d directions | only one run reached the least norm: %d entries | ||Ihat_cons|| mean %.4f, max %.4f -> i_cons.mat\n', ...
        sz(4), sz(2), sz(3), nnz(NOK == 1), mean(nmin), max(nmin));

function c = cvec_(u, sg)
% pair sums inside the sign sector sg: the two equal-sign entries exp(u), the odd one closes c_x + c_y + c_z = 0
c = zeros(3, 1);  o = find(sg ~= mode(sg));  m = find(sg == mode(sg));
c(m) = sg(m) .* exp(u(:));  c(o) = -sum(c(m));
end

function J = jfun_(u, sg, r_hat)
% J = 1/2 (c_x^2 + c_y^2 + c_z^2) + 1/72 ||C^-1 r_hat||^2
c = cvec_(u, sg);
J = 0.5 * sum(c.^2) + sum((r_hat ./ c).^2) / 72;
end

function [J, gJ, D, nf] = cpos_J_(c, L, T, Uc, rhat, Dp)
% Dp (optional, 3 x N_dir): Newton start per direction (relay: the parent point's dPhihat); [] = the linear-term start
% c_position (1).pdf Sec. 6: J = 1/2 sum over directions {||dPhihat||^2 + ||c||^2}; dPhihat (3 x N_dir) from
% fhat_k = 1/4 C' L_k C + 1/2 C' L_k T d + 1/4 d' T' L_k T d = rhat_k: Newton from the linear-term solution, then
% lsqnonlin (from the same start) for the directions Newton leaves; a direction still without a solution -> J = Inf
% (nf = their number).  gJ = dJ/dc by implicit differentiation: dd/dc = -(dF/dd)^-1 dF/dc,
% dF_k/dd = 1/2 C' L_k T + 1/2 d' T' L_k T,  dF_k/dc = 1/2 Uc' L_k (C + T d).
c = c(:);  C = Uc * c;  N_dir = size(rhat, 2);  g = zeros(3, 1);  B = zeros(3);  A = cell(1, 3);
for k = 1:3
    g(k) = C.' * L{k} * C / 4;  B(k, :) = (C.' * L{k} * T) / 2;  A{k} = T.' * L{k} * T;
end
Fun = @(D) [g(1) + B(1,:)*D + sum(D.*(A{1}*D), 1)/4 - rhat(1,:); ...
            g(2) + B(2,:)*D + sum(D.*(A{2}*D), 1)/4 - rhat(2,:); ...
            g(3) + B(3,:)*D + sum(D.*(A{3}*D), 1)/4 - rhat(3,:)];
D0 = B \ (rhat - g);  D = D0;                               % linear-term start (exact at pbar = 0)
if nargin > 5 && ~isempty(Dp), D = Dp; end
for it = 1:30
    Jac = jac_(D, B, A);
    dD = reshape(pagemldivide(Jac, reshape(Fun(D), 3, 1, N_dir)), 3, N_dir);
    dD(~isfinite(dD)) = 0;  D = D - dD;
    if max(abs(dD(:))) < 1e-13, break, end
end
bad = find(~all(isfinite(D), 1) | vecnorm(Fun(D)) > 1e-10);
if ~isempty(bad)
    o = optimoptions('lsqnonlin', 'Display', 'off', 'FunctionTolerance', 1e-20, 'StepTolerance', 1e-15, ...
                     'OptimalityTolerance', 1e-16, 'MaxIterations', 400);
    for n = bad
        fn = @(d) [g(1) + B(1,:)*d + d.'*A{1}*d/4; g(2) + B(2,:)*d + d.'*A{2}*d/4; g(3) + B(3,:)*d + d.'*A{3}*d/4] - rhat(:, n);
        d0 = D(:, n);  d0(~isfinite(d0)) = 0;
        d_ = lsqnonlin(fn, d0, [], [], o);                                % from where Newton stopped (parent / linear start)
        if norm(fn(d_)) > 1e-10, d_ = lsqnonlin(fn, D0(:, n), [], [], o); end   % then from the linear-term start
        D(:, n) = d_;
    end
end
res = vecnorm(Fun(D));  nf = nnz(~all(isfinite(D), 1) | res > 1e-9);
J  = 0.5 * sum(sum(D.^2, 1) + sum(c.^2));
gJ = N_dir * c;
if nf > 0, J = Inf;  gJ = zeros(3, 1);  return, end
Jac = jac_(D, B, A);  Fc = zeros(3, 3, N_dir);              % dF/dc per direction
for k = 1:3, Fc(k, :, :) = reshape(0.5 * (Uc.' * L{k} * (C + T * D)), 1, 3, N_dir); end
dDdc = -pagemldivide(Jac, Fc);                               % 3 x 3 x N_dir
dDdc(:, :, ~all(isfinite(reshape(dDdc, 9, N_dir)), 1)) = 0;   % singular dF/dd (fold): that direction left out of gJ
gJ = gJ + sum(squeeze(sum(reshape(D, 3, 1, N_dir) .* dDdc, 1)), 2);
end

function Jac = jac_(D, B, A)
% dF/dd per direction (3 x 3 x N_dir): row k = B(k,:) + 1/2 (A_k d)'
N_dir = size(D, 2);  Jac = zeros(3, 3, N_dir);
for k = 1:3, Jac(k, :, :) = reshape(B(k, :).' + A{k} * D / 2, 1, 3, N_dir); end
end
