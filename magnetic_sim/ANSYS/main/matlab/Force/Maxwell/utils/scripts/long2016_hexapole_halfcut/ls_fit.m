%% ls_fit.m -- least-squares inverse model (Long et al., 2021), symbols as in C:\Users\Kuo\Downloads\LS_inverse_model.pdf.
%  For every direction (phi_j, theta_j) and coil i the unit-force optimal current is fitted by a quadratic form
%      Ihat_opt(phi_j, theta_j, pbar)_i ~ P' D_LS_i(phi_j, theta_j) P,   P = [xbar ybar zbar 1]',  pbar = p / l_hat   (Sec. 3)
%  linear in the 10 coefficients d_i = [d44 d14 d24 d34 d12 d13 d23 d11 d22 d33]':
%      tau(pbar) = [1, 2xbar, 2ybar, 2zbar, 2xbar*ybar, 2xbar*zbar, 2ybar*zbar, xbar^2, ybar^2, zbar^2]'
%  Stacking the N_p positions (Sec. 4):  T = [tau(pbar_1)'; ...; tau(pbar_N_p)']  (N_p x 10),  y_i = Ihat_opt at them (N_p x 1)
%  Normal equations (Sec. 5):  d_i = (T'T)^-1 T' y_i;  (T'T)^-1 T' is shared by every direction and coil (Sec. 6).
%  Labels y_i: nn_label_sgn.mat (Gurobi-certified I_opt_hat), re-aligned in sign within each direction to the label
%  at pbar = 0 (the whole 6-vector flips; the force is unchanged).
%  LAB = 'chain' (REG = 'cube' only): labels from chain.mat instead, i.e. the position continuation from the centre's
%  global minimum (one branch per direction); output ls_fit_chain.mat.
%  Actuator frame throughout; coil i = pole P1..P6.
%  Directions j = 1..N_j: the distinct directions of the label grid (theta = 360 dropped, phi = +-90 kept once at theta = 0),
%  ordered phi fastest, then theta.
%  REG selects the fitting positions:  'cube' = |x|, |y|, |z| <= A (default A = 45 um: the 90 um cube of Long's
%                                               dissertation p. 49, all octants)
%                                      'all'  = every grid node, |p| <= 150 um (N_p = 515)
%                                      'oct'  = first octant x, y, z >= 0
%  Output: utils/data/long2016_hexapole_halfcut/ls_fit[_all|_oct|_chain].mat
%    D_LS            4 x 4 x N_j x N_I   [A]    D_LS_i(phi_j, theta_j), entry (:,:,j,i)
%    d               10 x N_j x N_I      [A]    the same coefficients as d_i
%    phi_j, theta_j  N_j x 1             [deg]  directions
%    pbar            N_p x 3                    fitting positions (normalised)
%    l_hat           [um]
%    e_ij            N_j x N_p           [%]    relative force error at position pbar_i, direction j:
%                    e_ij = || fhat_LS(pbar_i, theta_j, phi_j) - fhat_d(pbar_i, theta_j, phi_j) || / || fhat_d || * 100 %
%                    fhat_d = r_hat(phi_j, theta_j) (unit force, 1 pN); fhat_LS_k = Ihat_LS' Q_k Ihat_LS with
%                    Ihat_LS = P' D_LS P (all six coils) and Q_k = F H_I' L_k(pbar) F H_I as in nn_label.m
%  LOSS = 'phys' (C:\Users\Kuo\Downloads\LS_inverse_model (1).pdf, Sec. 5-6, 9): for every direction j,
%      d(phi_j, theta_j) = argmin_d sum_n [ || fhat(pbar_n, Ihat_LS) - r_hat ||^2 + LAM || Ihat_LS ||^2 ]
%    Levenberg-Marquardt on the residual e_j(d) (9 N_p x 1) with the Jacobian of Sec. 6, started from the label
%    fit d^(0) = (T'T)^-1 T' y_i (labels LAB, sign-aligned).  The polynomial variable of tau is u = p / A (half side,
%    A = 45 um) instead of pbar (same function space; T'T well conditioned); Q_k stays at pbar = p / l_hat.
%    D_LS, d then refer to u.  Output ls_fit_phys.mat, adds eta (N_j x N_p) = ||Ihat_opt||^2 / ||Ihat_LS||^2
%    (Ihat_opt: Gurobi global minimum of nn_label_sgn.mat), LAM, e_ij0 (force error of the label fit d^(0)).
if ~exist('REG','var'), REG = 'cube'; end
if ~exist('A','var'), A = 45; end
if ~exist('LAB','var'), LAB = 'grb'; end
if ~exist('LOSS','var'), LOSS = 'label'; end
if ~exist('LAM','var'), LAM = 0.01; end
clearvars -except REG A LAB LOSS LAM;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
addpath(fullfile(CALROOT, 'function'));
MODEL = 'long2016_hexapole_halfcut';  GEOM = 'tip40um';
DAT  = fullfile(CALROOT, 'utils', 'data', MODEL);
S    = load(fullfile(DAT, 'nn_label_sgn.mat'), 'P_act', 'ph', 'th', 'IOPT', 'l_hat', 'MATNAME');
l_hat = S.l_hat;
[N_I, N_phi, N_theta, ~] = size(S.IOPT);

switch REG
    case 'cube', n = all(abs(S.P_act) <= A + 1e-9, 2);     OUT = 'ls_fit.mat';
    case 'all',  n = true(size(S.P_act, 1), 1);           OUT = 'ls_fit_all.mat';
    case 'oct',  n = all(S.P_act >= -1e-9, 2);             OUT = 'ls_fit_oct.mat';
    otherwise,   error('ls_fit:REG', 'REG = ''cube'' | ''all'' | ''oct''');
end
pbar = S.P_act(n, :) / l_hat;
N_p  = size(pbar, 1);
pt   = pbar;                                                  % variable of tau
if strcmp(LOSS, 'phys'), pt = pbar * l_hat / A; end            % u = p / A
xbar = pt(:,1);  ybar = pt(:,2);  zbar = pt(:,3);

% ---- distinct directions j ----
[PH, TH] = ndgrid(S.ph, S.th);
U = TH < 360 & (abs(PH) < 90 | TH == 0);
phi_j = PH(U);  theta_j = TH(U);  N_j = numel(phi_j);

% ---- Sec. 3-4: T (N_p x 10) ----
T = [ones(N_p,1), 2*xbar, 2*ybar, 2*zbar, 2*xbar.*ybar, 2*xbar.*zbar, 2*ybar.*zbar, xbar.^2, ybar.^2, zbar.^2];
assert(N_p >= 10 && rank(T) == 10, 'ls_fit:rank', 'T must have full column rank (N_p = %d, rank %d)', N_p, rank(T));

% ---- Sec. 5: d_i = (T'T)^-1 T' y_i for every direction and coil at once ----
I = reshape(S.IOPT(:,:,:,n), N_I, N_phi*N_theta, N_p);
I = I(:, U(:), :);                                           % N_I x N_j x N_p
I_opt = I;                                                   % global minimum (for eta)
if strcmp(LAB, 'chain')
    assert(strcmp(REG, 'cube'), 'ls_fit:LAB', 'LAB = ''chain'' needs REG = ''cube''');
    C = load(fullfile(DAT, 'chain.mat'), 'I_chain', 'pbar', 'phi_j', 'theta_j', 'A');
    [hit, loc] = ismembertol(pbar, C.pbar, 1e-12, 'ByRows', true);
    assert(all(hit) && C.A == A && isequal(C.phi_j, phi_j) && isequal(C.theta_j, theta_j), 'ls_fit:chain', 'chain.mat does not match');
    I = C.I_chain(:, :, loc);
    OUT = 'ls_fit_chain.mat';
end
% sign alignment within each direction: Ihat_opt and -Ihat_opt give the same force, so the label at every position
% is flipped as a whole 6-vector when it points against the label at the centre pbar = 0 of the same direction
p0 = find(all(pbar == 0, 2));
assert(isscalar(p0), 'ls_fit:centre', 'the fitting positions must contain pbar = 0');
I = I .* sign(sum(I .* I(:, :, p0), 1) + (sum(I .* I(:, :, p0), 1) == 0));
y = reshape(permute(I, [3 2 1]), N_p, []);                    % column = y_i(phi_j, theta_j), order (j, i)
d = (T.'*T) \ (T.'*y);                                        % 10 x (N_j N_I)
Ihat_LS = reshape(T*d, N_p, N_j, N_I);                      % P' D_LS P at the fitting positions

% ---- force model Q_k(pbar_n) (as nn_label.m) ----
mgB = 0.0451;  UF = 1e3;                                       % M g_B [A um^2/mT]; unit factor -> pN
cal = load_flux_calib(MODEL, GEOM, S.MATNAME, 'current');
assert(abs(cal.l_m - l_hat) < 1e-9, 'ls_fit:lhat', 'l_hat differs from the label run');
gF  = UF * mgB / (2*cal.l_m) * cal.gB_m^2;                     % F g_I [pN/A^2]
FH  = sqrt(gF) * cal.Mbar_m;                                   % F H_I [sqrt(pN)/A]
L   = build_L(pbar * l_hat, cal.l_m, cal.e_m, cal.Pc_base);    % N_p x 6 x 6 x 3
Q   = zeros(6, 6, 3, N_p);
for p = 1:N_p
    for k = 1:3, Q_k = FH.' * squeeze(L(p,:,:,k)) * FH;  Q(:,:,k,p) = (Q_k + Q_k.')/2; end
end
fhat_d = [cosd(phi_j).*cosd(theta_j), cosd(phi_j).*sind(theta_j), sind(phi_j)];   % N_j x 3

% ---- LOSS = 'phys': Sec. 5-6, one nonlinear least-squares problem per direction ----
e_ij0 = [];
if strcmp(LOSS, 'phys')
    e_ij0 = force_err(Ihat_LS, Q, fhat_d);                   % label fit d^(0)
    d0 = reshape(d, 10, N_j, N_I);
    opt = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt', 'SpecifyObjectiveGradient',true, ...
                       'FunctionTolerance',1e-14, 'StepTolerance',1e-14, 'MaxIterations',2000, 'Display','off');
    t0 = tic;
    for j = 1:N_j
        dj = lsqnonlin(@(x) e_j(x, T, Q, fhat_d(j,:).', LAM), reshape(d0(:,j,:), [], 1), [], [], opt);
        d0(:,j,:) = reshape(dj, 10, 1, N_I);
        if mod(j, 500) == 0, fprintf('  %d / %d directions, %.0f s\n', j, N_j, toc(t0)); end
    end
    d = reshape(d0, 10, []);
    Ihat_LS = reshape(T*d, N_p, N_j, N_I);
    OUT = 'ls_fit_phys.mat';
end

% ---- Sec. 3: fill D_LS_i from d_i = [d44 d14 d24 d34 d12 d13 d23 d11 d22 d33] ----
d44 = d(1,:);  d14 = d(2,:);  d24 = d(3,:);  d34 = d(4,:);  d12 = d(5,:);
d13 = d(6,:);  d23 = d(7,:);  d11 = d(8,:);  d22 = d(9,:);  d33 = d(10,:);
D_LS = reshape([d11; d12; d13; d14;  d12; d22; d23; d24;  d13; d23; d33; d34;  d14; d24; d34; d44], 4, 4, N_j, N_I);
d    = reshape(d, 10, N_j, N_I);

% ---- relative force error e_ij; Sec. 9 eta = ||Ihat_opt||^2 / ||Ihat_LS||^2 ----
e_ij = force_err(Ihat_LS, Q, fhat_d);
eta  = squeeze(vecnorm(I_opt, 2, 1)).^2 ./ squeeze(vecnorm(permute(Ihat_LS, [3 2 1]), 2, 1)).^2;   % N_j x N_p

% self-check: P' D_LS P reproduces the label column of one direction and coil through T d
Pn = [pt, ones(N_p,1)];  jc = 7;  ic = 4;
assert(max(abs(sum((Pn*D_LS(:,:,jc,ic)).*Pn, 2) - T*d(:,jc,ic))) < 1e-12, 'ls_fit:D', 'D_LS fill does not match d');
a = find(S.ph == phi_j(jc));  b = find(S.th == theta_j(jc));
if strcmp(LAB, 'grb') && strcmp(LOSS, 'label')
    assert(isequal(abs(squeeze(S.IOPT(ic, a, b, n))), abs(y(:, (ic-1)*N_j + jc))), 'ls_fit:y', 'label column order mismatch');
end

save(fullfile(DAT, OUT), 'D_LS', 'd', 'phi_j', 'theta_j', 'pbar', 'l_hat', 'e_ij', 'eta', 'e_ij0', 'REG', 'A', 'LAB', 'LOSS', 'LAM', ...
     'N_p', 'N_j', 'N_I');
fprintf('REG %s: N_p = %d, N_j = %d, N_I = %d, cond(T''T) = %.3g -> %s\n', REG, N_p, N_j, N_I, cond(T.'*T), fullfile(DAT, OUT));
fprintf('e_ij [%%]: median %.3g, 95th pct %.3g, max %.3g | eta: median %.4f, min %.4f\n', median(e_ij(:)), prctile(e_ij(:), 95), ...
        max(e_ij(:)), median(eta(:)), min(eta(:)));

% ---------------------------------------------------------------------------------------
function e_ij = force_err(Ihat_LS, Q, fhat_d)
% e_ij = || fhat(pbar_n, Ihat_LS) - r_hat || / || r_hat || * 100 %,  N_j x N_p
[N_p, N_j, ~] = size(Ihat_LS);
e_ij = zeros(N_j, N_p);
for p = 1:N_p
    X = squeeze(Ihat_LS(p, :, :));                             % N_j x N_I
    fhat_LS = zeros(N_j, 3);
    for k = 1:3, fhat_LS(:,k) = sum((X*Q(:,:,k,p)) .* X, 2); end
    e_ij(:, p) = vecnorm(fhat_LS - fhat_d, 2, 2) ./ vecnorm(fhat_d, 2, 2) * 100;
end
end

function [e, Jac] = e_j(dv, T, Q, r, LAM)
% Sec. 6: e_j(d) = [fhat(pbar_n, Ihat_LS) - r_hat (n = 1..N_p); sqrt(LAM) Ihat_LS (n = 1..N_p)],  d = [d_1; ...; d_6]
N_p = size(T, 1);
I_LS = T * reshape(dv, 10, 6);                                 % N_p x 6, row n = Ihat_LS(pbar_n)'
e = zeros(9*N_p, 1);  Jac = zeros(9*N_p, 60);
for n = 1:N_p
    In = I_LS(n,:).';
    for k = 1:3
        QI = Q(:,:,k,n) * In;
        e(3*(n-1)+k) = In.' * QI - r(k);
        Jac(3*(n-1)+k, :) = 2 * kron(QI.', T(n,:));          % 2 [Q_k Ihat_LS]_i tau(pbar_n)'
    end
    rows = 3*N_p + 6*(n-1) + (1:6);
    e(rows) = sqrt(LAM) * In;
    Jac(rows, :) = sqrt(LAM) * kron(eye(6), T(n,:));          % sqrt(LAM) tau(pbar_n)'
end
end
