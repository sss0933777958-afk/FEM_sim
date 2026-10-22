%% ls_fit.m -- least-squares inverse model (Long et al., 2021), symbols as in C:\Users\Kuo\Downloads\LS_inverse_model.pdf.
%  For every direction (phi_j, theta_j) and coil i the unit-force optimal current is fitted by a quadratic form
%      Ihat_opt(phi_j, theta_j, pbar)_i ~ P' D_LS_i(phi_j, theta_j) P,   P = [xbar ybar zbar 1]',  pbar = p / l_hat   (Sec. 3)
%  linear in the 10 coefficients d_i = [d44 d14 d24 d34 d12 d13 d23 d11 d22 d33]':
%      tau(pbar) = [1, 2xbar, 2ybar, 2zbar, 2xbar*ybar, 2xbar*zbar, 2ybar*zbar, xbar^2, ybar^2, zbar^2]'
%  Stacking the N_p positions (Sec. 4):  T = [tau(pbar_1)'; ...; tau(pbar_N_p)']  (N_p x 10),  y_i = Ihat_opt at them (N_p x 1)
%  Normal equations (Sec. 5):  d_i = (T'T)^-1 T' y_i;  (T'T)^-1 T' is shared by every direction and coil (Sec. 6).
%  Labels y_i: nn_label_sgn.mat (Gurobi-certified I_opt_hat, sign-aligned so that it is continuous in pbar).
%  Actuator frame throughout; coil i = pole P1..P6.
%  Directions j = 1..N_j: the distinct directions of the label grid (theta = 360 dropped, phi = +-90 kept once at theta = 0),
%  ordered phi fastest, then theta.
%  REG selects the fitting positions:  'cube' = |x|, |y|, |z| <= A (default A = 45 um: the 90 um cube of Long's
%                                               dissertation p. 49, all octants)
%                                      'all'  = every grid node, |p| <= 150 um (N_p = 515)
%                                      'oct'  = first octant x, y, z >= 0
%  Output: utils/data/long2016_hexapole_halfcut/ls_fit[_all|_oct].mat
%    D_LS            4 x 4 x N_j x N_I   [A]    D_LS_i(phi_j, theta_j), entry (:,:,j,i)
%    d               10 x N_j x N_I      [A]    the same coefficients as d_i
%    phi_j, theta_j  N_j x 1             [deg]  directions
%    pbar            N_p x 3                    fitting positions (normalised)
%    l_hat           [um]
%    e_ij            N_j x N_p           [%]    relative force error at position pbar_i, direction j:
%                    e_ij = || fhat_LS(pbar_i, theta_j, phi_j) - fhat_d(pbar_i, theta_j, phi_j) || / || fhat_d || * 100 %
%                    fhat_d = r_hat(phi_j, theta_j) (unit force, 1 pN); fhat_LS_k = Ihat_LS' Q_k Ihat_LS with
%                    Ihat_LS = P' D_LS P (all six coils) and Q_k = F H_I' L_k(pbar) F H_I as in nn_label.m
if ~exist('REG','var'), REG = 'cube'; end
if ~exist('A','var'), A = 45; end
clearvars -except REG A;
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
xbar = pbar(:,1);  ybar = pbar(:,2);  zbar = pbar(:,3);

% ---- distinct directions j ----
[PH, TH] = ndgrid(S.ph, S.th);
U = TH < 360 & (abs(PH) < 90 | TH == 0);
phi_j = PH(U);  theta_j = TH(U);  N_j = numel(phi_j);

% ---- Sec. 3-4: T (N_p x 10) ----
T = [ones(N_p,1), 2*xbar, 2*ybar, 2*zbar, 2*xbar.*ybar, 2*xbar.*zbar, 2*ybar.*zbar, xbar.^2, ybar.^2, zbar.^2];
assert(N_p >= 10 && rank(T) == 10, 'ls_fit:rank', 'T must have full column rank (N_p = %d, rank %d)', N_p, rank(T));

% ---- Sec. 5: d_i = (T'T)^-1 T' y_i for every direction and coil at once ----
I = reshape(S.IOPT(:,:,:,n), N_I, N_phi*N_theta, N_p);
y = reshape(permute(I(:, U(:), :), [3 2 1]), N_p, []);       % column = y_i(phi_j, theta_j), order (j, i)
d = (T.'*T) \ (T.'*y);                                        % 10 x (N_j N_I)
Ihat_LS = reshape(T*d, N_p, N_j, N_I);                      % P' D_LS P at the fitting positions

% ---- Sec. 3: fill D_LS_i from d_i = [d44 d14 d24 d34 d12 d13 d23 d11 d22 d33] ----
d44 = d(1,:);  d14 = d(2,:);  d24 = d(3,:);  d34 = d(4,:);  d12 = d(5,:);
d13 = d(6,:);  d23 = d(7,:);  d11 = d(8,:);  d22 = d(9,:);  d33 = d(10,:);
D_LS = reshape([d11; d12; d13; d14;  d12; d22; d23; d24;  d13; d23; d33; d34;  d14; d24; d34; d44], 4, 4, N_j, N_I);
d    = reshape(d, 10, N_j, N_I);

% ---- relative force error e_ij (force model as nn_label.m) ----
mgB = 0.0451;  UF = 1e3;                                       % M g_B [A um^2/mT]; unit factor -> pN
cal = load_flux_calib(MODEL, GEOM, S.MATNAME, 'current');
assert(abs(cal.l_m - l_hat) < 1e-9, 'ls_fit:lhat', 'l_hat differs from the label run');
gF  = UF * mgB / (2*cal.l_m) * cal.gB_m^2;                     % F g_I [pN/A^2]
FH  = sqrt(gF) * cal.Mbar_m;                                   % F H_I [sqrt(pN)/A]
L   = build_L(pbar * l_hat, cal.l_m, cal.e_m, cal.Pc_base);    % N_p x 6 x 6 x 3
fhat_d = [cosd(phi_j).*cosd(theta_j), cosd(phi_j).*sind(theta_j), sind(phi_j)];   % N_j x 3
e_ij = zeros(N_j, N_p);
for p = 1:N_p
    X = squeeze(Ihat_LS(p, :, :));                             % N_j x N_I
    fhat_LS = zeros(N_j, 3);
    for k = 1:3
        Q_k = FH.' * squeeze(L(p,:,:,k)) * FH;  Q_k = (Q_k + Q_k.')/2;
        fhat_LS(:,k) = sum((X*Q_k) .* X, 2);
    end
    e_ij(:, p) = vecnorm(fhat_LS - fhat_d, 2, 2) ./ vecnorm(fhat_d, 2, 2) * 100;
end

% self-check: P' D_LS P reproduces the label column of one direction and coil through T d
Pn = [pbar, ones(N_p,1)];  jc = 7;  ic = 4;
assert(max(abs(sum((Pn*D_LS(:,:,jc,ic)).*Pn, 2) - T*d(:,jc,ic))) < 1e-12, 'ls_fit:D', 'D_LS fill does not match d');
a = find(S.ph == phi_j(jc));  b = find(S.th == theta_j(jc));
assert(isequal(squeeze(S.IOPT(ic, a, b, n)), y(:, (ic-1)*N_j + jc)), 'ls_fit:y', 'label column order mismatch');

save(fullfile(DAT, OUT), 'D_LS', 'd', 'phi_j', 'theta_j', 'pbar', 'l_hat', 'e_ij', 'REG', 'A', 'N_p', 'N_j', 'N_I');
fprintf('REG %s: N_p = %d, N_j = %d, N_I = %d, cond(T''T) = %.3g -> %s\n', REG, N_p, N_j, N_I, cond(T.'*T), fullfile(DAT, OUT));
fprintf('e_ij [%%]: median %.3g, 95th pct %.3g, max %.3g\n', median(e_ij(:)), prctile(e_ij(:), 95), max(e_ij(:)));
