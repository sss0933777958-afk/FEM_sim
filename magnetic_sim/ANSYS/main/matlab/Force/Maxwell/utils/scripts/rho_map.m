% rho_map.m -- minimum-norm unit current over all directions at one position (actuator frame).
% Problem (reference/Force_model/Quantitative measure, steps 2-3), for every direction r_hat(phi, theta):
%   I_opt_hat(phi, theta, p_i) = argmin ||I_hat||^2   s.t.   I_hat' * Q_k * I_hat = r_k,  k = 1, 2, 3
%   Q_k = F H_I' * L_k(p_i / l_hat, e_hat) * F H_I
%   rho(p_i, r_hat) = ||f_hat|| / ||I_opt_hat||^2 = 1 / ||I_opt_hat||^2   [pN/A^2]
% Solved with fmincon (sqp, analytic gradients).  Scan order: phi = -90 .. 90 (one latitude ring per phi), theta =
% 0 .. 360 inside a ring, starting at (phi, theta) = (-90, 0) from I_INIT.  Multistart at every grid point: the
% previous grid point's solution (chain) plus NS random starts; among the feasible results (exitflag > 0 and
% ||f_hat - r_hat|| < TOLC) the smallest ||I_hat||^2 is kept.  phi = +-90 is a single direction: solved once, copied.
% Global-optimality certificate: Lagrangian L = ||I_hat||^2 + sum_k lambda_k (I_hat' Q_k I_hat - r_k) (same sign
% convention as fmincon, lambda = lambda.eqnonlin); if E + sum_k lambda_k Q_k is positive semidefinite then for every
% feasible I
%   ||I||^2 = I' (E + sum lambda_k Q_k) I - lambda' r_hat >= -lambda' r_hat = ||I_opt_hat||^2,
% i.e. the solution is the global minimum (sufficient condition).  MINEIG = smallest eigenvalue; CERT = MINEIG >= -TOLE.
% Check: f_hat = [I_opt_hat' Q_k I_opt_hat]_k is compared with r_hat (err = ||f_hat - r_hat||).
% Also: capacity = ((1/4pi) * int int rho^3 cos(phi) dtheta dphi)^(1/3), isotropy = min rho / max rho.
% Params: P_UM [um, actuator frame] (default origin), DA [deg], NS (random starts per point),
%         TAG = 'eighteen' (N55 eighteen calibration).
% Output: utils/data/rho_map[_<TAG>].mat
if ~exist('P_UM','var'), P_UM = [0 0 0]; end
if ~exist('DA','var'), DA = 5; end
if ~exist('NS','var'), NS = 20; end
if ~exist('TAG','var'), TAG = 'eighteen'; end
clearvars -except P_UM DA NS TAG;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(here));   % .../Force/Maxwell
addpath(fullfile(CALROOT, 'function'));
MODEL = 'long2016_hexapole_halfcut';  GEOM = 'tip40um';
MN = struct('eighteen','calib_current_maxwell_axshN55_R150_eighteen');
MATNAME = MN.(TAG);  SFX = '';  if ~strcmp(TAG, 'eighteen'), SFX = ['_' TAG]; end
mgB = 0.0451;  UF = 1e3;                                    % M g_B [A um^2/mT]; unit factor -> pN
I_INIT = [0 0 0 0 0 1].';  TOLC = 1e-6;  TOLE = 1e-8;       % initial guess; constraint / PSD tolerances
rng(0);                                                     % reproducible random starts

cal = load_flux_calib(MODEL, GEOM, MATNAME, 'current');
gF  = UF * mgB / (2*cal.l_m) * cal.gB_m^2;                  % F g_I [pN/A^2]
FH  = sqrt(gF) * cal.Mbar_m;                                % F H_I [sqrt(pN)/A]
L   = build_L(P_UM(:).', cal.l_m, cal.e_m, cal.Pc_base);    % 1 x 6 x 6 x 3
Q   = zeros(6,6,3);
for k = 1:3, Lk = squeeze(L(1,:,:,k));  Q(:,:,k) = FH.' * Lk * FH;  Q(:,:,k) = (Q(:,:,k) + Q(:,:,k).')/2; end

ph = -90:DA:90;  th = 0:DA:360;  NPH = numel(ph);  NTH = numel(th);
IOPT = zeros(6, NPH, NTH);  FHAT = zeros(3, NPH, NTH);  LAM = zeros(3, NPH, NTH);  RHO = nan(NPH, NTH);
ERR  = nan(NPH, NTH);  MINEIG = nan(NPH, NTH);  NFEAS = zeros(NPH, NTH);
opt = optimoptions('fmincon', 'Algorithm','sqp', 'SpecifyObjectiveGradient',true, ...
                   'SpecifyConstraintGradient',true, 'OptimalityTolerance',1e-10, 'ConstraintTolerance',1e-12, ...
                   'StepTolerance',1e-14, 'MaxIterations',1000, 'Display','off');
x0 = I_INIT;  t0 = tic;
for a = 1:NPH
    for b = 1:NTH
        if abs(ph(a)) == 90 && b > 1                         % pole: one direction
            IOPT(:,a,b) = IOPT(:,a,1);  FHAT(:,a,b) = FHAT(:,a,1);  LAM(:,a,b) = LAM(:,a,1);  RHO(a,b) = RHO(a,1);
            ERR(a,b) = ERR(a,1);  MINEIG(a,b) = MINEIG(a,1);  NFEAS(a,b) = NFEAS(a,1);  continue;
        end
        r  = [cosd(ph(a))*cosd(th(b)); cosd(ph(a))*sind(th(b)); sind(ph(a))];
        R0 = randn(6, NS);  X0 = [x0, 0.35 * R0 ./ vecnorm(R0)];   % chain start + NS random starts
        best = inf;
        for s = 1:size(X0, 2)
            [I, ~, ef, ~, lm] = fmincon(@(x) obj_(x), X0(:,s), [],[],[],[],[],[], @(x) con_(x, Q, r), opt);
            f = [I.'*Q(:,:,1)*I; I.'*Q(:,:,2)*I; I.'*Q(:,:,3)*I];
            if ef <= 0 || norm(f - r) >= TOLC, continue; end
            NFEAS(a,b) = NFEAS(a,b) + 1;
            if I.'*I < best
                best = I.'*I;  lam = lm.eqnonlin;
                IOPT(:,a,b) = I;  FHAT(:,a,b) = f;  LAM(:,a,b) = lam;  RHO(a,b) = 1/best;  ERR(a,b) = norm(f - r);
                MINEIG(a,b) = min(eig(eye(6) + (lam(1)*Q(:,:,1) + lam(2)*Q(:,:,2) + lam(3)*Q(:,:,3))));
            end
        end
        if isfinite(best), x0 = IOPT(:,a,b); end           % best solution seeds the next grid point
    end
end

% capacity (trapezoid in theta and phi, solid-angle weight cos(phi)) and isotropy
CERT = MINEIG >= -TOLE;                                    % certified global minimum
nofs = ~isfinite(RHO);                                      % no feasible start
cap  = (trapz(deg2rad(ph), trapz(deg2rad(th), RHO.^3, 2) .* cosd(ph(:))) / (4*pi))^(1/3);
iso  = min(RHO(:)) / max(RHO(:));
seam = abs(RHO(:,end) - RHO(:,1)) ./ RHO(:,1);              % theta = 360 vs 0, same direction

U = true(NPH, NTH);  U(:,end) = false;  U(abs(ph) == 90, 2:end) = false;   % unique directions
fprintf('p = [%g %g %g] um, %d x %d grid (%d directions), 1 + %d starts each, %.0f s\n', ...
        P_UM, NPH, NTH, nnz(U), NS, toc(t0));
fprintf('rho: %.4f .. %.4f pN/A^2 | capacity %.4f pN/A^2 | isotropy %.4f\n', min(RHO(:)), max(RHO(:)), cap, iso);
fprintf('no feasible start: %d | max ||f_hat - r_hat|| %.2e | feasible starts per point: min %d, median %g\n', ...
        nnz(nofs & U), max(ERR(:)), min(NFEAS(U)), median(NFEAS(U)));
fprintf('PSD certificate: %d / %d directions certified | min eig(E + sum lambda_k Q_k): %.2e\n', ...
        nnz(CERT & U), nnz(U), min(MINEIG(:)));
fprintf('seam |rho(360) - rho(0)| / rho(0): max %.2e\n', max(seam));
R_act = cal.R_act;  l_hat = cal.l_m;
save(fullfile(CALROOT, 'utils', 'data', ['rho_map' SFX '.mat']), 'P_UM', 'DA', 'TAG', 'ph', 'th', 'IOPT', 'FHAT', ...
     'LAM', 'RHO', 'ERR', 'MINEIG', 'CERT', 'NFEAS', 'NS', 'cap', 'iso', 'I_INIT', 'Q', 'R_act', 'l_hat', 'gF', 'mgB', 'UF', 'MODEL', 'MATNAME');

% ---------------------------------------------------------------------------------------
function [v, g] = obj_(x)
% ||I_hat||^2 and its gradient
    v = x.'*x;  g = 2*x;
end

function [c, ceq, gc, gceq] = con_(x, Q, r)
% equality constraints I_hat' Q_k I_hat - r_k = 0 (k = 1..3) and their gradients (6 x 3)
    c = [];  gc = [];
    ceq  = [x.'*Q(:,:,1)*x - r(1); x.'*Q(:,:,2)*x - r(2); x.'*Q(:,:,3)*x - r(3)];
    gceq = 2*[Q(:,:,1)*x, Q(:,:,2)*x, Q(:,:,3)*x];
end
