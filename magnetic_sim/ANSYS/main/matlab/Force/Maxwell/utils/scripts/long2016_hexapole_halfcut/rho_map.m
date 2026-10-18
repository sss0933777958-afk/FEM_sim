% rho_map.m -- minimum-norm unit current over all directions at one position (actuator frame).  Solver: rho_solve.m.
% Problem (reference/Force_model/Quantitative measure, steps 2-3), for every direction r_hat(phi, theta):
%   I_opt_hat(phi, theta, p_i) = argmin ||I_hat||^2   s.t.   I_hat' * Q_k * I_hat = r_k,  k = 1, 2, 3
%   Q_k = F H_I' * L_k(p_i / l_hat, e_hat) * F H_I
%   rho(p_i, r_hat) = ||f_hat|| / ||I_opt_hat||^2 = 1 / ||I_opt_hat||^2   [pN/A^2]
% Solved with fmincon (sqp, analytic gradients).  Scan order: phi = -90 .. 90 (one latitude ring per phi), theta =
% 0 .. 360 inside a ring, starting at (phi, theta) = (-90, 0) from I_INIT.  Multistart at every grid point: the
% previous grid point's solution (chain) plus NS random starts; among the feasible results (exitflag > 0 and
% ||f_hat - r_hat|| < TOLC) the smallest ||I_hat||^2 is kept.  phi = +-90 is a single direction: solved once, copied.
% Global-optimality certificate: Lagrangian L = ||I_hat||^2 - sum_k lambda_k (I_hat' Q_k I_hat - r_k) (document
% convention; fmincon uses +, so lambda = -lambda.eqnonlin); if E - sum_k lambda_k Q_k is positive semidefinite then
% for every feasible I
%   ||I||^2 = I' (E - sum lambda_k Q_k) I + lambda' r_hat >= lambda' r_hat = ||I_opt_hat||^2,
% i.e. the solution is the global minimum (sufficient condition).  MINEIG = smallest eigenvalue; CERT = MINEIG >= -TOLE.
% Stage 2 (GRB = true): every direction without the PSD certificate is solved once more with Gurobi (nonconvex
% QCQP, spatial branch and bound, NonConvex = 2), box |I_j| <= ||I_fmincon|| plus ball ||I||^2 <= ||I_fmincon||^2
% (cannot cut off the global minimum), start = fmincon solution.  The Gurobi point is polished by one local fmincon
% run (accurate multipliers), then lambda and MINEIG are recomputed.  The smaller ||I_hat||^2 of fmincon / polished
% Gurobi is kept; if the polish fails (PLSH false) the fmincon solution stays.  GLOB = CERT | (Gurobi status OPTIMAL, gap within MIPGap / MIPGapAbs): PSD is sufficient only, a duality gap leaves it negative
% even at the global minimum, so Gurobi's own bound is the certificate there.  Stage-1 results kept as *0 fields.
% Check: f_hat = [I_opt_hat' Q_k I_opt_hat]_k is compared with r_hat (err = ||f_hat - r_hat||); KKT = residual of
% I_hat = sum_k lambda_k Q_k I_hat.
% Also: capacity = ((1/4pi) * int int rho^3 cos(phi) dtheta dphi)^(1/3), isotropy = min rho / max rho.
% Params: P_UM [um, actuator frame] (default origin), DA [deg], NS (random starts per point),
%         TAG = 'eighteen' (N55 eighteen calibration), GRB (stage 2 on / off),
%         RELOAD (skip stage 1, take it from the saved output .mat), PTAG = position tag ('' = origin),
%         GTHR = Gurobi Threads (0 = Gurobi default).
% Output: utils/data/long2016_hexapole_halfcut/rho_map[_<TAG>][_<PTAG>].mat
if ~exist('P_UM','var'), P_UM = [0 0 0]; end
if ~exist('DA','var'), DA = 5; end
if ~exist('NS','var'), NS = 20; end
if ~exist('TAG','var'), TAG = 'eighteen'; end
if ~exist('GRB','var'), GRB = true; end
if ~exist('RELOAD','var'), RELOAD = false; end
if ~exist('PTAG','var'), PTAG = ''; end
if ~exist('GTHR','var'), GTHR = 0; end                  % Gurobi Threads (0 = all cores; 1 when many runs share the machine)
clearvars -except P_UM DA NS TAG GRB RELOAD PTAG GTHR;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
addpath(fullfile(CALROOT, 'function'), here);              % here: rho_solve.m (the solver, shared with cap_map.m)
MODEL = 'long2016_hexapole_halfcut';  GEOM = 'tip40um';
MN = struct('eighteen','calib_current_maxwell_axshN55_R150_eighteen');
MATNAME = MN.(TAG);  SFX = '';  if ~strcmp(TAG, 'eighteen'), SFX = ['_' TAG]; end
if ~isempty(PTAG), SFX = [SFX '_' PTAG]; end
mgB = 0.0451;  UF = 1e3;                                    % M g_B [A um^2/mT]; unit factor -> pN
I_INIT = [0 0 0 0 0 1].';                                   % first start of the first direction
OUT = fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut', ['rho_map' SFX '.mat']);
o = struct('DA',DA, 'NS',NS, 'GRB',GRB, 'GTHR',GTHR, 'TLIM',60, 'SEED',0, 'I_INIT',I_INIT);   % solver: rho_solve.m

if RELOAD
    S = load(OUT);                                          % stage 1 from file (the *0 copies if stage 2 ran)
    P_UM = S.P_UM;  DA = S.DA;  NS = S.NS;  Q = S.Q;  R_act = S.R_act;  l_hat = S.l_hat;  gF = S.gF;
    if isfield(S, 'IOPT0'), S0 = '0'; else, S0 = ''; end
    o.DA = DA;  o.NS = NS;
    o.S1 = struct('IOPT',S.(['IOPT' S0]), 'FHAT',S.(['FHAT' S0]), 'LAM',S.(['LAM' S0]), 'RHO',S.(['RHO' S0]), ...
                  'ERR',S.(['ERR' S0]), 'MINEIG',S.(['MINEIG' S0]), 'NFEAS',S.NFEAS);
    fprintf('stage 1 reloaded from %s\n', OUT);
else
    cal = load_flux_calib(MODEL, GEOM, MATNAME, 'current');
    gF  = UF * mgB / (2*cal.l_m) * cal.gB_m^2;              % F g_I [pN/A^2]
    FH  = sqrt(gF) * cal.Mbar_m;                            % F H_I [sqrt(pN)/A]
    L   = build_L(P_UM(:).', cal.l_m, cal.e_m, cal.Pc_base);   % 1 x 6 x 6 x 3
    Q   = zeros(6,6,3);
    for k = 1:3, Lk = squeeze(L(1,:,:,k));  Q(:,:,k) = FH.' * Lk * FH;  Q(:,:,k) = (Q(:,:,k) + Q(:,:,k).')/2; end
    R_act = cal.R_act;  l_hat = cal.l_m;
end
RS = rho_solve(Q, o);

% ---- report ----
NPH = numel(RS.ph);  NTH = numel(RS.th);
U = true(NPH, NTH);  U(:,end) = false;  U(abs(RS.ph) == 90, 2:end) = false;   % unique directions
seam = abs(RS.RHO(:,end) - RS.RHO(:,1)) ./ RS.RHO(:,1);    % theta = 360 vs 0, same direction
fprintf('p = [%g %g %g] um, %d x %d grid (%d directions), 1 + %d starts each, stage 1 %.0f s, stage 2 %.0f s\n', ...
        P_UM, NPH, NTH, nnz(U), NS, RS.t1, RS.t2);
fprintf('rho: %.4f .. %.4f pN/A^2 | capacity %.4f pN/A^2 | isotropy %.4f\n', min(RS.RHO(:)), max(RS.RHO(:)), RS.cap, RS.iso);
fprintf('no feasible start: %d | max ||f_hat - r_hat|| %.2e | feasible starts per point: min %d, median %g\n', ...
        nnz(~isfinite(RS.RHO) & U), max(RS.ERR(:)), min(RS.NFEAS(U)), median(RS.NFEAS(U)));
fprintf('PSD certificate: %d / %d directions certified | min eig(E - sum lambda_k Q_k): %.2e\n', ...
        nnz(RS.CERT & U), nnz(U), min(RS.MINEIG(:)));
fprintf('seam |rho(360) - rho(0)| / rho(0): max %.2e | max KKT residual %.2e\n', max(seam), max(RS.KKT(:)));
if GRB
    W = RS.GRBD & U;
    fprintf('Gurobi: %d directions | status %s | max gap %.2e | ||I||^2 fmincon - Gurobi: %.2e .. %.2e\n', ...
            nnz(W), strjoin(unique(RS.GSTAT(W)).', '/'), max(RS.GGAP(W)), min(RS.DOBJ(W)), max(RS.DOBJ(W)));
    fprintf('        polish failed: %d | replaced by Gurobi: %d | PSD after Gurobi: %d / %d | global (PSD or Gurobi): %d / %d\n', ...
            nnz(W & ~RS.PLSH), nnz(W & RS.RHO ~= RS.RHO0), nnz(W & RS.CERT), nnz(W), nnz(RS.GLOB & U), nnz(U));
end
RS = rmfield(RS, {'t1', 't2'});
save(OUT, '-struct', 'RS');
save(OUT, 'P_UM', 'DA', 'TAG', 'NS', 'I_INIT', 'Q', 'R_act', 'l_hat', 'gF', 'mgB', 'UF', 'MODEL', 'MATNAME', 'PTAG', '-append');
