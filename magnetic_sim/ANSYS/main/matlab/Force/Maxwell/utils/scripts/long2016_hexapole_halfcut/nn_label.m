%% nn_label.m -- training labels of the NN inverse model: I_opt_hat at every grid position and direction.
%  Positions: nn_grid.mat (P_act, actuator frame, 30 um cubic grid, |p| <= 150 um, 515 points).
%  Directions: rho_solve grid, phi = -90:DA:90, theta = 0:DA:360 (actuator frame).
%  Solver: rho_solve.m stage 1 (fmincon multistart = numerical Lagrange conditions + PSD certificate), the same
%  problem as rho_map.m:  I_opt_hat = argmin ||I_hat||^2  s.t.  I_hat' Q_k I_hat = r_k,  k = 1..3,
%  Q_k = F H_I' L_k(p / l_hat, e_hat) F H_I   (N55 eighteen calibration, actuator frame).
%  Everything is in the ACTUATOR frame: positions, directions r_hat, and I_opt_hat (pole order P1..P6).
%  Gurobi spot check (as cap_map.m): every GEV-th position (i = 1, 1+GEV, ...) also runs stage 2 on GMAX random
%  directions without the PSD certificate (Gurobi, TimeLimit TLIM, one thread); a smaller polished Gurobi
%  solution replaces the fmincon one.  Per checked position: NGRB directions sent, NOPT proved OPTIMAL,
%  NIMP improved (should be 0).
%
%  Workers: the positions are dealt round-robin to NW workers; worker WID solves positions WID, WID+NW, ...
%  and appends each finished position to its own file, so a stopped run resumes where it stopped.
%      one process:   WID = 1; NW = 1; run nn_label
%      NW processes:  matlab -batch "WID=w; NW=12; run('<this file>')" -singleCompThread   (w = 1..12)
%  NW must stay the same across a resumed run (it fixes which worker owns which position).
%  MERGE = true: collect the worker files into nn_label.mat (needs every position done) and delete them.
%
%  Output: utils/data/long2016_hexapole_halfcut/nn_label.mat
%    P_act   NP x 3                [um]   positions (actuator frame)
%    ph, th  1 x NPH, 1 x NTH      [deg]  direction grid
%    IOPT    6 x NPH x NTH x NP    [A]    I_opt_hat, unit force 1 pN along r_hat(ph, th)
%    RHO     NPH x NTH x NP        [pN/A^2]  1 / ||I_opt_hat||^2
%    ERR     NPH x NTH x NP        ||f_hat - r_hat||
%    MINEIG, CERT  NPH x NTH x NP  PSD certificate (smallest eigenvalue, certified)
%    NFEAS   NPH x NTH x NP        feasible starts
%    T1      NP x 1                [s]    stage-1 time per position
%    GCHK, NGRB, NOPT, NIMP  NP x 1        Gurobi spot check (checked?, sent, OPTIMAL, improved)
%    phi = +-90 is one direction, stored at every theta (copied by rho_solve).
%  Worker files (deleted by MERGE): utils/data/long2016_hexapole_halfcut/nn_label_w<WID>of<NW>.mat
if ~exist('WID','var'), WID = 1; end
if ~exist('NW','var'), NW = 1; end
if ~exist('DA','var'), DA = 5; end
if ~exist('NS','var'), NS = 40; end
if ~exist('MERGE','var'), MERGE = false; end
if ~exist('GEV','var'), GEV = 6; end                           % Gurobi check on every GEV-th position
if ~exist('GMAX','var'), GMAX = 60; end                        % directions per checked position
if ~exist('TLIM','var'), TLIM = 600; end                       % Gurobi TimeLimit per direction [s]
if ~exist('NMAX','var'), NMAX = Inf; end                       % stop after NMAX new positions (test runs)
clearvars -except WID NW DA NS MERGE NMAX GEV GMAX TLIM;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
addpath(fullfile(CALROOT, 'function'), here);                  % here: rho_solve.m
MODEL = 'long2016_hexapole_halfcut';  GEOM = 'tip40um';
MATNAME = 'calib_current_maxwell_axshN55_R150_eighteen';
DAT = fullfile(CALROOT, 'utils', 'data', MODEL);
G   = load(fullfile(DAT, 'nn_grid.mat'), 'P_act');
P_act = G.P_act;  NP = size(P_act, 1);
wfile = @(w) fullfile(DAT, sprintf('nn_label_w%02dof%02d.mat', w, NW));

if MERGE
    for w = 1:NW
        W = load(wfile(w));
        if w == 1
            ph = W.ph;  th = W.th;  NPH = numel(ph);  NTH = numel(th);
            IOPT = nan(6, NPH, NTH, NP);  RHO = nan(NPH, NTH, NP);  ERR = RHO;  MINEIG = RHO;
            CERT = false(NPH, NTH, NP);  NFEAS = zeros(NPH, NTH, NP);  T1 = nan(NP, 1);  DONE = false(NP, 1);
            GCHK = false(NP, 1);  NGRB = zeros(NP, 1);  NOPT = zeros(NP, 1);  NIMP = zeros(NP, 1);
        end
        k = W.idx;
        IOPT(:,:,:,k) = W.IOPT;  RHO(:,:,k) = W.RHO;  ERR(:,:,k) = W.ERR;  MINEIG(:,:,k) = W.MINEIG;
        CERT(:,:,k) = W.CERT;  NFEAS(:,:,k) = W.NFEAS;  T1(k) = W.T1;  DONE(k) = true;
        GCHK(k) = W.GCHK;  NGRB(k) = W.NGRB;  NOPT(k) = W.NOPT;  NIMP(k) = W.NIMP;
    end
    assert(all(DONE), 'nn_label:missing', '%d of %d positions not done yet', nnz(~DONE), NP);
    l_hat = W.l_hat;  gF = W.gF;
    save(fullfile(DAT, 'nn_label.mat'), 'P_act', 'ph', 'th', 'IOPT', 'RHO', 'ERR', 'MINEIG', 'CERT', 'NFEAS', 'T1', ...
         'GCHK', 'NGRB', 'NOPT', 'NIMP', 'GEV', 'GMAX', 'TLIM', 'DA', 'NS', 'MODEL', 'MATNAME', 'l_hat', 'gF', '-v7.3');
    for w = 1:NW, delete(wfile(w)); end
    fprintf('merged %d positions -> %s | certified %.1f %% | max err %.2g | stage-1 time %.1f h total\n', ...
            NP, fullfile(DAT, 'nn_label.mat'), 100*mean(CERT(:)), max(ERR(:)), sum(T1)/3600);
    return
end

% ---- force quadratic forms (same construction as rho_map.m) ----
mgB = 0.0451;  UF = 1e3;                                       % M g_B [A um^2/mT]; unit factor -> pN
cal = load_flux_calib(MODEL, GEOM, MATNAME, 'current');
gF  = UF * mgB / (2*cal.l_m) * cal.gB_m^2;                     % F g_I [pN/A^2]
FH  = sqrt(gF) * cal.Mbar_m;                                   % F H_I [sqrt(pN)/A]
l_hat = cal.l_m;
L   = build_L(P_act, cal.l_m, cal.e_m, cal.Pc_base);           % NP x 6 x 6 x 3, actuator frame

mine = WID:NW:NP;                                              % positions owned by this worker
if isfile(wfile(WID))
    W = load(wfile(WID));
else
    W = struct('idx', zeros(1,0), 'IOPT', [], 'RHO', [], 'ERR', [], 'MINEIG', [], 'CERT', [], 'NFEAS', [], 'T1', zeros(1,0), ...
               'GCHK', false(1,0), 'NGRB', zeros(1,0), 'NOPT', zeros(1,0), 'NIMP', zeros(1,0));
end
todo = setdiff(mine, W.idx, 'stable');
fprintf('worker %d / %d: %d positions, %d done, %d to go\n', WID, NW, numel(mine), numel(W.idx), numel(todo));
o0 = struct('DA',DA, 'NS',NS, 'GRB',false, 'SEED',0, 'RETRY',true, 'GTHR',1, 'TLIM',TLIM, 'GMAX',GMAX);
nnew = 0;
for i = todo
    if nnew >= NMAX, break; end
    Q = zeros(6,6,3);
    for k = 1:3, Lk = squeeze(L(i,:,:,k));  Q(:,:,k) = FH.' * Lk * FH;  Q(:,:,k) = (Q(:,:,k) + Q(:,:,k).')/2; end
    o = o0;  o.GRB = mod(i - 1, GEV) == 0;                     % Gurobi spot check on this position?
    S = rho_solve(Q, o);
    W.idx(end+1) = i;  W.T1(end+1) = S.t1;
    W.GCHK(end+1) = o.GRB;  W.NGRB(end+1) = nnz(S.GRBD);  W.NOPT(end+1) = nnz(S.GRBD & S.GSTAT == "OPTIMAL");
    W.NIMP(end+1) = nnz(S.RHO ~= S.RHO0);
    W.IOPT   = cat(4, W.IOPT, S.IOPT);    W.RHO   = cat(3, W.RHO, S.RHO);    W.ERR  = cat(3, W.ERR, S.ERR);
    W.MINEIG = cat(3, W.MINEIG, S.MINEIG);  W.CERT = cat(3, W.CERT, S.CERT);  W.NFEAS = cat(3, W.NFEAS, S.NFEAS);
    W.ph = S.ph;  W.th = S.th;  W.l_hat = l_hat;  W.gF = gF;  W.NW = NW;  W.WID = WID;
    save(wfile(WID), '-struct', 'W');
    nnew = nnew + 1;
    fprintf('[w%02d] %3d / %d  p = [%g %g %g] um  rho %.3f..%.3f pN/A^2  cert %.1f %%  %.0f s\n', WID, numel(W.idx), ...
            numel(mine), P_act(i,:), min(S.RHO(:)), max(S.RHO(:)), 100*mean(S.CERT(:)), S.t1);
end
