% cap_map.m -- capacity C_f(p) along rays from the origin, for the C_f iso-surfaces (actuator frame).
% Rays: phi = -90:DRAY:90, theta = 0:DRAY:360-DRAY over the whole sphere (phi = +-90 once each) -> ray ids 1..NR.
% Along a ray (unit vector u, position p = r u [um]):
%   coarse  r = RSTEP, 2 RSTEP, ... RMAX, stop once C_f > max(LEV) (points already on the ray are reused);
%           r = 0 uses the origin value of rho_map.mat
%   bisect  for every level L in LEV: halve the bracket [r_a, r_b] with C_f(r_a) < L <= C_f(r_b) NBIS times, then
%           the crossing RL is linear interpolation inside the last bracket (RSTEP = 50, NBIS = 4 -> 3.1 um bracket)
% Each point: rho_solve.m stage 1 only (fmincon 1 + NS random starts + the solution of the nearest point already
% done on this ray as an extra start, Lagrange PSD certificate; a direction without a feasible start is retried with
% new random starts until feasible), direction grid DA.  Gurobi spot check: on every SPOT-th ray (id - 1 divisible by
% SPOT) the point with C_f closest to LEV(2) is re-run with stage 2 on a random GMAX of its uncertified directions
% (Threads GTHR, TimeLimit TLIM); recorded: change of C_f, how many of them Gurobi proved OPTIMAL.
% Parallel: NW workers (separate MATLAB processes, -singleCompThread), worker WID takes rays perm(WID:NW:end) of a
% fixed random permutation; each worker writes WDIR/cap_w<WID>.mat after every ray and skips rays already there
% (restartable; a ray stored for other levels LEV is continued from its points); its spot-check rays go first.
% MERGE = true: collect all worker files -> utils/data/long2016_hexapole_halfcut/cap_map.mat.
% Params: WID, NW, WDIR, MERGE, DRAY [deg], RSTEP, RMAX [um], LEV [pN/A^2], NBIS, DA [deg], NS, SPOT, GTHR, TLIM [s],
%         GMAX
if ~exist('WID','var'), WID = 1; end
if ~exist('NW','var'), NW = 1; end
if ~exist('WDIR','var'), WDIR = fullfile(tempdir, 'cap_map'); end
if ~exist('MERGE','var'), MERGE = false; end
if ~exist('DRAY','var'), DRAY = 15; end
if ~exist('RSTEP','var'), RSTEP = 50; end
if ~exist('RMAX','var'), RMAX = 450; end
if ~exist('LEV','var'), LEV = [14 15 16]; end
if ~exist('NBIS','var'), NBIS = 4; end
if ~exist('DA','var'), DA = 5; end
if ~exist('NS','var'), NS = 10; end
if ~exist('SPOT','var'), SPOT = 6; end
if ~exist('GTHR','var'), GTHR = 1; end
if ~exist('TLIM','var'), TLIM = 600; end
if ~exist('GMAX','var'), GMAX = 60; end
clearvars -except WID NW WDIR MERGE DRAY RSTEP RMAX LEV NBIS DA NS SPOT GTHR TLIM GMAX;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
addpath(fullfile(CALROOT, 'function'), here);              % here: rho_solve.m
MODEL = 'long2016_hexapole_halfcut';  GEOM = 'tip40um';  MATNAME = 'calib_current_maxwell_axshN55_R150_eighteen';
mgB = 0.0451;  UF = 1e3;                                    % M g_B [A um^2/mT]; unit factor -> pN
OUT = fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut', 'cap_map.mat');

% ---- rays (actuator frame) ----
[PH, TH] = ndgrid(-90:DRAY:90, 0:DRAY:360-DRAY);
keep = ~(abs(PH) == 90 & TH > 0);                           % pole: one ray
PH = PH(keep);  TH = TH(keep);  NR = numel(PH);
UR = [cosd(PH).*cosd(TH), cosd(PH).*sind(TH), sind(PH)];   % NR x 3 unit vectors
if ~exist(WDIR, 'dir'), mkdir(WDIR); end

if MERGE
    RAY = [];
    for f = dir(fullfile(WDIR, 'cap_w*.mat')).'
        W = load(fullfile(f.folder, f.name), 'RAY');  RAY = [RAY, W.RAY];
    end
    [~, ix] = sort([RAY.id]);  RAY = RAY(ix);
    assert(isequal([RAY.id], 1:NR), 'cap_map:merge', '%d of %d rays present', numel(RAY), NR);
    assert(all(arrayfun(@(x) isequal(x.lev, LEV), RAY)), 'cap_map:merge', 'some rays were solved for other levels');
    RL = vertcat(RAY.RL);                                   % NR x numel(LEV) crossing radii [um] (NaN = none)
    SP = [RAY.spot];                                        % spot checks
    fprintf('merged %d rays | crossings found: %s of %d | spot checks %d, |dC_f| max %.2e, all global %d\n', NR, ...
            mat2str(sum(isfinite(RL), 1)), NR, numel(SP), max(abs([SP.dcap])), all([SP.allglob]));
    S0 = load(fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut', 'rho_map.mat'), 'R_act');  R_act = S0.R_act;
    save(OUT, 'RAY', 'RL', 'PH', 'TH', 'UR', 'LEV', 'DRAY', 'RSTEP', 'RMAX', 'NBIS', 'DA', 'NS', 'SPOT', 'TLIM', 'GMAX', ...
         'R_act', 'MODEL', 'MATNAME');
    fprintf('wrote %s\n', OUT);
    return
end

% ---- worker ----
cal = load_flux_calib(MODEL, GEOM, MATNAME, 'current');
gF  = UF * mgB / (2*cal.l_m) * cal.gB_m^2;                  % F g_I [pN/A^2]
FH  = sqrt(gF) * cal.Mbar_m;                                % F H_I [sqrt(pN)/A]
S0  = load(fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut', 'rho_map.mat'), 'cap');   CAP0 = S0.cap;   % C_f at the origin
perm = randperm(RandStream('mt19937ar', 'Seed', 1), NR);
MINE = perm(WID:NW:end);
WF = fullfile(WDIR, sprintf('cap_w%02d.mat', WID));
RAY = struct('id',{}, 'ph',{}, 'th',{}, 'r',{}, 'cap',{}, 'iso',{}, 'rho',{}, 'ncert',{}, 'nfeas',{}, 't',{}, ...
             'RL',{}, 'spot',{}, 'lev',{});
if exist(WF, 'file')
    W = load(WF, 'RAY');  RAY = W.RAY;
    if ~isfield(RAY, 'lev'), [RAY.lev] = deal([15 20 25]); end   % rays done before the levels were changed
end
% a ray is done when it was solved for the current levels; a ray solved for other levels is continued from its points
done = arrayfun(@(x) isequal(x.lev, LEV), RAY);
todo = setdiff(MINE, [RAY(done).id], 'stable');
todo = [todo(mod(todo - 1, SPOT) == 0), todo(mod(todo - 1, SPOT) ~= 0)];   % spot-check rays first
fprintf('worker %d / %d: %d rays (%d done, %d to go, %d of them continued) | levels %s\n', WID, NW, numel(MINE), ...
        numel(MINE) - numel(todo), numel(todo), nnz(ismember(todo, [RAY.id])), mat2str(LEV));
o = struct('DA',DA, 'NS',NS, 'GRB',false, 'RETRY',true);
tw = tic;
for id = todo
    tr = tic;  u = UR(id,:);
    j = find([RAY.id] == id, 1);
    if isempty(j)                                           % new ray; r = 0: origin
        R = struct('r',0, 'cap',CAP0, 'iso',NaN, 'rho',[NaN NaN], 'ncert',NaN, 'nfeas',NaN, 't',0);  sp = [];
    else                                                    % continue a ray solved for other levels
        R = struct('r',RAY(j).r, 'cap',RAY(j).cap, 'iso',RAY(j).iso, 'rho',RAY(j).rho, 'ncert',RAY(j).ncert, ...
                   'nfeas',RAY(j).nfeas, 't',RAY(j).t);  sp = RAY(j).spot;
    end
    n0 = numel(R.r);
    P1 = cell(1, n0);                                       % stage-1 results per point (for the extra start / spot)
    % coarse scan (points already on the ray are kept)
    for r = RSTEP:RSTEP:RMAX
        if any(R.cap > LEV(end)), break; end
        if ~any(abs(R.r - r) < 1e-9), [R, P1] = eval_(R, P1, r, u, id, cal, FH, o); end
    end
    % bisection per level
    RL = nan(1, numel(LEV));
    for q = 1:numel(LEV)
        for it = 1:NBIS
            [ra, rb] = bracket_(R, LEV(q));
            if isnan(ra) || rb - ra <= RSTEP / 2^NBIS + 1e-9, break; end
            [R, P1] = eval_(R, P1, (ra + rb)/2, u, id, cal, FH, o);
        end
        [ra, rb, ca, cb] = bracket_(R, LEV(q));
        if ~isnan(ra), RL(q) = ra + (LEV(q) - ca) / (cb - ca) * (rb - ra); end
    end
    % Gurobi spot check on the point closest to LEV(2) (a continued ray keeps its earlier check)
    if mod(id - 1, SPOT) == 0 && isempty(sp)
        [~, k] = min(abs(R.cap(2:end) - LEV(2)));  k = k + 1;
        Q = q_at_(R.r(k) * u, cal, FH);
        o2 = struct('DA',DA, 'NS',NS, 'GRB',true, 'GTHR',GTHR, 'TLIM',TLIM, 'S1',P1{k}, 'GMAX',GMAX, 'SEED',id, ...
                    'RETRY',true);                          % S1 empty (point from the file) -> stage 1 redone
        ts = tic;  G = rho_solve(Q, o2);
        U = true(size(G.RHO));  U(:,end) = false;  U(abs(G.ph) == 90, 2:end) = false;
        W2 = G.GRBD & U;
        sp = struct('r',R.r(k), 'cap1',R.cap(k), 'capg',G.cap, 'dcap',G.cap - R.cap(k), 'ngrb',nnz(W2), ...
                    'nopt',nnz(W2 & G.GSTAT == "OPTIMAL"), 'allglob',all(G.GSTAT(W2) == "OPTIMAL"), 'dobj',max([G.DOBJ(W2); -Inf]), ...
                    'nrepl',nnz(W2 & G.RHO ~= G.RHO0), 'npolf',nnz(W2 & ~G.PLSH), 't',toc(ts));
    end
    [R.r, ix] = sort(R.r);                                  % store the ray sorted by radius
    for f = {'cap','iso','ncert','nfeas','t'}, R.(f{1}) = R.(f{1})(ix); end
    R.rho = R.rho(ix,:);
    E = struct('id',id, 'ph',PH(id), 'th',TH(id), 'r',R.r, 'cap',R.cap, 'iso',R.iso, 'rho',R.rho, ...
               'ncert',R.ncert, 'nfeas',R.nfeas, 't',R.t, 'RL',RL, 'spot',sp, 'lev',LEV);
    if isempty(j), RAY(end+1) = E; else, RAY(j) = E; end   %#ok<SAGROW>
    save(WF, 'RAY');
    s = '';  if ~isempty(sp), s = sprintf(' | spot r=%.1f dC_f=%+.1e global %d (%d Gurobi, %.0f s)', ...
                                          sp.r, sp.dcap, sp.allglob, sp.ngrb, sp.t); end
    fprintf('[w%02d %s] ray %3d (ph %+3d th %3d): %2d points (%d new), RL = %s um, %.0f s%s | %d / %d, %.1f h\n', WID, ...
            datestr(now, 'HH:MM'), id, PH(id), TH(id), numel(R.r) - 1, numel(R.r) - n0, mat2str(RL, 4), toc(tr), s, ...
            nnz(arrayfun(@(x) isequal(x.lev, LEV), RAY)), numel(MINE), toc(tw)/3600);  %#ok<TNOW1,DATST>
end
fprintf('worker %d finished: %d rays\n', WID, numel(RAY));

% ---------------------------------------------------------------------------------------
function [R, P1] = eval_(R, P1, r, u, id, cal, FH, o)
% C_f at p = r u: stage 1 of rho_solve, extra start = solution at the nearest radius already done on this ray
    [~, k] = min(abs(R.r - r));
    o.X0G = [];  if ~isempty(P1{k}), o.X0G = P1{k}.IOPT; end
    o.SEED = id*1000 + numel(R.r);                          % reproducible, independent of the worker split
    tp = tic;  G = rho_solve(q_at_(r*u, cal, FH), o);      % o.RETRY: every direction is retried until feasible
    U = true(size(G.RHO));  U(:,end) = false;  U(abs(G.ph) == 90, 2:end) = false;
    R.r(end+1) = r;  R.cap(end+1) = G.cap;  R.iso(end+1) = G.iso;  R.rho(end+1,:) = [min(G.RHO(:)) max(G.RHO(:))];
    R.ncert(end+1) = nnz(G.CERT & U);  R.nfeas(end+1) = min(G.NFEAS(U));  R.t(end+1) = toc(tp);
    P1{end+1} = struct('IOPT',G.IOPT, 'FHAT',G.FHAT, 'LAM',G.LAM, 'RHO',G.RHO, 'ERR',G.ERR, ...
                       'MINEIG',G.MINEIG, 'NFEAS',G.NFEAS);
end

function [ra, rb, ca, cb] = bracket_(R, L)
% first bracket along the ray with C_f(ra) < L <= C_f(rb); NaN if the ray never reaches L (NaN points skipped)
    ok = isfinite(R.cap);  [rs, ix] = sort(R.r(ok));  cs = R.cap(ok);  cs = cs(ix);
    k = find(cs >= L, 1);
    if isempty(k) || k == 1, ra = NaN;  rb = NaN;  ca = NaN;  cb = NaN;  return, end
    ra = rs(k-1);  rb = rs(k);  ca = cs(k-1);  cb = cs(k);
end

function Q = q_at_(P, cal, FH)
% Q_k = F H_I' L_k(p / l_hat, e_hat) F H_I at the position P [um, actuator frame]
    L = build_L(P(:).', cal.l_m, cal.e_m, cal.Pc_base);    % 1 x 6 x 6 x 3
    Q = zeros(6,6,3);
    for k = 1:3, Lk = squeeze(L(1,:,:,k));  Q(:,:,k) = FH.' * Lk * FH;  Q(:,:,k) = (Q(:,:,k) + Q(:,:,k).')/2; end
end
