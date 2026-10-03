%% range_scan.m -- test NMAE vs data range, RBF against solid harmonics
%
%  For every radius R0 the FEM nodes with |p| <= R0 are split 80 / 20 (rng 0, the
%  same split for both models).  Each model is then tuned on the TRAIN part under
%  a two-stage rule:
%
%      stage 1   the gradient must be smooth:  d3(b.b)/ds^3 keeps a constant sign
%                along all three actuator axes over |s| <= R0  (301 samples)
%      stage 2   among the parameters that pass, take the one with the lowest
%                TEST NMAE
%
%  Parameters searched:  RBF  (rho, lambda) grid    sph  degree L = 1 .. 16
%  The training set is capped at NTRMAX = 4000 nodes so the RBF stays solvable; the
%  cap first binds at the radius where 80 % of the nodes first exceeds 4000, which the
%  run reports explicitly.  The test set is ALWAYS the full 20 %, never capped.
%
%  NMAE = sum_ij ||b_pred - B_FEM|| / sum_ij ||B_FEM|| * 100, the project's
%  vector-norm convention (calibration-transfer-matrix-output, appendix 3).
%
%  Writes data/long2016_hexapole_halfcut/.mat/range_scan.mat (after every radius,
%  so a long run can be inspected while it is still going).
%
%  Recovered from the 2026-09-15 inline run and kept verbatim except for the
%  data-loading preamble and this header.  [ADDED 2026-09-16]

clear;  clc;

%% ---- per-run knobs ---------------------------------------------------------
MODEL   = 'long2016_hexapole_halfcut';
GEOM    = 'tip40um';
VARIANT = 'maxwell';                  % No gap baseline
RL      = [25 30 40 50 60 70 80 100 120 150 180 200 250 300 350 400 450 500];
% [MODIFIED 2026-09-16] grids set by the user after the rho-dominance experiment.
%   rho 20:10:1000 -- linear, 99 points.  The upper bound was checked: capping rho at
%   1000 costs 0.00 / 0.78 / 0.19 / 0.00 % relative at R = 25 / 60 / 80 / 150, because
%   NMAE is very flat in rho up there (rho 850 and 1300 give the same value).  The
%   smooth window starts at rho = 140 .. 330 depending on radius, so 20 .. 130 is a
%   known-empty stretch kept as evidence.
%   lam 1e-10 .. 1e1 in decades -- 12 points.  NOTE (raised and overruled): the
%   smooth window can be as narrow as 1.26x (R = 450 at rho = 566), so a 10x step can
%   step over it; the dense rho grid may or may not compensate, since the window moves
%   with rho.  Watch for radii reporting "none smooth".
RHOL    = 20:10:1000;
LAML    = 10.^(-10:1:1);
% [MODIFIED 2026-09-16] L is swept exhaustively instead of to a fixed small cap: the
%   loop runs to LMAX and is cut short by K = (L+1)^2-1 <= ntr, which binds first at
%   every radius (L <= 1 at R = 25, <= 36 at R = 150, <= 62 at R = 500).  "1000 steps"
%   cannot apply to L -- there are fewer than 100 feasible degrees -- so every feasible
%   degree is tried and none is skipped, which is stricter than stopping after 1000.
LMAX    = 40;                         % highest harmonic degree tried
NTRMAX  = 4000;                       % training-set cap
FTRAIN  = 0.8;
SEED    = 0;
POLE    = [1 3 6];                    % x_a <- P1, y_a <- P3, z_a <- P6
SG      = [1 1 -1];                   % P6 sits on -z_a, so that axis is flipped
NQ      = 301;                        % samples of the smoothness test per axis

%% ---- paths -----------------------------------------------------------------
here = fileparts(mfilename('fullpath'));            % .../Force/Maxwell/temp_code
FMX  = fileparts(fileparts(here));                             % .../Force/Maxwell
MAIN = fileparts(fileparts(fileparts(FMX)));        % .../ANSYS/main
FLX  = fullfile(MAIN, 'matlab', 'Flux', 'Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT  = fullfile(FMX,'utils','data');

%% ---- FEM field: actuator frame, um, mT, all-source --------------------------
cfg = model_config(MODEL, GEOM);
raw = extract_maxwell_data(cfg, 'all', VARIANT);
ad  = build_actuator_data(raw, cfg);
Pa  = ad.Pa*1e6;  Ba = ad.Ba;  rr = vecnorm(Pa,2,2);
fprintf(['loaded %s / %s : %d nodes, |B| max %.4f mT' newline], MODEL, VARIANT, size(Pa,1), max(vecnorm(reshape(permute(Ba,[1 3 2]),[],3),2,2)));

%% ---- the scan --------------------------------------------------------------
RES = struct('R',{},'N',{},'ntr',{},'nte',{},'capped',{},'rbf',{},'sph',{},'nrpass',{},'nspass',{});
edgehit = {};                         % radii whose winner sits on a grid boundary

% [ADDED 2026-09-16] resume.  A full sweep is hours, and the .mat is written after
%   every radius, so an interrupted run picks up where it stopped instead of redoing
%   what is already there.  Only resumes when every knob matches, so a changed grid
%   silently starts over rather than mixing two searches in one file.
RESFILE = fullfile(DAT,'range_scan.mat');
if isfile(RESFILE)
    Z = load(RESFILE);
    same = isfield(Z,'RES') && isfield(Z,'LMAX') && isequal(Z.RHOL,RHOL) && ...
           isequal(Z.LAML,LAML) && isequal(Z.LMAX,LMAX) && isequal(Z.NTRMAX,NTRMAX) && ...
           isequal(Z.FTRAIN,FTRAIN) && isequal(Z.SEED,SEED);
    if same && ~isempty(Z.RES)
        RES = Z.RES;
        fprintf(['resuming: %d radii already done %s' newline], numel(RES), mat2str([RES.R]));
    elseif isfile(RESFILE)
        fprintf(['existing range_scan.mat uses different settings -- starting over' newline]);
    end
end
for R0 = RL
    if ~isempty(RES) && any([RES.R] == R0), continue, end     % already done, see resume
    in = rr <= R0;  P0 = Pa(in,:);  B0 = Ba(in,:,:);  N0 = nnz(in);
    rng(SEED);  ix = randperm(N0);  n1 = round(FTRAIN*N0);  i1 = ix(1:n1);  i2 = ix(n1+1:end);
    cap = n1 > NTRMAX;
    if cap, i1 = i1(1:NTRMAX);  end
    Ptr = P0(i1,:);  Btr = B0(i1,:,:);  ntr = numel(i1);
    Pte0 = P0(i2,:); Bte0 = B0(i2,:,:); nte = numel(i2);
    nB0 = sum(vecnorm(reshape(permute(Bte0,[1 3 2]),[],3),2,2));
    sA0 = linspace(-R0, R0, NQ).';
    AXQ = cell(1,3);  for a = 1:3, Pq = zeros(NQ,3);  Pq(:,a) = sA0;  AXQ{a} = Pq;  end

    % ---- RBF: Gaussian kernel, derivatives written out so the sweep stays cheap
    d2c = sum(Ptr.^2,2).';
    D2t = sum(Ptr.^2,2) + d2c - 2*(Ptr*Ptr.');  D2t(D2t<0) = 0;
    d2e = sum(Pte0.^2,2) + d2c - 2*(Pte0*Ptr.');
    Ytr0 = reshape(Btr, ntr, 18);
    bestR = struct('nm',inf,'rho',NaN,'lam',NaN);   nrp = 0;
    for rho = RHOL
        Phi = exp(-D2t/rho^2);  [V,d] = eig((Phi+Phi.')/2,'vector');  VtY = V.'*Ytr0;
        Kte = exp(-d2e/rho^2);
        PH = cell(1,3);  UU = cell(1,3);
        for a = 1:3
            d2q = sum(AXQ{a}.^2,2) + d2c - 2*(AXQ{a}*Ptr.');
            PH{a} = exp(-d2q/rho^2);  UU{a} = AXQ{a}(:,a) - Ptr(:,a).';
        end
        for lam = LAML
            W = V*(VtY ./ (d + lam));   ok = true;
            for a = 1:3                                   % stage 1: smoothness
                k = POLE(a);  Wk = W(:,(k-1)*3+(1:3));  ph = PH{a};  U = UU{a};
                b0 = ph*Wk;
                b1 = ((-2/rho^2)*(U.*ph))*Wk;
                b2 = ((-2/rho^2)*ph + (4/rho^4)*(U.^2.*ph))*Wk;
                b3 = ((12/rho^4)*(U.*ph) - (8/rho^6)*(U.^3.*ph))*Wk;
                v  = SG(a)*2*sum(3*b1.*b2 + b0.*b3, 2);
                g  = sign(v);  nz = find(g~=0);
                if sum(diff(g(nz))~=0) > 0, ok = false;  break, end
            end
            if ~ok, continue, end
            nrp = nrp + 1;
            Bp = reshape(Kte*W, nte, 3, 6);               % stage 2: test NMAE
            nm = sum(vecnorm(reshape(permute(Bp-Bte0,[1 3 2]),[],3),2,2)) / nB0 * 100;
            if nm < bestR.nm, bestR = struct('nm',nm,'rho',rho,'lam',lam);  end
        end
        clear Phi V d
    end

    % ---- solid harmonics: the only knob is the degree
    bestS = struct('nm',inf,'L',NaN,'K',NaN);   nsp = 0;
    for Lq = 1:LMAX
        if (Lq+1)^2-1 > ntr, break, end
        [btq,gtq,itq] = sph_field(R0, struct('P',Ptr,'B',Btr), struct('L',Lq,'Rn',R0,'quiet',true));
        ok = true;
        for a = 1:3                                       % stage 1: smoothness
            I = zeros(6,1);  I(POLE(a)) = 1;  d3 = gtq(AXQ{a}, I, 'd3');
            v = SG(a)*d3(:,a);  g = sign(v);  nz = find(g~=0);
            if sum(diff(g(nz))~=0) > 0, ok = false;  break, end
        end
        if ~ok, continue, end
        nsp = nsp + 1;
        Bp = zeros(nte,3,6);                              % stage 2: test NMAE
        for k = 1:6, I = zeros(6,1); I(k) = 1;  Bp(:,:,k) = btq(Pte0, I);  end
        nm = sum(vecnorm(reshape(permute(Bp-Bte0,[1 3 2]),[],3),2,2)) / nB0 * 100;
        if nm < bestS.nm, bestS = struct('nm',nm,'L',Lq,'K',itq.K);  end
    end

    RES(end+1) = struct('R',R0,'N',N0,'ntr',ntr,'nte',nte,'capped',cap, ...
                        'rbf',bestR,'sph',bestS,'nrpass',nrp,'nspass',nsp); %#ok<SAGROW>

    % a winner sitting on the edge of the grid means the grid, not the model, set it
    eh = {};
    if bestR.rho == RHOL(1),   eh{end+1} = 'rho min';  end %#ok<SAGROW>
    if bestR.rho == RHOL(end), eh{end+1} = 'rho max';  end %#ok<SAGROW>
    if bestR.lam == LAML(1),   eh{end+1} = 'lam min';  end %#ok<SAGROW>
    if bestR.lam == LAML(end), eh{end+1} = 'lam max';  end %#ok<SAGROW>
    if bestS.L   == LMAX,      eh{end+1} = 'L max';    end %#ok<SAGROW>
    if ~isempty(eh), edgehit{end+1} = sprintf('R %d: %s', R0, strjoin(eh,', '));  end %#ok<SAGROW>

    if cap, capmark = '*'; else, capmark = ' '; end
    fprintf('R %3d (N %5d, tr %5d%s, te %5d): ', R0, N0, ntr, capmark, nte);
    if isfinite(bestR.nm), fprintf('RBF rho %4g lam %.1e -> %8.4f %% (%3d pass) | ', bestR.rho, bestR.lam, bestR.nm, nrp);
    else, fprintf('RBF none smooth                        | ');  end
    if isfinite(bestS.nm), fprintf('sph L %2d (K %3d) -> %8.4f %% (%2d pass)\n', bestS.L, bestS.K, bestS.nm, nsp);
    else, fprintf('sph none smooth\n');  end
    save(fullfile(DAT,'range_scan.mat'), 'RES','RL','RHOL','LAML','NTRMAX','LMAX','FTRAIN','SEED');
end
fprintf(['done  (* = training set capped at %d)' newline], NTRMAX);

% keep the table in radius order and rebuild the edge check over EVERY radius, so a
% resumed run reports on the rows it inherited as well as the ones it just computed
[~,ord] = sort([RES.R]);   RES = RES(ord);
edgehit = {};
for q = 1:numel(RES)
    r = RES(q);  eh = {};
    if r.rbf.rho == RHOL(1),   eh{end+1} = 'rho min';  end %#ok<SAGROW>
    if r.rbf.rho == RHOL(end), eh{end+1} = 'rho max';  end %#ok<SAGROW>
    if r.rbf.lam == LAML(1),   eh{end+1} = 'lam min';  end %#ok<SAGROW>
    if r.rbf.lam == LAML(end), eh{end+1} = 'lam max';  end %#ok<SAGROW>
    if r.sph.L   == LMAX,      eh{end+1} = 'L max';    end %#ok<SAGROW>
    if ~isempty(eh), edgehit{end+1} = sprintf('R %d: %s', r.R, strjoin(eh,', '));  end %#ok<SAGROW>
end

% ---- summary table: every radius, and where the training cap starts biting --------
fprintf([newline 'SUMMARY  (80/20 split, rng %d; test set never capped)' newline], SEED);
fprintf('%6s %8s %8s %8s %6s | %10s %6s %9s | %10s %5s %5s\n', ...
        'R','N','train','test','cap','RBF NMAE','rho','lam','sph NMAE','L','K');
for q = 1:numel(RES)
    r = RES(q);
    if r.capped, cm = 'YES'; else, cm = '-'; end
    fprintf('%6d %8d %8d %8d %6s | ', r.R, r.N, r.ntr, r.nte, cm);
    if isfinite(r.rbf.nm), fprintf('%9.4f%% %6g %9.1e | ', r.rbf.nm, r.rbf.rho, r.rbf.lam);
    else,                  fprintf('%9s  %6s %9s | ', 'none', '-', '-');  end
    if isfinite(r.sph.nm), fprintf('%9.4f%% %5d %5d\n', r.sph.nm, r.sph.L, r.sph.K);
    else,                  fprintf('%9s  %5s %5s\n', 'none', '-', '-');  end
end
ic = find([RES.capped], 1);
if isempty(ic)
    fprintf(['training cap %d never binds: the largest training set is %d nodes' newline], ...
            NTRMAX, max([RES.ntr]));
else
    fprintf(['training cap %d first binds at R = %d um (80%% of %d nodes = %d > %d);' ...
             ' R <= %d um use the full 80%%' newline], NTRMAX, RES(ic).R, RES(ic).N, ...
            round(FTRAIN*RES(ic).N), NTRMAX, RES(ic-1).R);
end
if isempty(edgehit)
    fprintf(['grid check: no winner sits on a grid boundary -- the search space is wide enough' newline]);
else
    fprintf(['grid check: WINNER ON GRID EDGE, widen the search --' newline]);
    fprintf(['  %s' newline], edgehit{:});
end
