%% rbf_strat.m -- RBF parameter search at four data ranges, stratified sampling,
%%                 judged on the R <= 150 working region
%
%  SAMPLING.  The working region is never thinned; only the outside is.
%      R = 150            all 1771 nodes (already below the 5000 cap)
%      R = 250/350/450    1771 working-region nodes KEPT IN FULL
%                       + 3229 drawn uniformly at random from 150 < r <= R
%                       = 5000 per range
%  The FEM nodes are a uniform grid (N scales as R^3 to 0.3 %), so an equal-
%  probability draw keeps the outer density uniform too.  rng(0) per range, so the
%  draw repeats.  This is the point of the whole scheme: a plain 5000-node random
%  draw at R = 450 would leave only ~186 nodes inside R <= 150, i.e. 10 % of them.
%
%  SEARCH.  rho = 20:10:1000 (99) x lam = 1e-10 .. 1e1 by decades (12) = 1188 pairs
%  per range.  Cholesky per pair -- at N = 5000 that is ~0.55 s, cheaper than one
%  eigendecomposition amortised over 12 lambdas.
%
%  TWO-STAGE CRITERION.
%    1. smoothness -- d3(b.b)/ds^3 keeps one sign on all three actuator axes
%       (x_a <- P1, y_a <- P3, z_a <- P6) over |s| <= 150, THE WORKING REGION,
%       not each range's own extent
%    2. among the survivors, lowest NMAE
%
%  NMAE is always the same thing: the model's error over the 1771 nodes inside
%  R <= 150, all six excitations, vector-norm convention
%      NMAE = sum_ij ||b_pred - B_FEM|| / sum_ij ||B_FEM|| * 100
%  Note this is in-sample for every range (they all contain those nodes).
%
%  Set RUN_R before running to do a subset of the ranges; results accumulate in
%  rbf_strat.mat.   [ADDED 2026-09-17]

% RUN_R selects the ranges; RUN_NCAP the fitting-set cap (Inf = every node, which
% writes to rbf_full.mat instead of rbf_strat.mat so the two studies stay separate).
% [ADDED 2026-09-17] RUN_OUT overrides the output file name.  A long run on one
% range must not share a file with a run that is still appending ranges of its own:
% RES is read once at the start and written back whole, so two live runs on the same
% file silently drop each other's later ranges.  Give a long run its own file.
if exist('RUN_R','var'),    RUN_R_    = RUN_R;    else, RUN_R_    = [150 250 350 450]; end
if exist('RUN_NCAP','var'), RUN_NCAP_ = RUN_NCAP; else, RUN_NCAP_ = 5000;              end
if exist('RUN_OUT','var'),  RUN_OUT_  = RUN_OUT;  else, RUN_OUT_  = '';                end
if exist('RUN_RHO','var'),  RUN_RHO_  = RUN_RHO;  else, RUN_RHO_  = [];                end
if exist('RUN_LAM','var'),  RUN_LAM_  = RUN_LAM;  else, RUN_LAM_  = [];                end
if exist('RUN_NQ','var'),   RUN_NQ_   = RUN_NQ;   else, RUN_NQ_   = [];                end
clearvars -except RUN_R_ RUN_NCAP_ RUN_OUT_ RUN_RHO_ RUN_LAM_ RUN_NQ_;  clc;

MODEL='long2016_hexapole_halfcut';  GEOM='tip40um';  VARIANT='maxwell';
RL     = RUN_R_;
NCAP   = RUN_NCAP_;            % cap on the fitting set; Inf = no cap
RWORK  = 150;                  % working region: kept in full, and where everything is judged
RHOL   = 20:10:1000;
LAML   = 10.^(-10:1:1);
if ~isempty(RUN_RHO_), RHOL = RUN_RHO_(:).';  end    % [ADDED 2026-09-17] grid override
if ~isempty(RUN_LAM_), LAML = RUN_LAM_(:).';  end
POLE   = [1 3 6];   SG = [1 1 -1];   NQ = 1201;   NSKIP = 4;   SEED = 0;
if ~isempty(RUN_NQ_), NQ = RUN_NQ_;  end   % [ADDED 2026-09-17] gate sampling override
NB     = 1024;                 % block size for the in-place Phi build and Cholesky

here=fileparts(mfilename('fullpath'));  FMX=fileparts(fileparts(here));
MAIN=fileparts(fileparts(fileparts(FMX)));  FLX=fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT=fullfile(FMX,'utils','data');
if ~isempty(RUN_OUT_),   MATF = fullfile(DAT, RUN_OUT_);
elseif isinf(NCAP),      MATF = fullfile(DAT,'rbf_full.mat');
else,                    MATF = fullfile(DAT,'rbf_strat.mat');  end

cfg = model_config(MODEL, GEOM);
raw = extract_maxwell_data(cfg, 'all', VARIANT);
ad  = build_actuator_data(raw, cfg);
Pa  = ad.Pa*1e6;  Ba = ad.Ba;  rr = vecnorm(Pa,2,2);

iw  = rr <= RWORK;   Pw = Pa(iw,:);   Bw = Ba(iw,:,:);   Nw = nnz(iw);
nBw = sum(vecnorm(reshape(permute(Bw,[1 3 2]),[],3),2,2));
fprintf(['working region: %d nodes inside R <= %d um -- kept in full, and the NMAE set' newline], Nw, RWORK);

% query points for the smoothness gate: three axes, |s| <= RWORK
sA  = linspace(-RWORK, RWORK, NQ).';
AXQ = cell(1,3);
for a = 1:3, Pv = zeros(NQ,3);  Pv(:,a) = sA;  AXQ{a} = Pv;  end

RES = struct('R',{},'Nall',{},'nfit',{},'nout',{},'rho',{},'lam',{},'nmae',{},'npass',{},'ndisagree',{},'g',{});
if isfile(MATF)
    Z = load(MATF);
    % [MODIFIED 2026-09-17] NQ joins the resume key.  Without it a rerun at a new
    % gate density finds the range already listed, skips it, and leaves the caller
    % reading the OLD table while believing it is the new one -- silently, because
    % nothing about the file says which density produced it.
    if isfield(Z,'RES') && isequal(Z.RHOL,RHOL) && isequal(Z.LAML,LAML) && ...
       isequal(Z.NCAP,NCAP) && isfield(Z,'NQ') && isequal(Z.NQ,NQ)
        RES = Z.RES;
        fprintf(['resuming: already have R = %s' newline], mat2str([RES.R]));
    end
end

for R0 = RL
    if ~isempty(RES) && any([RES.R] == R0), fprintf(['R=%d already done, skipping' newline], R0);  continue, end
    t0 = tic;
    Nall = nnz(rr <= R0);
    % ---- stratified fitting set ---------------------------------------------
    if Nall <= NCAP
        Pf = Pw;  Bf = Bw;  nout = 0;                     % R = 150: nothing to thin
        if R0 > RWORK
            io = find(rr > RWORK & rr <= R0);
            Pf = [Pw; Pa(io,:)];  Bf = cat(1, Bw, Ba(io,:,:));  nout = numel(io);
        end
    else
        io = find(rr > RWORK & rr <= R0);
        rng(SEED);  pick = io(randperm(numel(io), NCAP - Nw));
        Pf = [Pw; Pa(pick,:)];   Bf = cat(1, Bw, Ba(pick,:,:));   nout = numel(pick);
    end
    nfit = size(Pf,1);
    fprintf(['R=%3d: %6d nodes in range -> fitting on %d (%d working + %d outer)' newline], ...
            R0, Nall, nfit, Nw, nout);

    % ---- distances, once per range ------------------------------------------
    % [MODIFIED 2026-09-17] D2 and eye(nfit) are gone.  At nfit = 47,693 each of
    % D2, Phi, eye and the `Phi + lam*In` temporary is 18.2 GB, so the original
    % four-copy path needs ~73 GB on a 68 GB machine and cannot run at all past
    % about N = 30,000.  Phi is now written block by block straight into one
    % preallocated array, lam goes onto the diagonal in place, and a blocked
    % Cholesky overwrites the lower triangle -- peak is ONE Phi.  Same scheme as
    % utils/rbf_field.m, which is where it was worked out.  The cost is that Phi
    % must be rebuilt for every lam (the factorisation destroys it); that is one
    % exp() over nfit^2 elements, small next to the N^3/3 of the Cholesky itself.
    d2c = sum(Pf.^2,2).';
    d2w = sum(Pw.^2,2) + d2c - 2*(Pw*Pf.');                 % working region <- fit nodes
    d2q = cell(1,3);  U = cell(1,3);
    for a = 1:3
        d2q{a} = sum(AXQ{a}.^2,2) + d2c - 2*(AXQ{a}*Pf.');
        U{a}   = AXQ{a}(:,a) - Pf(:,a).';
    end
    Y  = reshape(Bf, nfit, 18);
    nb = min(NB, nfit);
    Phi = zeros(nfit, nfit);                                % the single big array
    oLT = struct('LT',true);                                % lower-triangular solve
    oTR = struct('LT',true, 'TRANSA',true);                 % and its transpose
    fprintf(['       Phi is %.1f GB; peak is one copy of it' newline], nfit^2*8/2^30);

    % [ADDED 2026-09-17] the whole grid is recorded, not just the winner: one row
    % per (rho,lam) with the gate verdict and the NMAE.  A 1188-pair run is days
    % long and the map of where the gate passes is worth as much as the argmin.
    np_  = numel(RHOL)*numel(LAML);
    TAB  = struct('rho',nan(np_,1), 'lam',nan(np_,1), 'pass',false(np_,1), ...
                  'nmae',nan(np_,1), 'secs',nan(np_,1));
    it   = 0;
    best = struct('nm',inf,'rho',NaN,'lam',NaN,'g',[]);   npass = 0;   ndis = 0;
    for rho = RHOL
        trho = tic;
        Kw  = exp(-d2w/rho^2);
        PH  = cell(1,3);
        for a = 1:3, PH{a} = exp(-d2q{a}/rho^2);  end
        for lam = LAML
            tpair = tic;   it = it + 1;
            TAB.rho(it) = rho;   TAB.lam(it) = lam;
            % ---- Phi, one column block at a time, straight into place --------
            % Written out here rather than called as a helper ON PURPOSE: this is a
            % SCRIPT, and MATLAB's in-place optimisation only applies inside a
            % function, so passing an 18 GB Phi to a helper would copy it and undo
            % the whole point.  Indexed assignment in the script body does not copy.
            for j0 = 1:nb:nfit
                j  = j0:min(j0+nb-1, nfit);
                d2 = sum(Pf.^2,2) + sum(Pf(j,:).^2,2).' - 2*(Pf*Pf(j,:).');
                d2(d2 < 0) = 0;
                Phi(:,j) = exp(-d2/rho^2);
            end
            Phi(1:nfit+1:end) = 1 + lam;                    % phi(0) = 1, ridge in place
            % ---- blocked Cholesky, overwriting the lower triangle -------------
            for k0 = 1:nb:nfit
                j = k0:min(k0+nb-1, nfit);
                Phi(j,j) = chol(Phi(j,j), 'lower');
                r = j(end)+1:nfit;
                if ~isempty(r)
                    Phi(r,j) = Phi(r,j) / Phi(j,j).';
                    for c0 = k0+nb : nb : nfit
                        c  = c0:min(c0+nb-1, nfit);
                        rr = c(1):nfit;
                        Phi(rr,c) = Phi(rr,c) - Phi(rr,j) * Phi(c,j).';
                    end
                end
            end
            % linsolve reads only the lower triangle, so the stale upper half left
            % over from the build is never touched (and is overwritten next pair).
            W = linsolve(Phi, linsolve(Phi, Y, oLT), oTR);
            % stage 1: smoothness on the three axes, |s| <= RWORK
            ok = true;  okc = true;
            for a = 1:3
                k = POLE(a);  Wk = W(:,(k-1)*3+(1:3));  ph = PH{a};  Ua = U{a};
                b0 = ph*Wk;
                b1 = ((-2/rho^2)*(Ua.*ph))*Wk;
                b2 = ((-2/rho^2)*ph + (4/rho^4)*(Ua.^2.*ph))*Wk;
                b3 = ((12/rho^4)*(Ua.*ph) - (8/rho^6)*(Ua.^3.*ph))*Wk;
                v  = SG(a)*2*sum(3*b1.*b2 + b0.*b3, 2);
                if nsc_(v)             > 0, ok  = false;  end        % dense verdict
                if nsc_(v(1:NSKIP:end))> 0, okc = false;  end        % coarse, same data
                if ~ok && ~okc, break, end
            end
            if ok ~= okc, ndis = ndis + 1;  end      % densities disagreed here
            TAB.pass(it) = ok;   TAB.secs(it) = toc(tpair);
            if ~ok, continue, end
            npass = npass + 1;
            % stage 2: NMAE on the 1771 working-region nodes
            Bp = reshape(Kw*W, Nw, 3, 6);
            nm = sum(vecnorm(reshape(permute(Bp-Bw,[1 3 2]),[],3),2,2)) / nBw * 100;
            TAB.nmae(it) = nm;
            if nm < best.nm
                k = POLE(1);  Wk = W(:,(k-1)*3+(1:3));  ph = PH{1};  Ua = U{1};
                b0 = ph*Wk;   b1 = ((-2/rho^2)*(Ua.*ph))*Wk;
                best = struct('nm',nm,'rho',rho,'lam',lam,'g',2*sum(b0.*b1,2));
            end
        end
        clear Kw PH
        % [ADDED 2026-09-17] save after every rho.  At R = 450 one rho is ~53 min
        % and the whole grid is ~4 days; writing only at the end would put the
        % entire run at the mercy of one crash.
        RESp = [RES, struct('R',R0,'Nall',Nall,'nfit',nfit,'nout',nout, ...
                            'rho',best.rho,'lam',best.lam,'nmae',best.nm, ...
                            'npass',npass,'ndisagree',ndis,'g',best.g)];
        save(MATF, 'RESp','RES','TAB','RHOL','LAML','NCAP','RWORK','SEED','sA','Nw','NQ','it');
        fprintf(['       rho %g done (%d of %d pairs, %d pass so far) [%.0f s/rho, %.1f h elapsed]' newline], ...
                rho, it, np_, npass, toc(trho), toc(t0)/3600);
    end
    clear Phi

    RES(end+1) = struct('R',R0,'Nall',Nall,'nfit',nfit,'nout',nout, ...
                        'rho',best.rho,'lam',best.lam,'nmae',best.nm,'npass',npass, ...
                        'ndisagree',ndis,'g',best.g); %#ok<SAGROW>
    fprintf(['       gate density check: %d of %d pairs judged differently at %d vs %d points' newline], ...
            ndis, numel(RHOL)*numel(LAML), numel(1:NSKIP:NQ), NQ);
    if isfinite(best.nm)
        fprintf(['       -> rho %g, lam %.0e : NMAE(R<=%d) = %.4f %%  (%d of %d pass)  [%.0f s]' newline], ...
                best.rho, best.lam, RWORK, best.nm, npass, numel(RHOL)*numel(LAML), toc(t0));
    else
        fprintf(['       -> NO pair passes the gate  [%.0f s]' newline], toc(t0));
    end
    save(MATF, 'RES','TAB','RHOL','LAML','NCAP','RWORK','SEED','sA','Nw','NQ');
    clear d2w d2q U Y Pf Bf
end
fprintf(['done' newline]);

% ---- local helper -------------------------------------------------------------
function n = nsc_(v)
%NSC_  number of sign changes in v, ignoring exact zeros
    g = sign(v);  nz = find(g~=0);  n = sum(diff(g(nz)) ~= 0);
end
