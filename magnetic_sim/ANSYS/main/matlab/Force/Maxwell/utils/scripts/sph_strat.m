%% sph_strat.m -- solid harmonics at four data ranges, same stratified sampling and
%%                 same criterion as rbf_strat.m, so the two are directly comparable
%
%  SAMPLING -- identical to rbf_strat.m, same seed, same draw:
%      R = 150            all 1771 nodes
%      R = 250/350/450    1771 working-region nodes KEPT IN FULL
%                       + 3229 drawn uniformly from 150 < r <= R   (rng(0))
%                       = 5000 per range
%  The draw is reproduced line for line rather than shared through a helper, so the
%  two searches cannot drift apart.
%
%  SEARCH.  Degree L is swept from 1 to LMAX = 30, or until K = (L+1)^2-1 exceeds
%  the fitting-set size, whichever comes first.  "At least 1000 steps further"
%  cannot apply to L: even the mathematical ceiling is only 69 degrees at N = 5000.
%  The cap of 30 is an evidence-based cut, not an arbitrary one --
%      * the earlier 36-range sweep (cap 24) never selected an L above 21
%      * R = 150 here swept all 41 admissible degrees and chose L = 9, with only
%        9 degrees passing the gate at all
%      * L = 60..69 alone consumed most of the runtime of a single range, because
%        the design matrix is 15000 x 4899 there
%  Raising the cap back to the ceiling costs 2-3 hours and, on this evidence,
%  cannot change the answer.  [user's call, 2026-09-17]
%
%  TWO-STAGE CRITERION, same as the RBF run:
%    1. d3(b.b)/ds^3 keeps one sign on all three actuator axes over |s| <= 150
%       (x_a <- P1, y_a <- P3, z_a <- P6)
%    2. among the survivors, lowest NMAE on the 1771 nodes inside R <= 150,
%       six excitations, vector-norm convention
%
%  Also stores d(b.b)/dx_a along x_a over |x_a| <= 150 for the figure, and the
%  eighteen-parameter curve from the flux calibration for comparison.
%
%  GATE SAMPLING DENSITY.  The gate is evaluated on NQ = 1201 points per axis
%  (0.25 um spacing).  Because 1201 = 4*300+1, taking every 4th point gives exactly
%  the 301-point grid used before, so both densities are scored from the SAME
%  derivative evaluation at no extra cost, and the run reports whether they ever
%  disagree.  That turns the sampling density from a hidden assumption into a
%  checked one: if the coarse and fine verdicts always match, the density is not a
%  factor in any of the results.
%
%  Writes data/long2016_hexapole_halfcut/.mat/sph_strat.mat   [ADDED 2026-09-17]

% RUN_R / RUN_NCAP / RUN_LMAX mirror rbf_strat.m.  NCAP = Inf means every node and
% sends the results to sph_full.mat, keeping the capped study in sph_strat.mat.
if exist('RUN_R','var'),     RUN_R_     = RUN_R;     else, RUN_R_     = [150 250 350 450]; end
if exist('RUN_NCAP','var'),  RUN_NCAP_  = RUN_NCAP;  else, RUN_NCAP_  = 5000;              end
if exist('RUN_LMAX','var'),  RUN_LMAX_  = RUN_LMAX;  else, RUN_LMAX_  = 30;                end
% [ADDED 2026-09-17] RUN_LMAX = Inf means NO degree cap: climb until patience runs
% out.  RUN_OUT gives a long run its own file so it cannot overwrite a run that is
% still appending ranges to sph_full.mat.
if exist('RUN_OUT','var'),   RUN_OUT_   = RUN_OUT;   else, RUN_OUT_   = '';                end
if exist('RUN_NQ','var'),    RUN_NQ_    = RUN_NQ;    else, RUN_NQ_    = [];                end
clearvars -except RUN_R_ RUN_NCAP_ RUN_LMAX_ RUN_OUT_ RUN_NQ_;  clc;

MODEL='long2016_hexapole_halfcut';  GEOM='tip40um';  VARIANT='maxwell';
RL      = RUN_R_;
NCAP    = RUN_NCAP_;           % fitting-set cap; Inf = every node
RWORK   = 150;                 % working region: kept in full, and where all judging happens
SEED    = 0;                   % same seed as the RBF run
POLE    = [1 3 6];   SG = [1 1 -1];   NQ = 1201;   NSKIP = 4;
if ~isempty(RUN_NQ_), NQ = RUN_NQ_;  end   % [ADDED 2026-09-17] gate sampling override
NS      = 601;                 % samples of the plotted gradient curve
LMAX    = RUN_LMAX_;           % degree cap (patience usually stops first)
PATIENCE = 10;                 % consecutive failed degrees before the sweep gives up
CALNAME = 'calib_current_maxwell_R150_eighteen.mat';

here=fileparts(mfilename('fullpath'));  FMX=fileparts(fileparts(here));
MAIN=fileparts(fileparts(fileparts(FMX)));  FLX=fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT=fullfile(FMX,'utils','data');
if ~isempty(RUN_OUT_),   MATF = fullfile(DAT, RUN_OUT_);
elseif isinf(NCAP),      MATF = fullfile(DAT,'sph_full.mat');
else,                    MATF = fullfile(DAT,'sph_strat.mat');  end

cfg = model_config(MODEL, GEOM);
raw = extract_maxwell_data(cfg, 'all', VARIANT);
ad  = build_actuator_data(raw, cfg);
Pa  = ad.Pa*1e6;  Ba = ad.Ba;  rr = vecnorm(Pa,2,2);

iw  = rr <= RWORK;   Pw = Pa(iw,:);   Bw = Ba(iw,:,:);   Nw = nnz(iw);
nBw = sum(vecnorm(reshape(permute(Bw,[1 3 2]),[],3),2,2));
fprintf(['working region: %d nodes inside R <= %d um' newline], Nw, RWORK);

sA  = linspace(-RWORK, RWORK, NS).';         % the plotted x_a segment (NS points)
Pq  = [sA zeros(NS,2)];
Ix  = zeros(6,1);  Ix(1) = 1;                % P1 at 1 A
sG  = linspace(-RWORK, RWORK, NQ).';         % gate samples (same span)
AXQ = cell(1,3);
for a = 1:3, Pv = zeros(NQ,3);  Pv(:,a) = sG;  AXQ{a} = Pv;  end

RES = struct('R',{},'Nall',{},'nfit',{},'nout',{},'L',{},'K',{},'nmae',{},'npass',{},'Ltried',{},'ndisagree',{},'g',{});
for R0 = RL
    t0 = tic;
    Nall = nnz(rr <= R0);
    % ---- stratified fitting set: line for line the same as rbf_strat.m -------
    if Nall <= NCAP
        Pf = Pw;  Bf = Bw;  nout = 0;
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
    fprintf(['R=%3d: %6d in range -> fitting on %d (%d working + %d outer)' newline], ...
            R0, Nall, nfit, Nw, nout);

    best = struct('nm',inf,'L',NaN,'K',NaN,'g',[]);   npass = 0;   Ltried = 0;   ndis = 0;   miss = 0;
    % [MODIFIED 2026-09-17 user's call] NO degree cap when LMAX = Inf: L climbs
    % until PATIENCE consecutive failures, full stop.  Two guards remain, and
    % neither is a criterion -- they are the points past which the fit stops being
    % defined or stops fitting in memory:
    %   K <= nfit    below this the least-squares system is underdetermined
    %   JCAP         the design matrix is 3*nfit x K in doubles, and sph_field
    %                holds D, J and Js at once, so the real peak is ~3x that.  At
    %                nfit = 47,693 the cap below is hit at L ~ 82, far beyond any
    %                degree ever selected (the highest so far is 29).
    % [ADDED 2026-09-17] every degree is logged, not just the winner: L, K, the gate
    % verdict, the NMAE and the seconds.  With no cap the sweep's own trace is the
    % evidence for where it stopped and why.
    LTAB = zeros(0,5);
    Lq = 0;   JCAP = 8e9;                          % bytes, one copy of J
    while true
        Lq = Lq + 1;
        if isfinite(LMAX) && Lq > LMAX, break, end
        Kq = (Lq+1)^2 - 1;
        if Kq > nfit
            fprintf(['       stop: K = %d exceeds nfit = %d (underdetermined)' newline], Kq, nfit);
            break
        end
        if 3*nfit*Kq*8 > JCAP
            fprintf(['       stop: L = %d would need a %.1f GB design matrix' newline], ...
                    Lq, 3*nfit*Kq*8/2^30);
            break
        end
        Ltried = Lq;
        [bq, gq, iq] = sph_field(R0, struct('P',Pf,'B',Bf), ...
                                 struct('L',Lq,'Rn',R0,'quiet',true,'condmax',1500));
        % stage 1: smoothness on the three axes over the WORKING region
        ok = true;  okc = true;
        for a = 1:3
            I  = zeros(6,1);  I(POLE(a)) = 1;
            d3 = gq(AXQ{a}, I, 'd3');
            v  = SG(a)*d3(:,a);
            if nsc_(v)              > 0, ok  = false;  end       % dense verdict
            if nsc_(v(1:NSKIP:end)) > 0, okc = false;  end       % coarse, same data
            if ~ok && ~okc, break, end
        end
        if ok ~= okc, ndis = ndis + 1;  end
        LTAB(end+1,:) = [Lq, iq.K, double(ok), NaN, iq.secs]; %#ok<SAGROW>
        vstr = {'FAIL','pass'};
        fprintf(['       L = %3d (K = %5d): gate %s  [%.0f s, %.2f h elapsed]' newline], ...
                Lq, iq.K, vstr{ok+1}, iq.secs, toc(t0)/3600);
        % HIGHEST STILL-SMOOTH DEGREE, WITH PATIENCE (user's call, 2026-09-17).
        % L climbs from 1.  Every degree that passes the gate becomes the current
        % answer and resets the counter; a degree that fails only increments it.  The
        % sweep gives up after PATIENCE consecutive failures, so a single awkward
        % degree cannot end the search while higher ones are still smooth.
        % Not the lowest smooth degree (that is always L = 1 -- a linear potential
        % gives a constant b, so the gradient is zero and can never change sign), and
        % not the lowest-NMAE smooth degree (the earlier rule, which ran to L = 29).
        % NMAE is computed and reported but no longer selects anything.
        if ~ok
            miss = miss + 1;
            if miss >= PATIENCE && isfinite(best.nm), break, end
            continue
        end
        miss  = 0;
        npass = npass + 1;
        Bp = zeros(Nw,3,6);
        for k = 1:6, I = zeros(6,1); I(k) = 1;  Bp(:,:,k) = bq(Pw, I);  end
        nm = sum(vecnorm(reshape(permute(Bp-Bw,[1 3 2]),[],3),2,2)) / nBw * 100;
        G  = gq(Pq, Ix, 'du');
        LTAB(end,4) = nm;
        best = struct('nm',nm,'L',Lq,'K',iq.K,'g',G(:,1));
        % save on every degree that passes -- an uncapped sweep has no known end
        RESp = [RES, struct('R',R0,'Nall',Nall,'nfit',nfit,'nout',nout,'L',best.L, ...
                            'K',best.K,'nmae',best.nm,'npass',npass,'Ltried',Ltried, ...
                            'ndisagree',ndis,'g',best.g)];
        save(MATF, 'RESp','RES','LTAB','RL','NCAP','RWORK','SEED','sA','Nw','LMAX','PATIENCE','NQ');
    end

    RES(end+1) = struct('R',R0,'Nall',Nall,'nfit',nfit,'nout',nout,'L',best.L,'K',best.K, ...
                        'nmae',best.nm,'npass',npass,'Ltried',Ltried,'ndisagree',ndis,'g',best.g); %#ok<SAGROW>
    fprintf(['       gate density check: %d of %d degrees judged differently at %d vs %d points' newline], ndis, Ltried, numel(1:NSKIP:NQ), NQ);
    fprintf(['       -> L = %d (K = %d), NMAE(R<=%d) = %.4f %%  (%d of %d degrees pass)  [%.0f s]' newline], ...
            best.L, best.K, RWORK, best.nm, npass, Ltried, toc(t0));
    save(MATF, 'RES','LTAB','RL','NCAP','RWORK','SEED','sA','Nw','LMAX','PATIENCE','NQ');
end

% ---- the eighteen-parameter curve, for the figure and the relative errors -----
C   = load(fullfile(FLX,'data',MODEL,'.mat',CALNAME));
lum = C.l_hat*1e6;
assert(C.USE_BIAS == 1, 'sph_strat:notEighteen', 'calibration is not the 18-parameter one');
Lt  = build_L(Pq, lum, C.e, cfg.Pc_base);
Gc  = base_model(Lt, C.KI_bar, C.Fmap(:,1), C.gI_hat^2 / lum);
g18 = squeeze(Gc(1,:,1)).';
fprintf([newline 'eighteen parameters: %.4f .. %.4f mT^2/um (l_hat %.2f um, gI_hat %.4f)' newline], ...
        min(g18), max(g18), lum, C.gI_hat);

fprintf([newline '%5s %5s %6s %14s %14s %14s' newline], 'R','L','K','NMAE(R<=150)%','grad(+150)','rel.err vs 18');
for q = 1:numel(RES)
    gg = RES(q).g;
    fprintf('%5d %5d %6d %14.4f %14.4f %13.2f %%\n', RES(q).R, RES(q).L, RES(q).K, RES(q).nmae, ...
            gg(end), (gg(end)-g18(end))/g18(end)*100);
end

save(MATF, 'RES','LTAB','RL','NCAP','RWORK','SEED','sA','Nw','g18','lum','CALNAME','LMAX','PATIENCE','NQ');
fprintf([newline 'wrote %s' newline], MATF);

% ---- local helper -------------------------------------------------------------
function n = nsc_(v)
%NSC_  number of sign changes in v, ignoring exact zeros
    g = sign(v);  nz = find(g~=0);  n = sum(diff(g(nz)) ~= 0);
end
