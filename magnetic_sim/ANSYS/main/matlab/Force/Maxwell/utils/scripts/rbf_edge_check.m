%% rbf_edge_check.m -- is the RBF winner of range_scan.m actually the best?
%
%  range_scan.m searches a COARSE grid (rho steps of 1.33 .. 2x, lambda steps of
%  3.16x) and reports the best pair that passes the smoothness gate.  Two things
%  that search cannot rule out:
%
%      (a) the winner sits ON or NEXT TO a grid edge, so the true optimum may lie
%          outside the box      -- R = 50, 60, 200 pick lambda two or three steps
%                                  off the lower edge; R = 60 picks rho one step
%                                  off the upper edge
%      (b) the optimum falls BETWEEN two grid points
%
%  This script re-searches a local box around each radius's winner:
%      rho   x [1/2, 1/1.41, 1, 1.41, 2]           (geometric, half-step refinement)
%      lam   x 10.^(-1.5 : 0.25 : 1.5)             (quarter-decade, 1.5 decades out)
%  Both may leave the original grid, which is the point.  Same two-stage rule and
%  the same 80/20 split as the scan, so the numbers are directly comparable.
%
%  Reports, per radius: the coarse winner, the refined winner, and the improvement.
%  An improvement of a few 1e-4 % is noise-level; anything larger means the coarse
%  grid was actually leaving accuracy on the table.
%
%  [ADDED 2026-09-16]

clear;  clc;

MODEL   = 'long2016_hexapole_halfcut';
GEOM    = 'tip40um';
VARIANT = 'maxwell';
FRHO    = [1/2 1/sqrt(2) 1 sqrt(2) 2];
FLAM    = 10.^(-1.5 : 0.25 : 1.5);
POLE    = [1 3 6];   SG = [1 1 -1];   NQ = 301;
FTRAIN  = 0.8;       SEED = 0;

here = fileparts(mfilename('fullpath'));
FMX  = fileparts(fileparts(here));
MAIN = fileparts(fileparts(fileparts(FMX)));
FLX  = fullfile(MAIN, 'matlab', 'Flux', 'Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT  = fullfile(FMX,'utils','data');

S = load(fullfile(DAT,'range_scan.mat'));
cfg = model_config(MODEL, GEOM);
raw = extract_maxwell_data(cfg, 'all', VARIANT);
ad  = build_actuator_data(raw, cfg);
Pa  = ad.Pa*1e6;  Ba = ad.Ba;  rr = vecnorm(Pa,2,2);

fprintf('%5s | %8s %6s %9s | %8s %6s %9s | %10s\n', ...
        'R','coarse','rho','lam','refined','rho','lam','gain');
REF = struct('R',{},'coarse',{},'fine',{},'rho',{},'lam',{});
for q = 1:numel(S.RES)
    r  = S.RES(q);   R0 = r.R;
    if ~isfinite(r.rbf.nm), fprintf('%5d | (no coarse winner)\n', R0);  continue, end

    in = rr <= R0;  P0 = Pa(in,:);  B0 = Ba(in,:,:);  N0 = nnz(in);
    rng(SEED);  ix = randperm(N0);  n1 = round(FTRAIN*N0);  i1 = ix(1:n1);  i2 = ix(n1+1:end);
    if n1 > S.NTRMAX, i1 = i1(1:S.NTRMAX);  end
    Ptr = P0(i1,:);  Btr = B0(i1,:,:);  ntr = numel(i1);
    Pte0 = P0(i2,:); Bte0 = B0(i2,:,:); nte = numel(i2);
    nB0 = sum(vecnorm(reshape(permute(Bte0,[1 3 2]),[],3),2,2));
    sA0 = linspace(-R0, R0, NQ).';
    AXQ = cell(1,3);  for a = 1:3, Pq = zeros(NQ,3);  Pq(:,a) = sA0;  AXQ{a} = Pq;  end

    d2c = sum(Ptr.^2,2).';
    D2t = sum(Ptr.^2,2) + d2c - 2*(Ptr*Ptr.');  D2t(D2t<0) = 0;
    d2e = sum(Pte0.^2,2) + d2c - 2*(Pte0*Ptr.');
    Ytr0 = reshape(Btr, ntr, 18);

    best = struct('nm',inf,'rho',NaN,'lam',NaN);
    for rho = r.rbf.rho * FRHO
        Phi = exp(-D2t/rho^2);  [V,d] = eig((Phi+Phi.')/2,'vector');  VtY = V.'*Ytr0;
        Kte = exp(-d2e/rho^2);
        PH = cell(1,3);  UU = cell(1,3);
        for a = 1:3
            d2q = sum(AXQ{a}.^2,2) + d2c - 2*(AXQ{a}*Ptr.');
            PH{a} = exp(-d2q/rho^2);  UU{a} = AXQ{a}(:,a) - Ptr(:,a).';
        end
        for lam = r.rbf.lam * FLAM
            if min(d) + lam <= 0, continue, end        % guard: shifted system must stay PD
            W = V*(VtY ./ (d + lam));   ok = true;
            for a = 1:3
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
            Bp = reshape(Kte*W, nte, 3, 6);
            nm = sum(vecnorm(reshape(permute(Bp-Bte0,[1 3 2]),[],3),2,2)) / nB0 * 100;
            if nm < best.nm, best = struct('nm',nm,'rho',rho,'lam',lam);  end
        end
        clear Phi V d
    end
    REF(end+1) = struct('R',R0,'coarse',r.rbf.nm,'fine',best.nm,'rho',best.rho,'lam',best.lam); %#ok<SAGROW>
    fprintf('%5d | %7.4f%% %6g %9.1e | %7.4f%% %6.0f %9.1e | %+9.4f\n', ...
            R0, r.rbf.nm, r.rbf.rho, r.rbf.lam, best.nm, best.rho, best.lam, best.nm - r.rbf.nm);
    save(fullfile(DAT,'rbf_edge_check.mat'), 'REF','FRHO','FLAM');
end
g = [REF.fine] - [REF.coarse];
fprintf(['\nlargest improvement %.4f %% at R = %d   |   median %.5f %%' newline], ...
        min(g), REF(find(g==min(g),1)).R, median(g));
