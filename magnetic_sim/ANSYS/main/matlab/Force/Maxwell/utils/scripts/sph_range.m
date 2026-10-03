%% sph_range.m -- solid harmonics fitted at every range R = 150:10:500, scored on
%%                 the fixed R <= 150 working region
%
%  FOR EACH RANGE R:
%    1. fit on ALL nodes with |p| <= R, six excitations at once, for L = 1 .. LMAX
%    2. SMOOTHNESS GATE -- keep only the degrees whose gradient has d3(b.b)/ds^3 of
%       constant sign on all three actuator axes over |s| <= R (the range's own
%       extent, same definition the earlier scans used)
%    3. among the survivors take the lowest NMAE, evaluated ALWAYS on the same
%       fixed set: every node with |p| <= 150 um (1771 of them, six excitations)
%
%  NMAE = sum_ij ||b_pred - B_FEM|| / sum_ij ||B_FEM|| * 100, the project's
%  vector-norm convention.
%
%  NOTE THIS IS IN-SAMPLE.  Every range R >= 150 contains the whole evaluation set
%  inside its own fitting data, so the number says "how well does a model trained on
%  a wider region still describe the working region", not "how well does it
%  generalise".  A rise with R therefore means the extra data is pulling the fit
%  away from the working region, not that it is overfitting.
%
%  Also stores d(b.b)/dx_a along the x_a axis (P1 excited, 1 A) over |s| <= 150 for
%  every range, for the curve plot.
%
%  Writes data/long2016_hexapole_halfcut/.mat/sph_range.mat   [ADDED 2026-09-16]

clear;  clc;

MODEL='long2016_hexapole_halfcut';  GEOM='tip40um';  VARIANT='maxwell';
RL    = 150:10:500;            % the ranges
LMAX  = 24;                    % highest degree tried (K = 624; winners so far <= 16)
REVAL = 150;                   % the fixed evaluation region [um]
SPAN  = 150;   NS = 601;       % gradient curve segment
POLE  = [1 3 6];  SG = [1 1 -1];  NQ = 301;   % axes for the smoothness gate

here=fileparts(mfilename('fullpath'));  FMX=fileparts(fileparts(here));
MAIN=fileparts(fileparts(fileparts(FMX)));  FLX=fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT=fullfile(FMX,'utils','data');

cfg = model_config(MODEL, GEOM);
raw = extract_maxwell_data(cfg, 'all', VARIANT);
ad  = build_actuator_data(raw, cfg);
Pa  = ad.Pa*1e6;  Ba = ad.Ba;  rr = vecnorm(Pa,2,2);

% the fixed evaluation set
ie   = rr <= REVAL;   Pe = Pa(ie,:);   Be = Ba(ie,:,:);   Ne = nnz(ie);
nBe  = sum(vecnorm(reshape(permute(Be,[1 3 2]),[],3),2,2));
fprintf(['evaluation set: %d nodes inside R <= %d um' newline], Ne, REVAL);

s  = linspace(-SPAN, SPAN, NS).';
Pq = [s zeros(NS,2)];                       % x_a axis, for the curves
Ix = zeros(6,1);  Ix(1) = 1;                % P1 at 1 A

RES = struct('R',{},'N',{},'L',{},'K',{},'nmae',{},'npass',{},'g',{},'Ltried',{},'nmae_all',{});
for R0 = RL
    t0 = tic;
    in = rr <= R0;   P0 = Pa(in,:);   B0 = Ba(in,:,:);   N0 = nnz(in);
    AXQ = cell(1,3);  sA = linspace(-R0, R0, NQ).';
    for a = 1:3, Pv = zeros(NQ,3);  Pv(:,a) = sA;  AXQ{a} = Pv;  end

    best = struct('nm',inf,'L',NaN,'K',NaN,'g',[]);
    npass = 0;   Ltried = 0;   nm_all = nan(1,LMAX);
    for Lq = 1:LMAX
        if (Lq+1)^2-1 > N0, break, end
        Ltried = Lq;
        [bq, gq, iq] = sph_field(R0, struct('P',P0,'B',B0), ...
                                 struct('L',Lq,'Rn',R0,'quiet',true));
        % NMAE on the fixed evaluation set (computed for every L, for the record)
        Bp = zeros(Ne,3,6);
        for k = 1:6, I = zeros(6,1); I(k) = 1;  Bp(:,:,k) = bq(Pe, I);  end
        nm = sum(vecnorm(reshape(permute(Bp-Be,[1 3 2]),[],3),2,2)) / nBe * 100;
        nm_all(Lq) = nm;
        % smoothness gate over the range's own extent
        ok = true;
        for a = 1:3
            I = zeros(6,1);  I(POLE(a)) = 1;
            d3 = gq(AXQ{a}, I, 'd3');
            v  = SG(a)*d3(:,a);   g = sign(v);   nz = find(g~=0);
            if sum(diff(g(nz))~=0) > 0, ok = false;  break, end
        end
        if ~ok, continue, end
        npass = npass + 1;
        if nm < best.nm
            G = gq(Pq, Ix, 'du');
            best = struct('nm',nm,'L',Lq,'K',iq.K,'g',G(:,1));
        end
    end

    RES(end+1) = struct('R',R0,'N',N0,'L',best.L,'K',best.K,'nmae',best.nm, ...
                        'npass',npass,'g',best.g,'Ltried',Ltried,'nmae_all',nm_all); %#ok<SAGROW>
    if isfinite(best.nm)
        fprintf(['R=%3d (%6d nodes, L tried 1..%2d): L=%2d K=%3d  NMAE(R<=150) = %8.4f %%  (%d of %d pass)  [%.0f s]' newline], ...
                R0, N0, Ltried, best.L, best.K, best.nm, npass, Ltried, toc(t0));
    else
        fprintf(['R=%3d (%6d nodes): NO degree passes the gate  [%.0f s]' newline], R0, N0, toc(t0));
    end
    save(fullfile(DAT,'sph_range.mat'), 'RES','RL','LMAX','REVAL','SPAN','s','Ne');
end
fprintf(['done' newline]);
