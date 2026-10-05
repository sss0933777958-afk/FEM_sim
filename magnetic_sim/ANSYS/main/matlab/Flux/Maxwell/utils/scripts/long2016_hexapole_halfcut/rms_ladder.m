%% rms_ladder.m -- calibration RMS vs number of calibration points (long2016, current)
%  Question: how does the calibration error RMS = sqrt(J/(3*6*N)) change along the
%  axes-shell ladder, for the single-parameter and the eighteen-parameter model?
%  Input : Maxwell baseline field (model_config long2016_hexapole_halfcut / tip40um, variant maxwell)
%  Output: RLIST = []     -> utils/data/<MODEL>/rms_ladder.mat    (one R: Nr, N, J, RMS, l_hat for both models)
%          RLIST = vector -> utils/data/<MODEL>/rms_ladder_R.mat  (per R: Nc, RMS_inf, l_hat, g_I at Nc)
%  Plot  : plot/long2016_hexapole_halfcut/current/rms_ladder.m          (rms_ladder.mat)
%          figures/paper_fig_plot/plot/plot_ell_gain_2panel.m           (rms_ladder_R.mat)
%
%  Convergence (R sweep):
%    RMS_inf = RMS at the start of the first window of KW consecutive step changes < TOLS %
%    Nc      = start of the first window of KC consecutive rungs with RMS - RMS_inf <= TOLC * RMS_inf
%  Each R's ladder is extended until RMS_inf is found for both models (cap NRMAX).

MODEL   = 'long2016_hexapole_halfcut';
GEOM    = 'tip40um';
VARIANT = 'maxwell';
R       = 150e-6;          % sampling sphere radius [m] (single-R mode)
l0      = 0.5e-3;          % l_hat initial guess [m]
NR      = 1:400;           % N = 6*Nr + 1 -> 7 ... 2401 (single-R mode)
RLIST   = (40:20:500)*1e-6;          % [] = single-R ladder; e.g. (40:20:500)*1e-6 = R sweep
TOLS    = 0.01;            % RMS_inf: step-change threshold [%]
KW      = 10;              % RMS_inf: consecutive steps
TOLC    = 0.10;            % Nc: RMS within 10 % of RMS_inf
KC      = 10;              % Nc: consecutive rungs inside the band
NRMAX   = 1000;            % R sweep: ladder cap

here = fileparts(mfilename('fullpath'));            % .../utils/scripts/<model>
CALROOT = fileparts(fileparts(fileparts(here)));       % .../Flux/Maxwell
addpath(fullfile(CALROOT, 'function'), fullfile(CALROOT, 'common_path'));

cfg = model_config(MODEL, GEOM);
raw = extract_maxwell_data(cfg, 'all', VARIANT);
ad  = build_actuator_data(raw, cfg);
opt = struct('model',MODEL, 'geom',GEOM, 'variant',VARIANT, 'frame','actuator', 'raw',raw);

if isempty(RLIST)
    %% ---- single R: full ladder ---------------------------------------------
    nq = numel(NR);
    N  = nan(nq,1);   J = nan(nq,2);   l_hat = nan(nq,2);   % columns: single, eighteen
    BRMS = nan(nq,1);                                    % field RMS at the calibration points [mT]
    t0 = tic;
    for q = 1:nq
        [P, Bs, gi] = conv_design_ws(NR(q), R, opt);
        N(q) = gi.npts_kept;
        BRMS(q) = sqrt(sum(Bs(:).^2) / (3*6*N(q)));      % same points and divisor as RMS
        for m = 1:2
            [~, l_hat(q,m), J(q,m)] = fitting(P, Bs, ad.Pc_base, l0, m == 2);
        end
        fprintf('Nr=%3d N=%3d | RMS single %.5f  eighteen %.5f mT | %.0f s\n', ...
                NR(q), N(q), sqrt(J(q,1)/(18*N(q))), sqrt(J(q,2)/(18*N(q))), toc(t0));
    end
    RMS = sqrt(J ./ (3*6*N));      % [mT] per residual: 3 components x 6 excitations x N points

    out = fullfile(CALROOT, 'utils', 'data', MODEL, 'rms_ladder.mat');
    Nr = NR(:);
    save(out, 'MODEL','GEOM','VARIANT','R','l0','Nr','N','J','RMS','l_hat','BRMS');
    fprintf('saved %s\n', out);
    return
end

%% ---- R sweep: Nc per R, then l_hat and g_I at Nc ---------------------------
F = zeros(6, cfg.N_I);
for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end
ws = warning('off', 'MATLAB:nearlySingularMatrix');

nR = numel(RLIST);
Nc = nan(nR,2);   Nr_c = nan(nR,2);   RMS_inf = nan(nR,2);   N_inf = nan(nR,2);
l_hat = nan(nR,2);   gI = nan(nR,2);   KI = nan(6,6,nR,2);   e_hat = cell(nR,2);
lad = cell(nR,1);                                    % per-R ladder: [Nr N RMS_single RMS_eighteen]
out = fullfile(CALROOT, 'utils', 'data', MODEL, 'rms_ladder_R.mat');
t0 = tic;
for k = 1:nR
    Rk = RLIST(k);
    [P_ev, B_ev] = cfg.select_ball(ad, Rk);
    L = nan(NRMAX, 4);   iw = nan(1,2);
    for q = 1:NRMAX
        evalc('[P, Bs, gi] = conv_design_ws(q, Rk, opt);');    % silence the iron-exclusion report
        L(q,1:2) = [q, gi.npts_kept];
        for m = 1:2
            if ~isnan(iw(m)), continue; end          % this model already converged
            [~, ~, Jq] = fitting(P, Bs, ad.Pc_base, l0, m == 2);
            L(q,2+m) = sqrt(Jq / (18*gi.npts_kept));
            iw(m) = rms_window(L(1:q,2+m), TOLS, KW);
        end
        if all(isfinite(iw)), break; end
    end
    L = L(1:q,:);   lad{k} = L;
    for m = 1:2
        if isnan(iw(m)), continue; end
        r = L(1:iw(m)+KW, 2+m);   % rungs up to the window end are all filled
        RMS_inf(k,m) = r(iw(m));   N_inf(k,m) = L(iw(m),2);
        inb = r - RMS_inf(k,m) <= TOLC*RMS_inf(k,m);
        ic  = find(movsum(double(inb), [0 KC-1], 'Endpoints','discard') == KC, 1);
        Nr_c(k,m) = L(ic,1);   Nc(k,m) = L(ic,2);
        evalc('[P, Bs] = conv_design_ws(Nr_c(k,m), Rk, opt);');
        [e, lh] = fitting(P, Bs, ad.Pc_base, l0, m == 2);
        [KI(:,:,k,m), gI(k,m)] = solve_current(lh, e, ad.Pc_base, P, Bs, F, [], P_ev, B_ev);
        l_hat(k,m) = lh;   e_hat{k,m} = e;
    end
    fprintf(['R=%3.0f um | single Nc=%4d l=%.1f g=%.3f | eighteen Nc=%4d l=%.1f g=%.3f ' ...
             '| ladder %d | %.0f s\n'], Rk*1e6, Nc(k,1), l_hat(k,1)*1e6, gI(k,1), ...
            Nc(k,2), l_hat(k,2)*1e6, gI(k,2), q, toc(t0));
    R_um = RLIST(:)*1e6;                               % save after every R (partial results kept)
    save(out, 'MODEL','GEOM','VARIANT','l0','TOLS','KW','TOLC','KC','NRMAX', ...
         'R_um','Nc','Nr_c','RMS_inf','N_inf','l_hat','gI','KI','e_hat','lad');
end
warning(ws);
fprintf('saved %s\n', out);

% ---------------------------------------------------------------------------
function i0 = rms_window(r, tol, K)
% First i such that the K step changes r(i)->r(i+1) ... r(i+K-1)->r(i+K) are all < tol %.
    i0 = NaN;
    if numel(r) < K + 1, return; end
    ch = abs(diff(r)) ./ r(1:end-1) * 100;
    for i = 1:numel(ch)-K+1
        if all(ch(i:i+K-1) < tol), i0 = i;  return; end
    end
end
