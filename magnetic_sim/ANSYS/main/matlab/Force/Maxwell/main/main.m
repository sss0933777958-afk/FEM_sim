%% main.m -- Force calibration driver (self-consistency check, either base)
%  load_flux_calib (flux-side parameters) -> build_L -> base_model (generate
%  f_m, expanded form) -> fitting_force (minimise, merged form) -> solve_force
%  (gauge) -> compare against the truth -> save .mat -> emit_force_results (PDF).
%
%  WHAT THIS RUN DOES.  The "measurement" f_m is generated from the flux
%  calibration itself, then every parameter is thrown away and refitted from a
%  deliberately naive starting point. Recovering (l, e, B g, Mbar) is the pass
%  criterion. Note this validates the cost function, the optimiser and the
%  identifiability of the parameters -- it does NOT validate the L formula,
%  because both sides share build_L and an error there would cancel. L was
%  checked separately against finite differences of grad(b.b).
%
%  THE TWO BASES.  BASE = 'current' | 'voltage' is a pure relabelling of the
%  same model; nothing in the arithmetic changes:
%
%      f_ijk = c * (Mbar*u_j)' * L_k * (Mbar*u_j)          expanded (measurement)
%      f_ijk =     (H_hat*u_j)' * L_k * (H_hat*u_j)        merged   (model)
%
%      current   Mbar = K_I bar,  u = Fmap [A],   c = F g_I [pN/A^2]
%      voltage   Mbar = D bar,    u = V    [mV],  c = F g_V [pN/mV^2]
%
%  and H_hat = sqrt(F g) * Mbar in both. The excitation u is taken from the
%  flux calibration (Fmap or V), not assumed here.
%
%  UNITS.  This package works in um / mT / pN throughout, per the project unit
%  reference; the excitation is A (current) or mV (voltage). The flux .mat
%  stores l_hat in metres; load_flux_calib converts it on the way in.
%
%  Edit the per-run block below to switch base / model / source calibration.

% Batch hook, same contract as the flux main.m: the caller may set MAIN_OVR
% (struct) and its fields overwrite the per-run block below.
if exist('MAIN_OVR','var'), MAIN_OVR_ = MAIN_OVR; else, MAIN_OVR_ = struct(); end
clearvars -except MAIN_OVR_;  clc;

%% ---- per-run knobs ---------------------------------------------------------
% [ADDED 2026-09-21 user decision] Where the measurement f_m comes from:
%   'selfcheck'  f_m is SYNTHESISED from a flux calibration (the original mode,
%                everything below unchanged) -> every parameter has a truth to
%                compare against, so the run ends in a PASS/FAIL table.
%   'sph'        f_m is the REAL measurement, f = (1/2)*mgB*grad(b.b), with b.b
%                taken from a spherical-harmonic model of the FEM field. There
%                is NO truth to compare against; the run instead walks the
%                axes-shell ladder and stops on the flux-side convergence
%                criterion applied to [l_hat, F g].
FSRC     = 'sph';       % 'selfcheck' | 'sph'
BASE     = 'current';   % 'current' | 'voltage'
MODEL    = 'long2016_hexapole_halfcut';
GEOM     = 'tip40um';
MATNAME  = '';          % flux-side source; '' = the per-base default below
% ---- FSRC = 'sph' only ------------------------------------------------------
VARIANT  = 'maxwell';            % [2026-09-22] baseline export
% Degree of the harmonic model. It is PER DATASET, not a constant: the ceiling
% is set by the sph_strat smoothness gate (d3(b.b)/ds^3 keeping one sign on the
% three actuator axes), and the two exports do not agree on it. Measured
% 2026-09-21 on the same 1771 nodes:
%     baseline  L <= 9 passes,  NMAE_all 0.1067 % at L = 9
%     b_0.02    L <= 7 passes,  NMAE_all 0.7406 % at L = 7
% The b_0.02 figure is a floor, not a fit failure: it does not fall with L
% (0.744 % at L=5, 0.697 % at L=25), i.e. that export carries ~0.7 % of field
% that no harmonic can represent. Set LDEG to the gate ceiling of whatever
% VARIANT is selected above.
LDEG     = 9;           % baseline: 9   |   b_0.02: 7
% Excitation set. 'pairs21' = the six unit vectors plus all C(6,2) = 15 pairs.
% The field superposes (linear material), so every pair is read off the same six
% fitted harmonic coefficient sets -- no extra fit, no extra solve. The pairs
% earn their place on the FORCE side: f is quadratic in H_hat*u, so a pair
% carries the cross term 2*h_j' L h_k and fixes the RELATIVE column signs of
% H_hat, which a one-at-a-time set leaves undetermined.
EXC      = 'pairs21';   % 'pairs21' | 'singles6'

mgB      = 0.0451;      % particle magnetisation constant [A.um^2/mT], from m = mgB*b
UF       = 1e3;         % force unit factor: A.um.mT = 1e-9 N = 1e3 pN, so UF puts f in pN

R_select = 150;         % sampling sphere radius [um]
NR_AXSH  = [];          % axes-shell layers; [] = inherit from the source calib
USE_BIAS = [];          % [] = inherit from the source calib

% Initial guess (user decision, 2026-08-30): H_hat = I, l = l_design, e = 0.
% H0_SCALE multiplies that identity. It exists because the natural size of
% H_hat = sqrt(F g)*Mbar differs by orders of magnitude between the bases
% (current ~1, voltage ~1e-3), and the fit is known to be far more sensitive to
% the scale of H0 than to the number of unknowns. 1 keeps the note's guess.
H0_SCALE = 1;
l0       = 500;         % l_design [um];  e0 = 0 is built into fitting_force

TOL_PCT  = 0.1;         % self-consistency pass threshold [%]  (FSRC='selfcheck')

% ---- convergence criterion, FSRC = 'sph' only -------------------------------
% Same two-stage judge as the flux main.m, same numbers (user's call: do not
% touch the criterion). Series = [l_hat, F g]; both must converge and the later
% index wins. F g is used rather than B g because it is what the gauge gives
% directly -- B g needs mgB, and the criterion should not lean on an outside
% constant. The K_I bar style physical gate is NOT applied here (user's call).
%   (1) steady value  : first window of CONV_KS levels whose step-to-step
%                       relative change is below CONV_TOL_S %  -> value at the
%                       window head
%   (2) convergence   : first window of CONV_KC levels with at most CONV_KOUT
%                       of them outside steady*(1 +/- CONV_TOL_C %) -> the head
%                       of that window is N_c
CONV_TOL_S = 0.01;      % [%]
CONV_KS    = 20;
CONV_TOL_C = 0.2;       % [%]
CONV_KC    = 10;
CONV_KOUT  = 1;
% Safety stop, NOT a design limit. Lower than the flux side's 800 because a
% force level costs ~2 s (21 excitations x 54 unknowns) and Nr = 120 already
% means 721 sample points.
CONV_NDMAX = 120;

%% ---- apply batch overrides -------------------------------------------------
ovf_ = fieldnames(MAIN_OVR_);
for k_ = 1:numel(ovf_)
    eval([ovf_{k_} ' = MAIN_OVR_.(ovf_{k_});']);
    fprintf(['[MAIN_OVR] ' ovf_{k_} newline]);
end

BASE = validatestring(BASE, {'current','voltage'}, 'main', 'BASE');
if isempty(MATNAME)
    switch BASE
        case 'current', MATNAME = 'calib_current_maxwell_axshN13_R150_eighteen';
        case 'voltage', MATNAME = 'calib_voltage_maxwell_axshN13_soff3mm_R150_eighteen';
    end
end
if isempty(mgB)
    error('main:mgB', ['mgB is not set: give the particle magnetisation ' ...
                       'constant [A.um^2/mT] in the per-run block above']);
end
validateattributes(mgB, {'numeric'}, {'scalar','real','positive','finite'}, 'main', 'mgB');
H0 = H0_SCALE * eye(6);

%% ---- paths (self-relative) -------------------------------------------------
here = fileparts(mfilename('fullpath'));      % .../Force/Maxwell/main
FMX  = fileparts(here);                       % .../Force/Maxwell
addpath(fullfile(FMX, 'function'));

SELF = strcmpi(FSRC, 'selfcheck');
if ~SELF, assert(strcmpi(FSRC,'sph'), 'main:fsrc', 'FSRC is ''selfcheck'' or ''sph'''); end

if SELF
%% ---- (1) flux-side parameters, measurement naming --------------------------
cal = load_flux_calib(MODEL, GEOM, MATNAME, BASE);   % also puts Flux/function on the path
if isempty(USE_BIAS), USE_BIAS = logical(cal.meta.USE_BIAS); end
u   = cal.u_m;                                % excitation, column j is u_j

%% ---- (2) sampling points, actuator frame, um -------------------------------
if isempty(NR_AXSH)
    assert(strcmp(cal.meta.sampler, 'axes_shells'), 'main:sampler', ...
           ['source calib used sampler ''%s''; NR_AXSH cannot be inherited, ' ...
            'set it explicitly'], cal.meta.sampler);
    NR_AXSH = (cal.meta.npts - 1) / 6;
    assert(NR_AXSH == round(NR_AXSH) && NR_AXSH >= 1, 'main:nrAxsh', ...
           'npts=%g is not 6*Nr+1', cal.meta.npts);
end
% conv_design_ws returns [P, Bstack, info]; the info struct is the third output.
[~, ~, si] = conv_design_ws(NR_AXSH, R_select, struct('points_only',true,'R_act',cal.R_act,'quiet',true));
P    = si.P_act;                              % actuator frame, um
npts = size(P, 1);
fprintf('[ws] axes-shell Nr=%d -> %d points, r = %.1f ... %.1f um%s', ...
        NR_AXSH, npts, si.h, R_select, newline);

%% ---- (3) generate the measurement (EXPANDED form) --------------------------
% c_m IS F g,m -- the whole prefactor. The square root that turns it into the
% merged form lives on H_hat = sqrt(F g)*Mbar, not here.
c_m  = UF * (mgB / (2 * cal.l_m)) * cal.gB_m^2;
L_m  = build_L(P, cal.l_m, cal.e_m, cal.Pc_base);
f_m  = base_model(L_m, cal.Mbar_m, u, c_m);
fmax = max(abs(f_m), [], 'all');
fprintf('[f_m] c_m = F g,m = %.6g   peak |f_m| = %.4g pN   size=%s%s', ...
        c_m, fmax, mat2str(size(f_m)), newline);

%% ---- (4) minimise (MERGED form: c = 1, H_hat is the unknown) ---------------
[l_hat, e_hat, H_hat, J, fout] = fitting_force(f_m, P, cal.Pc_base, l0, H0, u, USE_BIAS);
fprintf('[fit] %d unknowns, J=%.4e, rms=%.4e pN, %.1f s%s', ...
        fout.nvar, J, fout.rms, fout.elapsed, newline);

%% ---- (5) gauge to F g and Mbar, then derive B g ----------------------------
[Mbar_hat, gF_hat, gB_hat, sout] = solve_force(H_hat, mgB, l_hat, UF, BASE);
lb = sout.labels;
if sout.nflip > 0
    fprintf('[gauge] flipped column(s) %s to enforce diag(%s) > 0%s', ...
            mat2str(sout.flipped), lb.M_name, newline);
end

%% ---- (6) self-consistency verdict ------------------------------------------
d.l  = abs(l_hat  - cal.l_m ) / abs(cal.l_m ) * 100;
d.gF = abs(gF_hat - c_m     ) / abs(c_m     ) * 100;   % the fitted gain
d.gB = abs(gB_hat - cal.gB_m) / abs(cal.gB_m) * 100;   % the derived one
d.M  = norm(Mbar_hat - cal.Mbar_m, 'fro') / norm(cal.Mbar_m, 'fro') * 100;
d.e  = max(abs(e_hat - cal.e_m));                      % dimensionless, absolute

fprintf('%s---- self-consistency, base=%s (tol %.3g%%) ----%s', newline, BASE, TOL_PCT, newline);
fprintf('  l      %12.4f um  vs %12.4f um   d = %8.3g %%   %s%s', ...
        l_hat, cal.l_m, d.l, verdict(d.l < TOL_PCT), newline);
fprintf('  %-6s %12.6g     vs %12.6g      d = %8.3g %%   %s   [%s, fitted]%s', ...
        lb.gF_name, gF_hat, c_m, d.gF, verdict(d.gF < TOL_PCT), lb.gF_unit, newline);
fprintf('  %-6s %12.6g     vs %12.6g      d = %8.3g %%   %s   [%s, derived via mgB]%s', ...
        lb.gB_name, gB_hat, cal.gB_m, d.gB, verdict(d.gB < TOL_PCT), lb.gB_unit, newline);
fprintf('  %-6s Frobenius                              d = %8.3g %%   %s%s', ...
        lb.M_name, d.M, verdict(d.M < TOL_PCT), newline);
if USE_BIAS
    fprintf('  e      max abs deviation                    d = %8.3g     %s%s', ...
            d.e, verdict(d.e < 1e-4), newline);
end
fprintf('  %s checks: diag>0 %d  diag-dominant %d  offdiag<0 %d  max|rowsum| %.3g%s', ...
        lb.M_name, sout.chk.diag_pos, sout.chk.diag_dominant, sout.chk.offdiag_neg, ...
        sout.chk.rowsum_max, newline);

else
%% ---- (1s) FEM field -> ONE harmonic model ----------------------------------
% The coefficients are fitted once, on the six single-pole exports. Every one of
% the 21 excitations is then read off the same fit: sph_field's grad(Pq, I, ...)
% forms c(I) = sum_k I_k c_k internally, which IS the superposition. No refit.
here_ = fileparts(mfilename('fullpath'));
FMX_  = fileparts(here_);
MAINR = fileparts(fileparts(fileparts(FMX_)));          % .../ANSYS/main
CALR  = fullfile(MAINR, 'matlab', 'Flux', 'Maxwell');
addpath(fullfile(CALR,'function'), fullfile(CALR,'utils'), fullfile(CALR,'common_path'));
addpath(fullfile(FMX_, 'utils'));                       % sph_field

cfg = model_config(MODEL, GEOM);
% The model superposes RAW per-ampere fields, so the all-source presentation
% flip must be a no-op for this dataset. Assert rather than assume.
assert(all(cfg.s_source == 1), 'main:sSource', ...
       ['cfg.s_source = %s: raw is not all-source, so b = sum_k I_k b_k would ' ...
        'superpose sign-flipped fields'], mat2str(cfg.s_source));
raw = extract_maxwell_data(cfg, 'all', VARIANT);
ad  = build_actuator_data(raw, cfg);                    % actuator frame, mT, iron removed
assert(isequal(ad.F, eye(6)), 'main:Fmap', ...
       'the coil->pole map is not identity here; u would address the wrong poles');
inb = find(ad.r2 < (R_select*1e-6)^2);
Pn  = ad.Pa(inb,:) * 1e6;                               % um, actuator frame
Bn  = ad.Ba(inb,:,:);                                   % mT, page k = pole Pk
[~, gq, sinfo] = sph_field(R_select, struct('P',Pn,'B',Bn), ...
                           struct('L',LDEG,'quiet',true,'condmax',1500));
fprintf('[sph] L=%d K=%d on %d nodes (r<=%g um), NMAE_all %.4f %%%s', ...
        sinfo.L, sinfo.K, numel(inb), R_select, sinfo.NMAE_all, newline);

switch lower(EXC)
    case 'singles6', u = eye(6);
    case 'pairs21'
        pr = nchoosek(1:6, 2);
        u  = [eye(6), zeros(6, size(pr,1))];
        for q = 1:size(pr,1), u(pr(q,:), 6+q) = 1; end
    otherwise, error('main:exc', 'EXC is ''pairs21'' or ''singles6''');
end
if isempty(USE_BIAS), USE_BIAS = true; end
fprintf('[exc] u = %s  (%s)%s', mat2str(size(u)), EXC, newline);

% `cal` is stubbed so the record written in (7) keeps one shape across both
% modes. There is no truth in this mode: every *_m field is NaN.
cal = struct('Pc_base',ad.Pc_base, 'R_act',ad.R_act, 'u_m',u, ...
             'l_m',NaN, 'e_m',nan(17,1), 'gB_m',NaN, 'Mbar_m',nan(6,6), ...
             'src',sprintf('%s | %s', VARIANT, cfg.fld_dir), ...
             'meta',struct('VARIANT',VARIANT, 'SOFF',NaN, 'USE_BIAS',USE_BIAS, ...
                           'sampler','axes_shells', 'npts',NaN));
c_m = NaN;

%% ---- (2s-5s) walk the ladder, stop on the flux convergence criterion -------
if isempty(NR_AXSH)
    [~,~,wi] = conv_design_ws([], R_select, struct('ladder',CONV_NDMAX));
    LADW = wi.ladder;                                   % Nr = 1, 2, 3, ...
else
    LADW = NR_AXSH(:);                                  % fixed design, no loop
end
SER = nan(numel(LADW), 2);   nc_hit = numel(LADW) == 1;   t_all = tic;
for q = 1:numel(LADW)
    [l_hat, e_hat, H_hat, J, fout, P, npts, Mbar_hat, gF_hat, gB_hat, sout] = ...
        run_level(LADW(q), R_select, cal, gq, u, mgB, UF, l0, H0, USE_BIAS, BASE);
    SER(q,:) = [l_hat, gF_hat];
    fprintf('  Nr=%3d  N=%4d   l_hat %9.3f um   %s %10.6f   rms %9.3g pN%s', ...
            LADW(q), npts, l_hat, sout.labels.gF_name, gF_hat, fout.rms, newline);
    if nc_hit, break; end
    iE = judge(SER(1:q,1), CONV_TOL_S, CONV_KS, CONV_TOL_C, CONV_KC, CONV_KOUT);
    iG = judge(SER(1:q,2), CONV_TOL_S, CONV_KS, CONV_TOL_C, CONV_KC, CONV_KOUT);
    % BOTH series must converge. max() would silently ignore a NaN and stop the
    % ladder on l_hat alone -- the trap already documented in the flux main.m.
    if ~all(isfinite([iE iG])), continue; end
    NR_AXSH = LADW(max(iE, iG));                        % N_c is the judge window HEAD
    [l_hat, e_hat, H_hat, J, fout, P, npts, Mbar_hat, gF_hat, gB_hat, sout] = ...
        run_level(NR_AXSH, R_select, cal, gq, u, mgB, UF, l0, H0, USE_BIAS, BASE);
    fprintf(['[conv] l_hat converged at level %d, gain at level %d -> ' ...
             'N_c: Nr=%d, %d points (swept %d levels, %.0f s)%s'], ...
            iE, iG, NR_AXSH, npts, q, toc(t_all), newline);
    nc_hit = true;   break
end
if ~nc_hit
    error('main:notConverged', ...
          'no convergence within %d ladder levels (last Nr=%d)', CONV_NDMAX, LADW(end));
end
fmax = fout.fmax;
lb   = sout.labels;
d.l = NaN;  d.gF = NaN;  d.gB = NaN;  d.M = NaN;  d.e = NaN;   % no truth in this mode

fprintf('%s---- calibration, base=%s, source=%s, L=%d ----%s', ...
        newline, BASE, VARIANT, LDEG, newline);
fprintf('  l_hat  %12.4f um%s', l_hat, newline);
fprintf('  %-6s %12.6g      [%s, fitted]%s', lb.gF_name, gF_hat, lb.gF_unit, newline);
fprintf('  %-6s %12.6g      [%s, derived via mgB]%s', lb.gB_name, gB_hat, lb.gB_unit, newline);
fprintf('  rms residual %.4g pN   peak |f| %.4g pN   (%d points x %d excitations)%s', ...
        fout.rms, fmax, npts, size(u,2), newline);
fprintf('  %s checks: diag>0 %d  diag-dominant %d  offdiag<0 %d  max|rowsum| %.3g%s', ...
        lb.M_name, sout.chk.diag_pos, sout.chk.diag_dominant, sout.chk.offdiag_neg, ...
        sout.chk.rowsum_max, newline);
end

%% ---- (7) save (self-describing) --------------------------------------------
rec = struct('base',BASE, ...
    'l_hat',l_hat, 'e_hat',e_hat, 'H_hat',H_hat, 'Mbar_hat',Mbar_hat, ...
    'gF_hat',gF_hat, 'gB_hat',gB_hat, ...
    'J',J, 'rms',fout.rms, 'nvar',fout.nvar, 'elapsed',fout.elapsed, 'fmax',fmax, ...
    'Jhist',fout.Jhist, ...
    'l_m',cal.l_m, 'e_m',cal.e_m, 'gB_m',cal.gB_m, 'gF_m',c_m, 'Mbar_m',cal.Mbar_m, ...
    'c_m',c_m, 'scale',sout.scale, ...
    'd_l_pct',d.l, 'd_gF_pct',d.gF, 'd_gB_pct',d.gB, 'd_M_pct',d.M, 'd_e_abs',d.e, ...
    'mgB',mgB, 'UF',UF, 'model',MODEL, 'GEOM',GEOM, 'MATNAME',MATNAME, 'src',cal.src, ...
    'VARIANT',cal.meta.VARIANT, 'SOFF',cal.meta.SOFF, ...
    'R_select',R_select, 'NR_AXSH',NR_AXSH, 'npts',npts, 'USE_BIAS',USE_BIAS, ...
    'u',u, 'l0',l0, 'H0',H0, 'H0_SCALE',H0_SCALE, 'Pc_base',cal.Pc_base, 'P',P, ...
    'flipped',sout.flipped, 'h11',sout.h11, 'chk',sout.chk, 'labels',lb, ...
    'FSRC',FSRC);

% [ADDED 2026-09-21] the sph mode carries what the selfcheck mode has no use for:
% which harmonic model produced the measurement, and the whole ladder series, so
% the convergence can be re-plotted without re-running anything.
if ~SELF
    rec.MATNAME = '';         % no flux .mat was read in this mode
    rec.LDEG   = LDEG;        rec.EXC  = EXC;
    rec.NMAE_sph = sinfo.NMAE_all;   rec.K_sph = sinfo.K;
    rec.ladder = LADW(:);     rec.SER  = SER;    % SER(:,1)=l_hat, SER(:,2)=F g
    rec.conv   = struct('TOL_S',CONV_TOL_S, 'KS',CONV_KS, 'TOL_C',CONV_TOL_C, ...
                        'KC',CONV_KC, 'KOUT',CONV_KOUT, 'NDMAX',CONV_NDMAX);
end

% The stem is built once, here, and carried in the record: emit_force_results
% names the PDF from rec.stem rather than re-deriving it, so the two can never
% drift apart.
rec.stem = force_stem(rec);
tag      = 'single';  if USE_BIAS, tag = 'eighteen'; end
matdir   = fullfile(FMX, 'data', MODEL, '.mat');
if ~exist(matdir, 'dir'), mkdir(matdir); end
matfile = fullfile(matdir, [rec.stem '_' tag '.mat']);
save(matfile, '-struct', 'rec');
fprintf('saved %s%s', matfile, newline);

%% ---- (8) PDF ---------------------------------------------------------------
% emit_force_results lays out the self-consistency comparison (it prints l_m,
% B g,m, Mbar_m next to the fitted ones). In sph mode those truths do not exist,
% so the PDF is left for a later, purpose-written emitter rather than fed NaNs.
if SELF
    emit_force_results(matfile);
else
    fprintf('[pdf] skipped: emit_force_results is the self-consistency layout%s', newline);
end

%% ---- local -----------------------------------------------------------------
function s = verdict(ok)
    if ok, s = 'PASS'; else, s = '<<< FAIL'; end
end

function [l_hat, e_hat, H_hat, J, fout, P, npts, Mbar_hat, gF_hat, gB_hat, sout] = ...
         run_level(Nr, R_select, cal, gq, u, mgB, UF, l0, H0, USE_BIAS, BASE)
%RUN_LEVEL  One ladder level: sample points -> measurement -> fit -> gauge.
%   The measurement is the real one:  f = (1/2) * mgB * grad(b.b), with grad
%   taken analytically from the harmonic model (sph_field's 'du' returns
%   d(b.b)/dx_j). UF puts it in pN.
    [~, ~, si] = conv_design_ws(Nr, R_select, ...
                    struct('points_only',true, 'R_act',cal.R_act, 'quiet',true));
    P    = si.P_act;                  % um, actuator frame
    npts = size(P, 1);
    NI   = size(u, 2);
    f_m  = zeros(3, npts, NI);
    for j = 1:NI
        f_m(:,:,j) = (0.5 * mgB * UF * gq(P, u(:,j), 'du')).';
    end
    [l_hat, e_hat, H_hat, J, fout] = fitting_force(f_m, P, cal.Pc_base, l0, H0, u, USE_BIAS);
    fout.fmax = max(abs(f_m), [], 'all');
    [Mbar_hat, gF_hat, gB_hat, sout] = solve_force(H_hat, mgB, l_hat, UF, BASE);
end

function i0 = judge(v, tolS, KS, tolC, KC, kout)
%JUDGE  Two-stage convergence criterion, the same ruler as the flux main.m
%   (that copy returns a ladder index too; the ladder here is Nr = 1,2,3,... so
%   the index IS Nr). Kept as a local on purpose: the two branches each hold
%   their own copy rather than sharing an engine file (user's call, 2026-08-24).
%     (1) steady value vs : first window of KS consecutive step-to-step changes
%                           all below tolS %, take the value at the window head
%     (2) convergence     : first window of KC levels with at most `kout` of
%                           them outside vs*(1 +/- tolC %), return its head
%   tolS / tolC are PERCENTAGES (0.01 means 0.01 %, not 1 %).
    if nargin < 6 || isempty(kout), kout = 0; end
    i0 = NaN;   v = v(:);
    if nnz(isfinite(v)) < KS + 1, return; end
    ch = abs(diff(v)) ./ abs(v(1:end-1)) * 100;
    vs = NaN;
    for i = 1:(numel(ch)-KS+1)
        w = ch(i:i+KS-1);
        if all(isfinite(w)) && all(w < tolS), vs = v(i+1);  break; end
    end
    if ~isfinite(vs), return; end
    inb = abs(v - vs)/abs(vs)*100 <= tolC;
    for i = 1:(numel(v)-KC+1)
        if nnz(~inb(i:i+KC-1)) <= kout, i0 = i;  return; end
    end
end

% Output stem:  <base>_R<R>_N<npts>[_<variant tag>][_soff<d>mm]
%   The model and single/eighteen are expressed by the folders, so they stay out
%   of the name (short-names rule). The sensor offset IS in the name for a
%   voltage run: the flux side keeps several offsets side by side and an
%   untagged name would let them overwrite each other.
function s = force_stem(rec)
    v = regexprep(rec.VARIANT, '^maxwell_?', '');     % drop the branch name
    v = regexprep(v, '_?axshN\d+', '');               % point count goes in _N<npts>
    if ~isempty(v), v = ['_' v]; end
    so = '';
    if strcmp(rec.base, 'voltage') && isfinite(rec.SOFF)
        so = ['_soff' strrep(sprintf('%gmm', rec.SOFF*1e3), '.', 'p')];
    end
    % [ADDED 2026-09-21] the harmonic degree goes in the name for an sph run:
    % it is the measurement model, so two runs that differ only in L are two
    % different calibrations and must not overwrite each other.
    lt = '';
    if isfield(rec,'LDEG') && ~isempty(rec.LDEG), lt = sprintf('_L%d', rec.LDEG); end
    s = sprintf('%s_R%d_N%d%s%s%s', rec.base, round(rec.R_select), rec.npts, v, so, lt);
end
