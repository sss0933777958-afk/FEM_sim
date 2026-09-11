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
BASE     = 'current';   % 'current' | 'voltage'
MODEL    = 'long2016_hexapole_halfcut';
GEOM     = 'tip40um';
MATNAME  = '';          % flux-side source; '' = the per-base default below

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

TOL_PCT  = 0.1;         % self-consistency pass threshold [%]

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
    'flipped',sout.flipped, 'h11',sout.h11, 'chk',sout.chk, 'labels',lb);

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
emit_force_results(matfile);

%% ---- local -----------------------------------------------------------------
function s = verdict(ok)
    if ok, s = 'PASS'; else, s = '<<< FAIL'; end
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
    s = sprintf('%s_R%d_N%d%s%s', rec.base, round(rec.R_select), rec.npts, v, so);
end
