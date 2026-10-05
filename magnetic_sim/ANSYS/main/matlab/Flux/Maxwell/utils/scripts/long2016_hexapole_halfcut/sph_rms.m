%% sph_rms.m -- spherical-harmonic field-fit RMS vs degree L (long2016, current, R <= 150 um)
%  Fit sph_field to every exported node with |p| <= R, all 6 excitations, for L = 1..LMAX, and record
%  RMS = sqrt( sum ||b_model - b_FEM||^2 / (3*6*Np) )   [mT]
%  Convergence: first L of ten consecutive steps whose relative change is < 0.01 %.
%  Output: utils/data/<MODEL>/sph_rms.mat

MODEL = 'long2016_hexapole_halfcut';   GEOM = 'tip40um';   VARIANT = 'maxwell';
R = 150;   LMAX = 40;   TOL = 0.01;   KWIN = 10;               % R [um], TOL [%]

here = fileparts(mfilename('fullpath'));   CAL = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
addpath(fullfile(CAL, 'function'), fullfile(CAL, 'common_path'), ...
        fullfile(fileparts(fileparts(CAL)), 'Force', 'Maxwell', 'function'));  % sph_field
cfg = model_config(MODEL, GEOM);
raw = extract_maxwell_data(cfg, 'all', VARIANT);
ad  = build_actuator_data(raw, cfg);                          % actuator frame, mT, all-source
FEM = struct('P', ad.Pa * 1e6, 'B', ad.Ba);                   % [um], [mT]

L   = (1:LMAX).';   RMS = nan(LMAX, 1);   Np = nan;
for q = 1:LMAX
    [~, ~, info] = sph_field(R, FEM, struct('L', L(q), 'quiet', true));
    Np = info.Np;
    RMS(q) = sqrt(sum(info.resid(:).^2) / (3 * 6 * Np));
    fprintf('L=%2d  K=%4d  RMS=%.6e mT\n', L(q), info.K, RMS(q));
end

ch = abs(diff(RMS)) ./ RMS(1:end-1) * 100;                   % step change [%], ch(q): L(q) -> L(q+1)
ok = ch < TOL;   Lc = NaN;
for q = 1:numel(ok) - KWIN + 1
    if all(ok(q:q+KWIN-1)), Lc = L(q); break; end
end
if isnan(Lc)
    fprintf('not converged up to L=%d\n', LMAX);
else
    fprintf('converged: L=%d  RMS=%.6e mT  (window L=%d..%d, Np=%d)\n', Lc, RMS(L == Lc), Lc, Lc + KWIN, Np);
end

out = fullfile(CAL, 'utils', 'data', MODEL, 'sph_rms.mat');
save(out, 'MODEL', 'R', 'L', 'RMS', 'Np', 'TOL', 'KWIN', 'Lc');
fprintf('saved %s\n', out);
