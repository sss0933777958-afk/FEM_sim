function cal = load_flux_calib(model, geom, matname, base)
%LOAD_FLUX_CALIB  Pull a flux-side calibration result into force-model naming.
%
%   cal = LOAD_FLUX_CALIB(model, geom, matname, base)
%     model   : 'long2016_hexapole_halfcut' | 'zhi_peng' | 'hung_hexapole' | ...
%     geom    : config geometry variant ('tip40um', 'R500', ...); '' for flat configs
%     matname : calib file stem under Flux/Maxwell/data/<model>/.mat/
%               Omit it (or pass '') to list the candidates for `base` and stop.
%     base    : 'current' (default) | 'voltage'
%
%   Returns the MEASUREMENT-side parameters. They carry the _m suffix so they can
%   never be confused with the fitted quantities (_hat) produced later by
%   fitting_force:
%
%     cal.l_m       l_hat   [um]     effective length (converted from the .mat, which stores m)
%     cal.e_m       e       17x1     dimensionless charge-grid offset
%     cal.gB_m      B g     scalar   field gain: B g_I [mT/A] or B g_V [mT/mV]
%     cal.Mbar_m    Mbar    6x6      gauged matrix, Mbar(1,1) = 5/6:
%                                      current -> K_I bar     voltage -> D bar
%     cal.u_m       u       6x6      excitation, column j is u_j:
%                                      current -> Fmap [A]    voltage -> V [mV]
%     cal.base                       which of the two the above are
%     cal.Pc_base           3x6      ideal charge grid, from model_config
%     cal.R_act             3x3      measure -> actuator rotation, for sampling
%     cal.src                        absolute path of the .mat actually read
%     cal.meta                       provenance struct (see below)
%
%   [MODIFIED 2026-09-02] Both bases are accepted. `base` selects which family of
%   files is listed and which fields are read, and the loaded record's own `base`
%   field is checked against it -- reading a voltage calibration as if it were a
%   current one would silently produce a wrong f_m.
%
%   Per .claude/rules/result-read-safety.md the resolved path and a fingerprint
%   are printed before the data is used, and an ambiguous request stops instead
%   of guessing which file was meant.

    if nargin < 2, geom    = ''; end
    if nargin < 3, matname = ''; end
    if nargin < 4 || isempty(base), base = 'current'; end
    base = validatestring(base, {'current','voltage'}, mfilename, 'base');

    % the field names that differ between the two bases; everything else is shared
    if strcmp(base, 'current')
        F = struct('M','KI_bar', 'g','gI_hat', 'u','Fmap');
    else
        F = struct('M','D_bar',  'g','gV_hat', 'u','V');
    end

    % ---- locate the Flux/Maxwell tree relative to this file -----------------
    here  = fileparts(mfilename('fullpath'));        % .../matlab/Force/Maxwell/function
    MROOT = fileparts(fileparts(fileparts(here)));   % .../matlab
    FLUX  = fullfile(MROOT, 'Flux', 'Maxwell');
    assert(isfolder(FLUX), 'load_flux_calib:noFlux', 'Flux/Maxwell not found at %s', FLUX);
    addpath(fullfile(FLUX, 'function'));             % model_config lives here

    datadir = fullfile(FLUX, 'data', model, '.mat');
    assert(isfolder(datadir), 'load_flux_calib:noData', ...
           'no calib folder for model ''%s'' (expected %s)', model, datadir);

    % ---- resolve the file, never guess between candidates -------------------
    cand = dir(fullfile(datadir, sprintf('calib_%s_*.mat', base)));
    if isempty(matname)
        list_candidates(datadir, cand, base);
        error('load_flux_calib:needMatname', ...
              'matname not given -- pick one of the candidates listed above');
    end
    if ~endsWith(matname, '.mat'), matname = [matname '.mat']; end
    src = fullfile(datadir, matname);
    if ~isfile(src)
        list_candidates(datadir, cand, base);
        error('load_flux_calib:noSuchMat', '%s does not exist', src);
    end

    rec = load(src);
    assert(isfield(rec, 'base') && strcmp(rec.base, base), ...
           'load_flux_calib:wrongBase', ...
           '%s is base=''%s'' but ''%s'' was asked for', ...
           matname, getfielddef(rec, 'base', '?'), base);
    for fn = {'l_hat', 'e', F.g, F.M, F.u}
        assert(isfield(rec, fn{1}), 'load_flux_calib:missingField', ...
               '%s has no field ''%s''', matname, fn{1});
    end

    % ---- ideal charge grid from the live config (never from backup/) --------
    cfg = model_config(model, geom);
    if isfield(cfg, 'Pc_base')
        Pc_base = cfg.Pc_base;
    elseif isfield(cfg, 'pole_tip_x_wp')                    % ntu_flat: derive from tips
        tipwp   = [cfg.pole_tip_x_wp; cfg.pole_tip_y_wp; cfg.pole_tip_z_wp];
        Pc_base = tipwp ./ vecnorm(tipwp);
    else
        error('load_flux_calib:noPcBase', ...
              'config for %s/%s provides neither Pc_base nor pole tips', model, geom);
    end
    validateattributes(Pc_base, {'numeric'}, {'size',[3 6],'real','finite'}, mfilename, 'Pc_base');

    % ---- pack, measurement-side naming --------------------------------------
    cal = struct();
    cal.l_m     = rec.l_hat * 1e6;                % .mat stores m; this package works in um
    cal.e_m     = rec.e(:);
    cal.gB_m    = rec.(F.g);
    cal.Mbar_m  = rec.(F.M);
    cal.u_m     = rec.(F.u);
    cal.base    = base;
    cal.Pc_base = Pc_base;
    if isfield(cfg,'R_act'), cal.R_act = cfg.R_act; else, cal.R_act = eye(3); end
    cal.src     = src;
    cal.meta    = struct( ...
        'model',    model, ...
        'geom',     geom, ...
        'VARIANT',  getfielddef(rec, 'VARIANT',  ''), ...
        'DATASET',  getfielddef(rec, 'DATASET',  ''), ...
        'R_select', getfielddef(rec, 'R_select', NaN), ...
        'USE_BIAS', getfielddef(rec, 'USE_BIAS', NaN), ...
        'sampler',  getfielddef(rec, 'sampler',  ''), ...
        'npts',     getfielddef(rec, 'npts',     NaN), ...
        'SOFF',     getfielddef(rec, 'SOFF_upper', NaN), ...
        'NMAE',     getfielddef(rec, 'NMAE',     NaN));

    assert(numel(cal.e_m) == 17, 'load_flux_calib:eSize', ...
           'e must be 17x1 (got %d) -- is this an eighteen-parameter calibration?', numel(cal.e_m));

    % ---- fingerprint (result-read-safety layer 1) ---------------------------
    fprintf('[load_flux_calib] %s\n', src);
    fprintf('    model=%s geom=%s variant=%s\n', model, geom, cal.meta.VARIANT);
    fprintf('    R=%g um  sampler=%s  npts=%g  USE_BIAS=%g\n', ...
            cal.meta.R_select*1e6, cal.meta.sampler, cal.meta.npts, cal.meta.USE_BIAS);
    fprintf('    base=%s   l_m=%.2f um   B g_m=%.6g   Mbar(1,1)=%.4f   max|e_m|=%.4g\n', ...
            base, cal.l_m, cal.gB_m, cal.Mbar_m(1,1), max(abs(cal.e_m)));
end

% ---- list what is available, so the caller can choose --------------------
function list_candidates(datadir, cand, base)
    fprintf('[load_flux_calib] %s-base calibrations in %s:\n', base, datadir);
    if isempty(cand)
        fprintf('    (none)\n');
    else
        [~, ord] = sort([cand.datenum], 'descend');
        NSHOW = 15;                                  % keep the listing readable
        if numel(ord) > NSHOW
            fprintf('    (%d files, showing the %d most recent)\n', numel(ord), NSHOW);
            ord = ord(1:NSHOW);
        end
        for k = ord(:).'
            fprintf('    %s   (%s)\n', cand(k).name, string(datetime(cand(k).datenum,'ConvertFrom','datenum','Format','yyyy-MM-dd HH:mm')));
        end
    end
end

function v = getfielddef(s, fn, dflt)
    if isfield(s, fn), v = s.(fn); else, v = dflt; end
end
