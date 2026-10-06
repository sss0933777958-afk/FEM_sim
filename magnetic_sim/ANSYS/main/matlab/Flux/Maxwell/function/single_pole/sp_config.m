function cfg = sp_config(model, geom)
%SP_CONFIG  Load config/<model>/<geom>/mt_constants.m for a single-pole model.
%   cfg = SP_CONFIG(model, geom)
%   Single-pole counterpart of hexapole/model_config.m. Differences:
%     - keeps the config's N_I (hexapole/model_config forces N_I = 6)
%     - requires strategy = 'single_pole'
%   Named sp_* so it never shadows the hexapole functions on the MATLAB path.

    here = fileparts(mfilename('fullpath'));              % .../function/single_pole
    CAL  = fileparts(fileparts(here));                    % .../Maxwell
    cfgdir = fullfile(CAL, 'config', model, geom);
    if ~isfile(fullfile(cfgdir, 'mt_constants.m'))
        error('sp_config:noConfig', 'no single-pole config at %s', cfgdir);
    end
    addpath(cfgdir, '-begin');
    clear mt_constants;                                   % drop a cached mt_constants from another config
    cfg = mt_constants();
    rmpath(cfgdir);

    cfg.model = model;
    cfg.geom  = geom;
    need = {'strategy','N_I','fld_dir','fld_files','WP','s_source','calib'};
    for k = 1:numel(need)
        assert(isfield(cfg, need{k}), 'config/%s/%s mt_constants.m lacks field ''%s''', model, geom, need{k});
    end
    assert(strcmp(cfg.strategy, 'single_pole'), ...
           'sp_config: strategy is ''%s'', expected ''single_pole''', cfg.strategy);
end
