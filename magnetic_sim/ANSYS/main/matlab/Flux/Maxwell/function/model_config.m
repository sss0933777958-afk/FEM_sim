function cfg = model_config(model, geom)
%MODEL_CONFIG  依 (model, geom) 呼叫 config 的 mt_constants.m，疊共用預設後回完整 cfg。
%   cfg = MODEL_CONFIG(model)          % geom='' → 用 flat config/<model>/mt_constants.m（hung/NTU）
%   cfg = MODEL_CONFIG(model, geom)    % 用 config/<model>/<geom>/mt_constants.m（幾何變體，如 tip40um/tip400um）
%   本函式只做 dispatch + 疊共用欄位（N_I=6、model、geom）+ sanity。
%
%   幾何解析：geom 非空且 config/<model>/<geom>/ 存在 → 載該子夾；否則有 flat config/<model>/mt_constants.m
%   就載 flat；兩者皆無（如 long2016 未給 geom）→ 報錯並列出可用 geom 子夾（逼明確指定、不靜默載錯幾何）。
%
%   🔴 R1：mt_constants 是同名函式（各 geom 子夾一份），MATLAB 會快取第一次解析。interp 路徑一個 session 內
%      會載兩個 geom（如 tip400um + tip40um），故載入前 `clear mt_constants` 強制依 '-begin' 重解析，否則拿到舊版。

    if nargin < 2, geom = ''; end
    here = fileparts(mfilename('fullpath'));              % .../function
    CAL  = fileparts(here);                               % .../Calibration_using_FEM_modeling
    base = fullfile(CAL, 'config', model);
    if ~isfolder(base)
        error('model_config:unknownModel', 'no config for model ''%s'' (expected folder %s)', model, base);
    end
    if ~isempty(geom) && isfolder(fullfile(base, geom))
        cfgdir = fullfile(base, geom);
    elseif isfile(fullfile(base, 'mt_constants.m'))
        cfgdir = base;  geom = '';                        % flat（hung/NTU）
    else
        d = dir(base);  subs = {d([d.isdir] & ~startsWith({d.name},'.')).name};
        error('model_config:needGeom', ...
              'model ''%s'' 無 flat config；請指定 geom 變體之一：{%s}', model, strjoin(subs, ', '));
    end
    addpath(cfgdir, '-begin');                            % 該 config 置最前、最優先
    clear mt_constants;                                  % R1：清同名快取，確保載 cfgdir 版本
    cfg = mt_constants();

    cfg.model = model;
    cfg.geom  = geom;
    cfg.N_I   = 6;
    cfg.select_ball = @select_ball;                      % 取點選球 helper（local；下游 cfg.select_ball(D,R)）

    % sanity：路由欄位齊備
    need = {'strategy','apdl_to_paper_idx','default_variant','regions'};
    for k = 1:numel(need)
        assert(isfield(cfg, need{k}), ...
            'config/%s/%s mt_constants.m 缺欄位 ''%s''', model, geom, need{k});
    end
end

% ---- 取點選球（併入 model_config；對外經 cfg.select_ball handle）----
function [P, Bstack, npts] = select_ball(D, R)
%SELECT_BALL  取 D 中 |p|<R 的節點，堆疊各 coil 場供擬合（model-agnostic）。
%   D = extract_ansys_data 回傳（D.r2=|p|² 旋轉不變、D.Ba N×3×N_I）；R = 球半徑 [m]（about WP）。
%   回 P(Np×3 actuator[m])、Bstack(3Np×N_I，欄 j = coil j 場堆成 [Bx;By;Bz;…])、npts。
    idx  = find(D.r2 < R^2);
    npts = numel(idx);
    P    = D.Pa(idx, :);
    N_I  = size(D.Ba, 3);
    Bstack = zeros(3*npts, N_I);
    for j = 1:N_I
        Bstack(:,j) = reshape(D.Ba(idx,:,j).', [], 1);
    end
end
