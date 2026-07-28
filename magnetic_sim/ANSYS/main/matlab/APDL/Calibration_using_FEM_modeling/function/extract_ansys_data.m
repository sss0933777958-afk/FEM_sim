function raw = extract_ansys_data(cfg, dataset, variant)
%EXTRACT_ANSYS_DATA  純資料提取：讀 6-coil FEM .dat → 原始 per-coil 場（measure frame）。
%   raw = EXTRACT_ANSYS_DATA(cfg, dataset, variant)
%   **只做「讀 .dat」**：迴圈 coil1..N_I 呼叫 import_ansys_data，回原始資料。
%   不濾鐵、不轉座標系、不翻號、不換單位、不建 F/R_act/Pc_base —— 那些是 to_actuator 的事。
%   cfg = model_config(model)（提供 model / N_I / default_variant / regions）。
%
%   raw：
%     .x .y .z   N×1      節點座標（measure frame [m]，取自 coil1；6 coil 共用同 mesh）
%     .B         N×3×N_I  各 coil 原始 B（measure frame，[T]；欄 = [bx by bz]）
%     .node_id   N×1      節點 id（coil1）
%     .model .variant .dataset  記錄用
%   需 ansys_path + import_ansys_data 在 path。
    if nargin < 2 || isempty(dataset), dataset = 'all'; end
    if nargin < 3 || isempty(variant), variant = cfg.default_variant; end
    model = cfg.model;
    N_I   = cfg.N_I;

    d1 = import_ansys_data(ansys_path(model,'data',variant,'coil1'), dataset, 'coil1');
    N  = numel(d1.x);
    B  = zeros(N, 3, N_I);
    B(:,:,1) = [d1.bx, d1.by, d1.bz];
    for k = 2:N_I
        cn = sprintf('coil%d', k);
        dk = import_ansys_data(ansys_path(model,'data',variant,cn), dataset, cn);
        B(:,:,k) = [dk.bx, dk.by, dk.bz];
    end

    raw = struct('x',d1.x, 'y',d1.y, 'z',d1.z, 'B',B, 'node_id',d1.node_id, ...
                 'model',model, 'variant',variant, 'dataset',dataset);
end
