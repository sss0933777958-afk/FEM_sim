function raw = extract_maxwell_data(cfg, dataset, variant)
%EXTRACT_MAXWELL_DATA  純資料提取：讀 N_I 個 Maxwell .fld → 原始 per-coil 場（Maxwell frame）。
%   raw = EXTRACT_MAXWELL_DATA(cfg, dataset, variant)
%   **只做「讀 .fld」**：迴圈 coil1..N_I 呼叫 import_maxwell_fld、組 raw。
%   不濾鐵、不轉座標系、不翻號、不換單位 —— 那些是 build_actuator_data 的事。
%   與 extract_ansys_data **同輸出介面** → 下游 build_actuator_data / main 直接吃（drop-in）。
%
%   .fld 檔位置由 cfg 提供（Maxwell config）：
%     cfg.fld_dir            .fld 所在資料夾（Maxwell 匯出根）
%     cfg.fld_files          （可選）1×N_I cell，各 coil 的 .fld 檔名；
%                            未給則預設 sprintf('coil%d.fld', k)
%     cfg.fld_variant_subdir （可選 true）→ 路徑再進 variant 子夾
%   ⚠ N_I 個 .fld 必須是「同一組匯出格點」（x/y/z 逐點相同）才能疊成 raw.B —— 內建檢查、不符即報錯。
%   ⚠ dataset 對 Maxwell 只是「標記」（一個 .fld = 一個場，無 all/wp/circuit 之分）。
%
%   raw（同 extract_ansys_data）：
%     .x .y .z   N×1      格點座標（Maxwell frame [m]，取自 coil1；N_I coil 共用同格）
%     .B         N×3×N_I  各 coil 原始 B（[T]；欄 = [bx by bz]）
%     .node_id   N×1      格點流水號（coil1）
%     .model .variant .dataset  記錄用
%   需 import_maxwell_fld 在 path。
    if nargin < 2 || isempty(dataset), dataset = 'all'; end
    if nargin < 3 || isempty(variant), variant = getdef(cfg,'default_variant',''); end
    model = cfg.model;
    N_I   = cfg.N_I;
    fdir  = fld_dir_for(cfg, variant);

    d1 = import_maxwell_fld(fullfile(fdir, coil_fname(cfg,1)));
    N  = numel(d1.x);
    B  = zeros(N, 3, N_I);
    B(:,:,1) = [d1.bx, d1.by, d1.bz];
    tol = 1e-9;
    for k = 2:N_I
        dk = import_maxwell_fld(fullfile(fdir, coil_fname(cfg,k)));
        assert(numel(dk.x)==N, ...
               'extract_maxwell_data:gridN', 'coil%d 格點數 %d ≠ coil1 %d（不同匯出格）', k, numel(dk.x), N);
        assert(max(abs(dk.x-d1.x))<tol && max(abs(dk.y-d1.y))<tol && max(abs(dk.z-d1.z))<tol, ...
               'extract_maxwell_data:gridXYZ', 'coil%d 格點座標與 coil1 不一致（.fld 須同一組匯出格點）', k);
        B(:,:,k) = [dk.bx, dk.by, dk.bz];
    end

    raw = struct('x',d1.x, 'y',d1.y, 'z',d1.z, 'B',B, 'node_id',d1.node_id, ...
                 'model',model, 'variant',variant, 'dataset',dataset);
end

% ---- local：解析 .fld 資料夾（可選 variant 子夾）----
function fdir = fld_dir_for(cfg, variant)
    base = cfg.fld_dir;
    if isfield(cfg,'fld_variant_subdir') && ~isempty(cfg.fld_variant_subdir) ...
            && cfg.fld_variant_subdir && ~isempty(variant)
        fdir = fullfile(base, variant);
    else
        fdir = base;
    end
end

% ---- local：第 k 個 coil 的 .fld 檔名 ----
function nm = coil_fname(cfg, k)
    if isfield(cfg,'fld_files') && ~isempty(cfg.fld_files)
        nm = cfg.fld_files{k};
    else
        nm = sprintf('coil%d.fld', k);
    end
end

% ---- local：cfg 取欄位或預設 ----
function v = getdef(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
