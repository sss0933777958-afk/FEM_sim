function data = import_maxwell_fld(fld_path)
%IMPORT_MAXWELL_FLD  讀一個 Maxwell「Export on a grid」.fld → 單一激發(coil)場。
%   data = IMPORT_MAXWELL_FLD(fld_path)
%   .fld 格式（ANSYS Electronics Desktop / Maxwell）：
%     檔頭數行（"Grid Output Min:..." / "X, Y, Z, Vector data ..."）—— 自動偵測略過，
%     之後每行 "x y z Bx By Bz"（座標 mm、B 特斯拉）。
%
%   回傳 struct（與 import_ansys_data 同介面，供下游 drop-in）：
%     node_id  N×1   格點流水號（合成，1..N；Maxwell 無節點 id）
%     x,y,z    N×1   座標 [m]（mm→m）
%     bx,by,bz N×1   B 分量 [T]
%     bsum     N×1   |B| [T]
%
%   ⚠ .fld 座標系 = Maxwell/STEP frame（= APDL raw frame，WP 在 z≈-12.71mm、SPH_OFST），
%     與 import_ansys_data 一致 → 下游 build_actuator_data 的 zoff=-SPH_OFST 照樣成立。
%   ⚠ 大檔（全域 ±60mm/0.5mm ≈ 2GB）讀起來慢又吃記憶體 —— 匯出時只給 WP/sensor/磁鐵區的
%     目標小區細格（見 memory project_long2016_maxwell_validation）。

    if nargin < 1 || isempty(fld_path)
        error('import_maxwell_fld:noPath', '需要 .fld 路徑');
    end
    fid = fopen(fld_path, 'r');
    if fid < 0, error('import_maxwell_fld:open', '開不了 %s', fld_path); end

    % --- 偵測檔頭行數（peek 前幾行，第一行能解出 ≥6 個數字即為資料）---
    nh = 0;
    for i = 1:8
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        if numel(sscanf(ln, '%f')) >= 6, break; end
        nh = nh + 1;
    end
    frewind(fid);

    % --- 串流讀 6 欄（x y z Bx By Bz）---
    C = textscan(fid, '%f%f%f%f%f%f', 'HeaderLines', nh, ...
                 'CollectOutput', true, 'MultipleDelimsAsOne', true);
    fclose(fid);
    M = C{1};
    if isempty(M)
        error('import_maxwell_fld:empty', '%s 沒解析到資料（檔頭偵測 nh=%d，請檢查格式）', fld_path, nh);
    end

    data.node_id = (1:size(M,1)).';
    data.x  = M(:,1) * 1e-3;   % mm → m
    data.y  = M(:,2) * 1e-3;
    data.z  = M(:,3) * 1e-3;
    data.bx = M(:,4);          % T
    data.by = M(:,5);
    data.bz = M(:,6);
    data.bsum = sqrt(M(:,4).^2 + M(:,5).^2 + M(:,6).^2);
    fprintf('  Maxwell .fld: %d 格點, |B|max=%.4f T  (%s)\n', ...
            size(M,1), max(data.bsum), fld_path);
end
