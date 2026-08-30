function ad = build_actuator_data(raw, cfg, opts)
%BUILD_ACTUATOR_DATA  raw（measure frame, Tesla）→ ad（mT, canonical all-source）。
%   ad = BUILD_ACTUATOR_DATA(raw, cfg)          校正用（預設：actuator frame + 濾鐵）
%   ad = BUILD_ACTUATOR_DATA(raw, cfg, opts)    繪圖用（可留鐵、可不旋轉）
%
%   把純讀的 raw 場轉成資料包。做四件事，其中兩件由 opts 控制：
%     (1) 換單位 ×1e3(T→mT)                         ── 一律做
%     (2) 套 canonical all-source 號誌 cfg.s_source ── 一律做
%     (3) 濾鐵 filter_iron_nodes                    ── opts.filter_iron
%     (4) 平移到原點 + 旋轉到 actuator frame        ── opts.frame
%
%   [ADDED 2026-08-30] opts（可選，省略時行為與先前逐字相同）：
%     .frame        'actuator'（預設）平移 z-SPH_OFST 後再乘 R_act ── 校正用
%                   'wp'              只平移、不旋轉               ── 疊 CAD 幾何的場圖用
%                   'measure'         都不做，維持 .fld 原始座標
%     .filter_iron  true（預設）只留空氣格點；false 保留鐵節點     ── 畫鐵芯內磁路要 false
%
%   起因：場圖腳本要的是「含鐵、WP frame」的場，拿不到就只好自己手刻 ×1e3 與號誌，
%   於是單位與 all-source 慣例散落在各腳本裡（`plot_p1p2_poles_3d` 甚至寫死 s=1 並註明待驗）。
%   開這兩個旗標之後，轉換只有這一份定義，繪圖端一律從這裡拿。
%
%   依 cfg.strategy 分支。回 struct
%     .Pa(N×3[m]) .r2(N×1) .Ba(N×3×M[mT]) .R_act(3×3) .Pc_base(3×6) .F(6×M) [.WP]
%     .frame .coils（記錄用）
%   ⚠ M = 讀進來的 coil 數（raw.coils 的長度），未給子集時 = N_I。
%   ⚠ .R_act 回的是**實際套用**的旋轉（frame≠actuator 時為 eye(3)），故恆有 Pa = R_act*P_meas；
%     模型自己的旋轉另存 .R_act_model。
%   ⚠ .r2 是「距該 frame 原點」的平方距離 —— frame='measure' 時原點是 CAD 原點、不是磁極共球心。
%   ⚠ 回傳的是「資料 struct」，不是 D̄ 轉移矩陣（那是 solve_voltage 的事）。
%   需 filter_iron_nodes（hex）在 path。
    if nargin < 3 || isempty(opts), opts = struct(); end
    frame   = lower(getdef(opts, 'frame', 'actuator'));
    do_iron = getdef(opts, 'filter_iron', true);
    assert(any(strcmp(frame, {'actuator','wp','measure'})), ...
           'build_actuator_data:frame', 'opts.frame 只能是 actuator | wp | measure（收到 ''%s''）', frame);

    N_I = cfg.N_I;
    % [ADDED 2026-08-30] 支援 coil 子集：raw.B 第 3 維依 raw.coils 的順序排。
    coils = getdef(raw, 'coils', 1:N_I);
    coils = coils(:).';
    assert(size(raw.B,3) == numel(coils), ...
           'build_actuator_data:coilsN', 'raw.B 第 3 維 %d ≠ raw.coils 長度 %d', size(raw.B,3), numel(coils));
    idx_paper = cfg.apdl_to_paper_idx(coils);                % 1×M：各欄激發的 paper 極
    M = numel(coils);
    F = zeros(6, M);
    for j = 1:M, F(idx_paper(j), j) = 1; end
    % 每激發欄的 all-source 號誌：該欄激發 paper 極 p，
    %   raw sink → cfg.s_source(p)=−1（翻）；source → +1（不翻）。
    col_sign = cfg.s_source(idx_paper);                      % 1×M

    switch cfg.strategy
        % ---------- long2016 / hung：magic-angle 六極 ----------
        case 'hex_magic'
            % 幾何（R_act / Pc_base）由 config「提供」——不再假設 baseline 魔術角（tip40um=canonical、
            % tip400um=尖端到尖端）。只驗「R_act 是正規正交旋轉」，不再驗 canonical。
            R_act   = cfg.R_act;
            Pc_base = cfg.Pc_base;
            assert(abs(det(R_act)-1) < 1e-6, 'R_act must be a proper rotation');
            assert(max(abs(R_act*R_act.' - eye(3)), [], 'all') < 1e-6, 'R_act must be orthonormal');

            % (3) 濾鐵：預設只留空氣格點；畫鐵芯內磁路時傳 filter_iron=false 全留
            if do_iron
                sel = filter_iron_nodes(raw.x, raw.y, raw.z, cfg, struct('visualize',false));
            else
                sel = true(numel(raw.x), 1);
            end
            % (4) frame：actuator = 平移 + 旋轉；wp = 只平移；measure = 都不動
            switch frame
                case 'actuator', zc = raw.z - cfg.SPH_OFST;  Rm = R_act;
                case 'wp',       zc = raw.z - cfg.SPH_OFST;  Rm = eye(3);
                case 'measure',  zc = raw.z;                 Rm = eye(3);
            end
            P_meas = [raw.x(sel), raw.y(sel), zc(sel)];
            Pa = (Rm * P_meas.').';
            r2 = sum(P_meas.^2, 2);

            Ba = zeros(sum(sel), 3, M);
            for k = 1:M
                Bk = (1e3 * col_sign(k)) * raw.B(sel,:,k);   % (1) T→mT + (2) all-source 翻 sink
                Ba(:,:,k) = (Rm * Bk.').';
            end
            ad = struct('Pa',Pa, 'r2',r2, 'Ba',Ba, 'R_act',Rm, 'R_act_model',R_act, ...
                        'Pc_base',Pc_base, 'F',F);

        % ---------- NTU：扁平板（不旋轉不濾鐵） ----------
        case 'ntu_flat'
            % 這條路徑本來就不濾鐵、不旋轉（只把座標平移到 WP），frame/filter_iron 無從套用。
            % 給了非預設值就明報錯，不要靜默忽略。
            assert(strcmp(frame,'actuator') && do_iron, ...
                   'build_actuator_data:ntuOpts', ...
                   'strategy=''ntu_flat'' 不支援 opts.frame / opts.filter_iron（它不濾鐵也不旋轉）');
            WP = cfg.WP(:).';
            tipwp = [cfg.pole_tip_x_wp; cfg.pole_tip_y_wp; cfg.pole_tip_z_wp];
            Pc_base = tipwp ./ vecnorm(tipwp);
            R_load = cfg.R_load; if isempty(R_load), R_load = 2e-3; end

            P1  = [raw.x, raw.y, raw.z] - WP;
            sel = sum(P1.^2, 2) < R_load^2;
            P_meas = P1(sel, :);
            r2  = sum(P_meas.^2, 2);

            Ba = zeros(sum(sel), 3, M);
            for k = 1:M
                Ba(:,:,k) = (1e3 * col_sign(k)) * raw.B(sel,:,k);
            end
            ad = struct('Pa',P_meas, 'r2',r2, 'Ba',Ba, 'R_act',eye(3), 'R_act_model',eye(3), ...
                        'Pc_base',Pc_base, 'F',F, 'WP',WP);

        otherwise
            error('build_actuator_data:unknownStrategy', 'unknown strategy ''%s''', cfg.strategy);
    end
    ad.frame = frame;                 % 記錄用：下游/除錯可確認拿到的是哪個座標系
    ad.coils = coils;
end

% ---- local：opts/raw 取欄位或預設 ----
function v = getdef(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
