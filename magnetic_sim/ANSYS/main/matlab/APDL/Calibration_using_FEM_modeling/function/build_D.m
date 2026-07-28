function D = build_D(raw, cfg)
%BUILD_D  raw（measure frame, Tesla）→ D（actuator frame, mT, canonical all-source）。
%   D = BUILD_D(raw, cfg)
%   依 cfg.strategy 分支；B 換單位 ×1e3(T→mT) 並套 canonical all-source
%   （逐激發欄 cfg.s_source：只翻「raw 是 sink」的極 → Ba 直接 all-source，下游不再補翻）。
%   回 .Pa(N×3[m]) .r2(N×1) .Ba(N×3×N_I[mT]) .R_act(3×3) .Pc_base(3×6) .F(6×N_I) [.WP]。
%   需 filter_iron_nodes（hex）在 path。
    N_I = cfg.N_I;
    F = zeros(6, N_I);
    for j = 1:N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end
    % 每激發欄的 all-source 號誌：coil j 激發 paper 極 p=apdl_to_paper_idx(j)，
    %   該極 raw sink → cfg.s_source(p)=−1（翻）；source → +1（不翻）。
    col_sign = cfg.s_source(cfg.apdl_to_paper_idx);          % 1×N_I

    switch cfg.strategy
        % ---------- long2016 / hung：magic-angle 六極 ----------
        case 'hex_magic'
            % 幾何（R_act / Pc_base）由 config「提供」——不再假設 baseline 魔術角（tip40um=canonical、
            % tip400um=尖端到尖端）。只驗「R_act 是正規正交旋轉」，不再驗 canonical。
            R_act   = cfg.R_act;
            Pc_base = cfg.Pc_base;
            assert(abs(det(R_act)-1) < 1e-6, 'R_act must be a proper rotation');
            assert(max(abs(R_act*R_act.' - eye(3)), [], 'all') < 1e-6, 'R_act must be orthonormal');

            air = filter_iron_nodes(raw.x, raw.y, raw.z, cfg, struct('visualize',false));
            zwp = raw.z - cfg.SPH_OFST;
            P_meas = [raw.x(air), raw.y(air), zwp(air)];
            Pa = (R_act * P_meas.').';
            r2 = sum(P_meas.^2, 2);

            Ba = zeros(sum(air), 3, N_I);
            for k = 1:N_I
                Bk = (1e3 * col_sign(k)) * raw.B(air,:,k);   % T→mT + all-source 翻 sink
                Ba(:,:,k) = (R_act * Bk.').';
            end
            D = struct('Pa',Pa, 'r2',r2, 'Ba',Ba, 'R_act',R_act, 'Pc_base',Pc_base, 'F',F);

        % ---------- NTU：扁平板（不旋轉不濾鐵） ----------
        case 'ntu_flat'
            WP = cfg.WP(:).';
            tipwp = [cfg.pole_tip_x_wp; cfg.pole_tip_y_wp; cfg.pole_tip_z_wp];
            Pc_base = tipwp ./ vecnorm(tipwp);
            R_load = cfg.R_load; if isempty(R_load), R_load = 2e-3; end

            P1  = [raw.x, raw.y, raw.z] - WP;
            sel = sum(P1.^2, 2) < R_load^2;
            P_meas = P1(sel, :);
            r2  = sum(P_meas.^2, 2);

            Ba = zeros(sum(sel), 3, N_I);
            for k = 1:N_I
                Ba(:,:,k) = (1e3 * col_sign(k)) * raw.B(sel,:,k);
            end
            D = struct('Pa',P_meas, 'r2',r2, 'Ba',Ba, 'R_act',eye(3), 'Pc_base',Pc_base, 'F',F, 'WP',WP);

        otherwise
            error('to_actuator:unknownStrategy', 'unknown strategy ''%s''', cfg.strategy);
    end
end
