function [air_mask, debug_info] = filter_iron_nodes(x_apdl, y_apdl, z_apdl, c, opts)
% FILTER_IRON_NODES  Geometric cone model to exclude iron nodes from WP data
%   [air_mask, debug_info] = filter_iron_nodes(x, y, z, c)
%   [air_mask, debug_info] = filter_iron_nodes(x, y, z, c, opts)
%
%   Inputs:
%     x_apdl, y_apdl, z_apdl - node coordinates in APDL frame [m] (Nx1)
%     c       - constants struct from mt_constants()
%     opts    - optional struct with fields:
%       .visualize  (false)  generate 3D scatter plot
%       .safety_r   (0)      extra safety sphere radius around each tip [m]
%
%   [MODIFIED 2026-08-14] safety_r default 100e-6 -> 0 (user decision).
%     The tips sit 500 um from the origin, so a 100 um safety sphere reached
%     inward to r = 400 um and punched six holes out of every sampling ball
%     with R > 400 um (measured: 0 nodes removed at R = 400, 70 at R = 420,
%     1452 = 2.23% at R = 500). Those nodes are AIR, not iron -- the cone
%     envelope below already covers the iron, and nodes in front of the tip
%     (s < 0) are outside the pole by construction. Keeping the sphere made
%     "a ball of radius R" untrue for R > 400 um. Set safety_r > 0 explicitly
%     if a particular study needs to stay clear of the tip singularity.
%
%   Outputs:
%     air_mask   - logical Nx1, true = air node (keep), false = iron (exclude)
%     debug_info - struct with per-pole exclusion counts

    if nargin < 5, opts = struct(); end
    if ~isfield(opts, 'visualize'), opts.visualize = false; end
    if ~isfield(opts, 'safety_r'),  opts.safety_r  = 0; end

    N = numel(x_apdl);

    % Convert to WP-centered coordinates
    x_wp = x_apdl;
    y_wp = y_apdl;
    z_wp = z_apdl - c.SPH_OFST;

    % Check each of 6 poles
    iron_mask = false(N, 1);
    pole_counts = zeros(1, 6);

    % [ADDED 2026-08-15] 錐體濾鐵需要 POLE_R / POLE_CONE_LEN / pole_axis。**平板型磁極**
    %   （zhi_peng：厚 686 µm 的極板，非圓錐）沒有這三個欄位 → 跳過錐體判定，只留 safety
    %   sphere，並明確警告。取樣半徑遠小於極尖距離時（R << R_norm）本來就取不到鐵，安全；
    %   要在近極尖區用本函式，必須先為該 model 定義板狀鐵件包絡（見其 config 的 TODO）。
    has_cone = isfield(c,'POLE_R') && isfield(c,'POLE_CONE_LEN') && isfield(c,'pole_axis');
    % [ADDED 2026-08-27 使用者要求] **平板型磁極的鐵件包絡**（zhi_peng）。原本缺錐體欄位時
    %   完全不濾（in_cone = false），對 R > IRON_SAFE_R 的取樣球會把鐵當空氣。
    %   包絡分兩段（config 的實測值，見 mt_constants「磁極板本體（平板；濾鐵包絡用）」）：
    %     ① 舌片 r_xy ∈ [R_norm_xy, POLE_WEDGE_R(2)]：厚 POLE_TIP_BAND（貼板的**單側**，
    %        下極貼板底、上極貼板頂），平面內是**與尖端 40 µm 圓角相切**的楔形（半角
    %        POLE_WEDGE_HALF）；圓角本身以圓柱補上。
    %     ② 外段矩形塊 r_xy ∈ POLE_BLOCK_R：寬 POLE_BLOCK_W、厚**整個板厚** POLE_PLATE_T。
    %   楔形半寬（u 自圓角軸心沿極方位量）：w(u) = R_tip/cos(beta) + u*tan(beta)
    %     —— 切線到圓心的垂距 = R_tip，故 u=0 處半寬是 R_tip/cos(beta) 而非 R_tip。
    has_plate = isfield(c,'POLE_TIP_BAND') && isfield(c,'POLE_WEDGE_HALF') && ...
                isfield(c,'pole_is_lower') && isfield(c,'R_norm_xy') && isfield(c,'R_norm_z');
    if ~has_cone && ~has_plate
        warning('filter_iron_nodes:noPoleGeom', ...
            ['cfg 既無錐體幾何（POLE_R/POLE_CONE_LEN）也無平板幾何（POLE_TIP_BAND/' ...
             'POLE_WEDGE_HALF）→ **未做濾鐵**，僅套 safety sphere = %g um。'], ...
            opts.safety_r*1e6);
    end
    if has_plate
        bet   = c.POLE_WEDGE_HALF * pi/180;
        r_ax  = c.R_norm_xy + c.POLE_TIP_R;          % 尖端圓角軸心的徑向位置
        wR    = c.POLE_WEDGE_R;                      % 舌片徑向範圍 [內, 外]
        hasBlk = isfield(c,'POLE_BLOCK_R') && isfield(c,'POLE_BLOCK_W') && isfield(c,'POLE_PLATE_T');
    end

    for i = 1:6
        % Tip position in WP frame
        tip = [c.pole_tip_x(i); c.pole_tip_y(i); c.pole_tip_z_wp(i)];

        % Vector from tip to each node
        vx = x_wp - tip(1);
        vy = y_wp - tip(2);
        vz = z_wp - tip(3);

        % Distance from tip
        dist = sqrt(vx.^2 + vy.^2 + vz.^2);

        % Safety sphere: exclude everything within safety_r of tip
        in_safety = dist < opts.safety_r;

        if has_cone
            % Project onto pole axis (positive = behind tip, toward base)
            ax = c.pole_axis(:, i);
            s = vx*ax(1) + vy*ax(2) + vz*ax(3);

            % Perpendicular distance from axis
            r_perp = sqrt(dist.^2 - s.^2);

            % Cone radius at distance s along axis
            r_cone = c.POLE_TIP_R + s * (c.POLE_R - c.POLE_TIP_R) / c.POLE_CONE_LEN;

            % Inside cone: positive projection AND within cone envelope
            in_cone = (s > 0) & (r_perp < r_cone) & (s < c.POLE_CONE_LEN);
        elseif has_plate
            % --- 平板包絡（見上方註解）---
            th  = atan2(c.pole_tip_y(i), c.pole_tip_x(i));
            u   = (x_wp - r_ax*cos(th))*cos(th) + (y_wp - r_ax*sin(th))*sin(th);   % 沿極方位
            v   = -(x_wp - r_ax*cos(th))*sin(th) + (y_wp - r_ax*sin(th))*cos(th);  % 橫向
            rxy = hypot(x_wp, y_wp);
            % 舌片的單側 z 帶：下極貼板底、上極貼板頂
            if c.pole_is_lower(i)
                inz = (z_wp >= -c.R_norm_z) & (z_wp <= -c.R_norm_z + c.POLE_TIP_BAND);
            else
                inz = (z_wp <=  c.R_norm_z) & (z_wp >=  c.R_norm_z - c.POLE_TIP_BAND);
            end
            w_half  = c.POLE_TIP_R/cos(bet) + max(u,0)*tan(bet);
            in_wedge = inz & (rxy <= wR(2)) & ...
                       ( ((u >= 0) & (abs(v) <= w_half)) | ...      % 楔形本體
                         ((u <  0) & (hypot(u,v) <= c.POLE_TIP_R)) );  % 尖端圓角
            if hasBlk
                inz_b = abs(z_wp) <= c.POLE_PLATE_T/2;              % 外段是整個板厚
                in_blk = inz_b & (rxy > c.POLE_BLOCK_R(1)) & (rxy <= c.POLE_BLOCK_R(2)) & ...
                         (abs(v) <= c.POLE_BLOCK_W/2) & (u > 0);
            else
                in_blk = false(N,1);
            end
            in_cone = in_wedge | in_blk;
        else
            in_cone = false(N, 1);          % 既無錐體也無平板幾何（見上方警告）
        end

        pole_mask = in_cone | in_safety;
        iron_mask = iron_mask | pole_mask;
        pole_counts(i) = sum(pole_mask);
    end

    air_mask = ~iron_mask;

    debug_info.total_nodes = N;
    debug_info.iron_count  = sum(iron_mask);
    debug_info.air_count   = sum(air_mask);
    debug_info.pole_counts = pole_counts;
    debug_info.pole_labels = c.pole_labels;

    fprintf('Iron exclusion: %d / %d nodes removed (%.2f%%)\n', ...
        debug_info.iron_count, N, 100*debug_info.iron_count/N);
    for i = 1:6
        fprintf('  %s: %d nodes\n', c.pole_labels{i}, pole_counts(i));
    end

    % Optional visualization
    if opts.visualize
        figure('Name', 'Iron Exclusion', 'Position', [100 100 800 600]);
        scatter3(x_wp(air_mask)*1e3, y_wp(air_mask)*1e3, z_wp(air_mask)*1e3, ...
            1, 'b', '.', 'DisplayName', 'Air');
        hold on;
        scatter3(x_wp(iron_mask)*1e3, y_wp(iron_mask)*1e3, z_wp(iron_mask)*1e3, ...
            8, 'r', 'filled', 'DisplayName', 'Iron (excluded)');
        scatter3(c.pole_tip_x*1e3, c.pole_tip_y*1e3, c.pole_tip_z_wp*1e3, ...
            80, 'g', 'filled', 'DisplayName', 'Pole tips');
        for i = 1:6
            text(c.pole_tip_x(i)*1e3, c.pole_tip_y(i)*1e3, c.pole_tip_z_wp(i)*1e3, ...
                ['  ' c.pole_labels{i}], 'FontSize', 10, 'FontWeight', 'bold');
        end
        hold off;
        xlabel('x_{wp} [mm]'); ylabel('y_{wp} [mm]'); zlabel('z_{wp} [mm]');
        title('Iron Node Exclusion (Geometric Cone Model)');
        legend('Location', 'best');
        axis equal; grid on;
        view(30, 25);
    end
end
