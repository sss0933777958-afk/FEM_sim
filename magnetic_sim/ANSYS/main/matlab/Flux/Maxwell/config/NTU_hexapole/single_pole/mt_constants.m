function c = mt_constants()
%MT_CONSTANTS  NTU single pole (Maxwell solver) per-model settings.
%   Geometry measured with OCC from the CAD source of truth:
%     CAD_model/NTU_hexapole/STEP/single_pole/single_pole.STEP   (pure analytic faces, mm)
%   Field: D:\Maxwell_sim\NTU\project\single_pole\single_pole.aedt -> export\single_pole\P1.fld
%
%   Structure follows config/long2016_hexapole_halfcut/tip40um/mt_constants.m, but this is a
%   SINGLE pole (one plate, one coil, one excitation), not a hexapole:
%     - no magic angle, no R_act rotation (pole axis is already +x), no 3x6 Pc_base
%     - origin = pole tip apex; l_hat = distance along x from the tip (user decision 2026-10-04)
%     - one excitation -> N_I = 1, apdl_to_paper_idx = 1
%   Do not copy values from the NTU hexapole (APDL) config or from long2016.

    % ===== (1) Routing =====
    c.strategy          = 'single_pole';             % translate to origin only (no rotation, no iron filter)
    c.apdl_to_paper_idx = 1;                         % one coil = P1 (identity, new model)
    c.default_variant   = 'maxwell';
    c.regions           = {'all'};
    c.R_load            = [];                        % no pre-filter ball; .fld box is already local
    c.N_I               = 1;                         % one excitation

    c.fld_dir           = 'D:\Maxwell_sim\NTU\export\single_pole';
    c.fld_files         = {'P1.fld'};
    %   Header: Min [-0.6 -0.6 -0.475] mm, Max [1.7 0.6 0.725] mm, step 0.02 mm (Maxwell global frame)
    %   -> 116 x 61 x 61 nodes. y = 0 and z = 0.125 (plate mid-plane) lie on grid planes;
    %   the tip x = 1.1306 mm does not (x grid = -0.6 + 0.02k).
    c.fld_files_variant.maxwell = c.fld_files;
    c.fld_variant_subdir = false;
    c.v_method          = 'grid';                    % regular .fld grid -> trilinear

    % ===== (2) Geometry (CAD measured, Maxwell global frame, [m]) =====
    % --- pole plate (steel, flat) ---
    c.PLATE_Z       = [0, 0.25e-3];                  % plate z range
    c.PLATE_T       = 0.25e-3;                       % thickness
    c.PLATE_ZMID    = 0.125e-3;                      % mid-plane = pole axis height
    c.PLATE_W       = 10e-3;                         % width (y = +-5 mm)
    c.PLATE_X       = [1.1306375e-3, 26.2776375e-3]; % x range (tip apex .. far end)
    c.POLE_TIP_R    = 5e-6;                          % tip fillet radius (centre x = 1.1356375 mm)
    c.pole_tip      = [1.1306375e-3; 0; c.PLATE_ZMID];   % tip apex on the mid-plane
    %   Outline in the x-y plane (symmetric about y = 0):
    %     tip wedge  half angle 7.1155 deg, x 1.135 .. 5.287 mm
    %     R10 blend  centres (4.0483461, +-10.4462517) mm
    %     body wedge half angle 16.8842 deg, x 6.953 .. 19.825 mm
    %     R5 end     around the post axis (21.2776375, 0) mm
    c.TIP_HALF_ANG  = 7.1155;                        % [deg]
    c.TIP_WEDGE_XMAX = 5.287e-3;                     % tip wedge ends here; iron filter is valid up to this x
    c.BODY_HALF_ANG = 16.8842;                       % [deg]

    % --- post (steel, part of the post+yoke solid) ---
    c.POST_XY       = [21.2776375e-3; 0];            % post axis = coil axis
    c.POST_R        = 2.5e-3;
    c.POST_Z        = [0.25e-3, 11.25e-3];

    % --- bushing (R3.5/R2.5 tube + R7.5 x 1 mm flange at z 6.05..7.05) ---
    c.BUSH_R        = [2.5e-3, 3.5e-3];              % [inner, outer]
    c.BUSH_Z        = [6.05e-3, 11.25e-3];

    % --- yoke bar on top of the post ---
    c.YOKE_Z        = [11.26e-3, 14.26e-3];
    c.YOKE_W        = 7e-3;                          % y = +-3.5 mm
    c.YOKE_X_END    = 36.2776375e-3;                 % R3.5 rounded end with R2.49 hole

    % --- coil (built in AEDT, not in the STEP) ---
    %   coil = Cyl(R7.5) - Cyl(R3.55), z 7.10..10.05 mm, axis = post axis (50 um off the bushing)
    %   Current1 = 50 A stranded, point out of terminal, section coil_Section1 (ZX, -x half)
    c.COIL_R        = [3.55e-3, 7.5e-3];
    c.COIL_Z        = [7.10e-3, 10.05e-3];

    % ===== (3) Calibration frame (user decision 2026-10-04) =====
    %   origin = pole tip apex (on the mid-plane); pole axis = +x (tip -> base).
    %   l_hat = distance along the x axis from the tip to the charge, i.e. charge at
    %   l_hat*(1, e_y, e_z) in this frame. Supersedes the 2026-10-02 frame (origin 500 um ahead of the tip).
    c.WP            = c.pole_tip;                    % origin in the Maxwell global frame [1.1306375 mm; 0; 0.125 mm]
    c.SPH_OFST      = 0;                             % no rotation; translation is done with c.WP
    c.pole_axis     = [1; 0; 0];                     % unit vector tip -> base
    c.pole_tip_wp   = [0; 0; 0];
    c.R_act         = eye(3);
    c.Pc_base       = [1; 0; 0];                     % charge direction for the single excitation
    c.pole_labels   = {'P1'};

    % Calibration region (user decision 2026-10-04): ball centred on the tip, air nodes only
    c.calib.center     = [0; 0; 0];                  % ball centre in the tip frame
    c.calib.R          = 150e-6;                     % ball radius
    c.calib.sampling   = 'axis';                     % Nr points x_k = -k*R/Nr (k = 1..Nr) on the air side, Bx/By/Bz used
    c.calib.fix_ex     = true;                       % e_x pinned to 0 (not identifiable with l_hat)
    c.calib.conv_tol   = 1e-4;                       % NMAE relative change < 0.01 % ...
    c.calib.conv_win   = 10;                         % ... for 10 consecutive levels

    % ===== (4) Sign convention + physical constants =====
    c.s_source = +1;                                 % P1.fld raw is source (user confirmed)
    c.N_c      = 50;                                 % turns; Maxwell 50 A = N_c * 1 A
    c.mu_0     = 4*pi*1e-7;
    c.k_m      = 1e-7;

    % ===== TODO =====
    %   (a) voltage path: no Hall sensor defined for the single pole (no S_hall / sensor geometry).
    %   (b) loaded by function/single_pole/sp_config.m (not hexapole/model_config.m, which forces N_I = 6).
end
