function c = mt_constants()
%MT_CONSTANTS  Zhi-peng pair pole (Maxwell solver) per-model settings.
%   Geometry source: CAD_model/zhi-Peng/STEP/zhi_peng_R500.STEP with only P1 and P2 kept (P3-P6
%   poles / posts / coils removed). Simulated in D:\Maxwell_sim\Zhi_peng\project\Pair_pole\Pair_pole.aedt,
%   which differs from the STEP in two places (2026-10-01 mesh fix):
%     - the P2 plate is thickened 0.5773 -> 0.5800 mm (mirror of P1; P2 tongue z 0.402..0.580)
%     - both coils rebuilt as primitives Cyl(R6.5) - Cyl(R1.55), centred on the true post axes
%   Values below follow the simulated (AEDT) geometry; STEP values are kept as *_cad.
%
%   Structure follows config/NTU_hexapole/single_pole/mt_constants.m:
%     - origin = P1 tip apex; pole axis +x (tip -> base, into P1); l_hat > 0 means a charge inside P1
%     - two poles, two excitations (P1.fld = coil 1 = P1, P2.fld = coil 2 = P2)
%   Not a magic-angle hexapole any more: no R_act rotation is used.
%
%   Differences from NTU (do not copy values across):
%     - the tip is a VERTICAL fillet edge (40 um radius, 178 um tall tongue), not a point;
%       the origin is put on the apex line at mid-tongue height
%     - P1 is a lower pole (tongue at the plate bottom), P2 an upper pole (tongue at the plate top):
%       the two tips are NOT coaxial -- P2 tip sits 0.8163 mm toward -x and 0.402 mm higher
%     - 70 turns, 1 A per turn in the AEDT winding (NTU: 50 A stranded = 50 turns x 1 A)

    % ===== (1) Routing =====
    c.strategy          = 'pair_pole';
    c.apdl_to_paper_idx = [1, 2];                    % identity (new model): P1.fld = P1, P2.fld = P2
    c.default_variant   = 'maxwell';
    c.regions           = {'all'};
    c.R_load            = [];
    c.N_I               = 2;                         % two excitations

    c.fld_dir           = 'D:\Maxwell_sim\Zhi_peng\export\Pair_pole';
    c.fld_files         = {'P1.fld', 'P2.fld'};
    %   Header (both): Min [-2 -2 -1.711] mm, Max [2 2 2.289] mm, step 0.02 mm (Maxwell global = STEP frame,
    %   z = 0 at the plate bottom) -> 201^3 nodes. 2026-10-01 export (second version; the first one was
    %   centred at z = -12.71 mm and is void).
    %   Checked 2026-10-05: 20 um in front of each tip at mid-tongue height,
    %     P1.fld: |B| 327.07 mT at the +x tip vs 161.26 mT at the -x tip, B_x < 0 (away from P1)
    %     P2.fld: |B| 324.53 mT at the -x tip vs 163.44 mT at the +x tip, B_x > 0 (away from P2)
    %   -> P1.fld excites P1, P2.fld excites P2, and both raw fields are source.
    c.fld_files_variant.maxwell = c.fld_files;
    c.fld_variant_subdir = false;
    c.v_method          = 'grid';

    % ===== (2) Geometry (Maxwell global frame, [m]) =====
    % --- tips (OCC on zhi_peng_R500.STEP, 40 um fillet = vertical cylinder) ---
    c.POLE_TIP_R    = 40e-6;                         % fillet radius
    c.TONGUE_T      = 178e-6;                        % tongue thickness (tip edge height)
    c.P1_TIP_X      = 0.40800e-3;                    % apex x (fillet centre 0.44800 mm)
    c.P1_TIP_Z      = [0, 0.178e-3];                 % P1 tongue at the plate bottom
    c.P2_TIP_X      = -0.40830e-3;                   % apex x (fillet centre -0.44830 mm)
    c.P2_TIP_Z      = [0.402e-3, 0.580e-3];          % AEDT (plate thickened to 0.58 mm)
    c.P2_TIP_Z_cad  = [0.3993e-3, 0.5773e-3];        % STEP
    c.TIP_GAP_X     = c.P1_TIP_X - c.P2_TIP_X;       % 816.3 um horizontal tip-to-tip
    c.pole_tip_P1   = [c.P1_TIP_X; 0; mean(c.P1_TIP_Z)];   % [0.408; 0; 0.089] mm
    c.pole_tip_P2   = [c.P2_TIP_X; 0; mean(c.P2_TIP_Z)];   % [-0.4083; 0; 0.491] mm

    % --- tongue / plate (from config/zhi_peng/R500, unchanged by the pair cut) ---
    c.TIP_HALF_ANG  = 7.5;                           % tongue wedge half angle [deg]
    c.TONGUE_R      = [0.408e-3, 15.0e-3];           % tongue radial range (tip -> rectangular block)
    c.BLOCK_R       = [15.0e-3, 25.0e-3];            % outer block radial range
    c.BLOCK_W       = 6.0e-3;                        % outer block width (+-3 mm)
    c.PLATE_T       = 0.580e-3;                      % plate thickness, both poles in the AEDT model

    % --- posts and coils (AEDT) ---
    c.POST_XY       = [21.5015694e-3, -21.5015694e-3; 0, 0];   % [P1, P2] post axes = coil axes
    c.POST_R        = 1.5e-3;
    c.COIL_R        = [1.55e-3, 6.5e-3];             % coil = Cyl(R6.5) - Cyl(R1.55)
    c.COIL_Z        = [0.58e-3, 5.18e-3];
    %   Winding1/2: 70 conductors, Current = if(coil_idx==k, 1A, 0A);
    %   terminals Coil1 point-out true, Coil2 point-out false.

    % ===== (3) Calibration frame =====
    %   origin = P1 tip apex at mid-tongue height; pole axis = +x (tip -> base, into P1).
    %   l_hat is the charge's x coordinate in this frame: l_hat > 0 inside P1.
    c.WP            = c.pole_tip_P1;                 % origin in the Maxwell global frame
    c.SPH_OFST      = 0;                             % no rotation; translation is done with c.WP
    c.pole_axis     = [1; 0; 0];
    c.pole_tip_wp   = [zeros(3,1), c.pole_tip_P2 - c.pole_tip_P1];   % [P1, P2] = [0; -816.3 um, 0, +402 um]
    c.R_act         = eye(3);
    c.pole_labels   = {'P1', 'P2'};

    % ===== (4) Sign convention + physical constants =====
    c.s_source = [+1, +1];                           % both raw fields are source (checked above)
    c.N_c      = 70;                                 % turns
    c.I_exc    = 1;                                  % AEDT winding current [A] (= 70 A-turn)
    c.mu_0     = 4*pi*1e-7;
    c.k_m      = 1e-7;

    % ===== TODO =====
    %   (a) calibration sampling (number and positions of the points) is set in the analysis scripts.
    %   (b) no loader yet: function/single_pole/sp_config.m requires strategy = 'single_pole'.
    %   (c) no iron mask: the tongue tip is a vertical edge; write one before sampling near the tips.
end
