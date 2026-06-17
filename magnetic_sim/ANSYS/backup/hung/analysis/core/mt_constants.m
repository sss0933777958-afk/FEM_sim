function c = mt_constants()
% MT_CONSTANTS  Hung hexapole design constants
%   WP center at origin (0,0,0) — no SPH_OFST offset

    c.R_norm    = 500e-6;                      % working radius [m]
    c.R_norm_xy = c.R_norm * sqrt(2/3);        % 408.2 um
    c.R_norm_z  = c.R_norm / sqrt(3);          % 288.7 um
    c.SPH_OFST  = 0;                           % Hung WP = origin

    % Pole angles
    c.TILT_UP = 35.0;       % degrees
    c.TILT_DN = 5.71;       % degrees
    c.magic_angle = 54.74;  % degrees

    % Pole geometry
    c.POLE_R         = 3.175e-3;
    c.POLE_TIP_R     = 40e-6;         % tip fillet radius [m]
    c.POLE_CONE_LEN  = 15.875e-3;
    c.POLE_TOTAL_LEN = 43.0e-3;

    % Pole naming (Paper convention P1-P6)
    c.pole_angles = [0, 180, 120, 300, 60, 240];   % azimuthal [deg]
    c.pole_labels = {'P1','P2','P3','P4','P5','P6'};
    c.pole_is_lower = [1, 0, 1, 0, 0, 1];          % P1-P6

    % APDL coil index -> Paper pole index (Hung: identity mapping)
    %   Coil1->P1, Coil2->P2, Coil3->P3, Coil4->P4, Coil5->P5, Coil6->P6
    c.apdl_to_paper = {'P1','P2','P3','P4','P5','P6'};
    c.apdl_to_paper_idx = [1, 2, 3, 4, 5, 6];

    % 6 pole tip positions (Hung: on sphere R=0.5mm at origin)
    fc_h = cosd(c.magic_angle) * c.R_norm;
    fc_r = sind(c.magic_angle) * c.R_norm;
    layer = [-1, 1, -1, 1, 1, -1];       % -1=lower, +1=upper
    c.tip_x = zeros(1,6);
    c.tip_y = zeros(1,6);
    c.tip_z = zeros(1,6);
    for i = 1:6
        c.tip_x(i) = cosd(c.pole_angles(i)) * fc_r;
        c.tip_y(i) = sind(c.pole_angles(i)) * fc_r;
        c.tip_z(i) = layer(i) * fc_h;
    end
    c.pole_tip_z_wp = c.tip_z;  % alias for fitting scripts

    % Physical constants for point-charge model
    c.N_c   = 70;                  % coil turns per pole
    c.mu_0  = 4*pi*1e-7;          % vacuum permeability [H/m]
    c.k_m   = 1e-7;               % mu_0 / (4*pi) [N/A^2]
    c.alpha = atan2(c.R_norm_xy, c.R_norm_z);  % ~54.74 deg [rad]

    % Yoke geometry (for wide-view plot)
    c.YOKE_RI = 38.0e-3;
    c.YOKE_RO = 62.5e-3;
end
