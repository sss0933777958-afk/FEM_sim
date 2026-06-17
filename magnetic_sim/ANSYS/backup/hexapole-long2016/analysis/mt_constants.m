function c = mt_constants()
% MT_CONSTANTS  Shared geometry constants for magnetic tweezers analysis
%   c = mt_constants() returns struct with all geometry parameters,
%   pole naming conventions, and coordinate transforms.
%   Must match values in APDL scripts.

    % Working radius
    c.R_norm    = 500e-6;                              % [m]
    c.R_norm_xy = c.R_norm * sqrt(2/3);                % ~408 um
    c.R_norm_z  = c.R_norm / sqrt(3);                  % ~289 um

    % Yoke and protrusion geometry
    c.PROT_H    = 7.0e-3;                              % protrusion height [m]
    c.PROT_R    = 10.0e-3 / 2;                         % protrusion radius [m]
    c.YOKE_IN_R  = 84e-3 / 2;                          % yoke inner radius [m]
    c.YOKE_OUT_R = 106e-3 / 2;                         % yoke outer radius [m]
    c.YOKE_MID_R = (c.YOKE_IN_R + c.YOKE_OUT_R) / 2;  % ~47.5 mm

    % Working point offset from APDL origin (yoke base center)
    c.SPH_OFST = -c.PROT_H - 6e-3 + c.R_norm_z;       % ~ -12.711 mm

    % Pole naming: Paper convention
    %   Paper index:  P1    P2     P3     P4     P5    P6
    %   Angle (deg):  0     180    120    300    60    240
    %   Layer:        Lower Upper  Lower  Upper  Upper Lower
    c.pole_angles = [0, 180, 120, 300, 60, 240];        % degrees, indexed P1-P6
    c.pole_labels = {'P1','P2','P3','P4','P5','P6'};

    % APDL coil index -> Paper pole name
    %   Coil1->P1, Coil2->P3, Coil3->P6, Coil4->P5, Coil5->P2, Coil6->P4
    c.apdl_to_paper = {'P1','P3','P6','P5','P2','P4'};

    % Pole tip positions in WP-centered coordinates [m]
    %   x_wp = R_norm_xy * cos(angle), y_wp = R_norm_xy * sin(angle)
    %   z_wp = -R_norm_z for lower, +R_norm_z for upper
    c.pole_tip_x = c.R_norm_xy * cosd(c.pole_angles);   % [m]
    c.pole_tip_y = c.R_norm_xy * sind(c.pole_angles);   % [m]
    c.pole_tip_z_wp = [-1, +1, -1, +1, +1, -1] * c.R_norm_z;  % [m] P1-P6

    % Cone geometry (for iron exclusion filter)
    c.POLE_TIP_R    = 40e-6;       % tip radius [m]
    c.POLE_R        = 3e-3;        % base radius [m] (half of 6 mm diameter)
    c.POLE_CONE_LEN = 15e-3;       % cone length [m]

    % Upper pole incline angle (from APDL lines 218-233)
    %   tip -> (R_norm_xy, 0, SPH_OFST + R_norm_z) in APDL coords
    %   end -> (YOKE_MID_R - 11.5mm, 0, YOKE_H + PROT_H + 5mm)
    %   INCLINE_ANG = atan2(dz, dr_horizontal)
    YOKE_H = 2e-3;
    end_upper_r = c.YOKE_MID_R - 11.5e-3;            % 36 mm
    end_upper_z = YOKE_H + c.PROT_H + 5e-3;          % 14 mm (APDL coords)
    tip_upper_z = -c.PROT_H - 6e-3 + 2*c.R_norm_z;   % APDL coords
    dz = end_upper_z - tip_upper_z;
    dr = end_upper_r - c.R_norm_xy;
    c.upper_incline = atan2(dz, dr);  % ~36.6 deg from horizontal

    % Pole axis directions (unit vectors, tip -> base, in WP frame)
    %   Lower: horizontal radial outward
    %   Upper: inclined upward at upper_incline
    c.pole_is_lower = [1, 0, 1, 0, 0, 1];  % P1-P6
    c.pole_axis = zeros(3, 6);
    for i = 1:6
        theta = c.pole_angles(i) * pi/180;
        if c.pole_is_lower(i)
            c.pole_axis(:,i) = [cos(theta); sin(theta); 0];
        else
            inc = c.upper_incline;
            c.pole_axis(:,i) = [cos(inc)*cos(theta); cos(inc)*sin(theta); sin(inc)];
        end
    end

    % Physical constants for point-charge model
    c.N_c   = 70;                  % coil turns per pole
    c.mu_0  = 4*pi*1e-7;          % vacuum permeability [H/m]
    c.k_m   = 1e-7;               % mu_0 / (4*pi) [N/A^2]
    c.alpha = atan2(c.R_norm_xy, c.R_norm_z);  % polar angle of tips on charge sphere (~54.74 deg)
end
