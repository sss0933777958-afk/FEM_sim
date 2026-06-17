function [Bx, By, Bz] = point_charge_model(p_wp, ell, R_a, I_vec, K_I, c)
% POINT_CHARGE_MODEL  Evaluate B-field from 6-pole point-charge model
%   [Bx, By, Bz] = point_charge_model(p_wp, ell, R_a, I_vec, K_I, c)
%
%   Inputs:
%     p_wp  - Nx3 array of field points in WP frame [m]
%     ell   - effective charge distance from WP center [m] (fitted parameter)
%     R_a   - lumped air reluctance [A/Wb] (fitted parameter)
%     I_vec - 6x1 current vector in paper pole order (P1-P6) [A]
%     K_I   - 6x6 flux distribution matrix (paper pole order)
%     c     - constants struct from mt_constants()
%
%   Outputs:
%     Bx, By, Bz - Nx1 B-field components in WP frame [T]
%
%   Model: Long 2016 dissertation Eq. 2.1-2.4
%     Charge positions: c_i = ell * (tip_direction_i)
%     Charge magnitudes: Q = (N_c / (mu_0 * R_a)) * K_I * I_vec
%     B(p) = sum_i k_m * Q(i) * (p - c_i) / |p - c_i|^3

    N = size(p_wp, 1);
    Bx = zeros(N, 1);
    By = zeros(N, 1);
    Bz = zeros(N, 1);

    % Place 6 charges on sphere of radius ell
    % Direction: same as pole tip direction from WP center
    alpha = c.alpha;  % polar angle (~54.74 deg)
    charge_pos = zeros(3, 6);
    for i = 1:6
        theta = c.pole_angles(i) * pi/180;
        z_sign = sign(c.pole_tip_z_wp(i));  % -1 for lower, +1 for upper
        charge_pos(:,i) = ell * [cos(theta)*sin(alpha); ...
                                  sin(theta)*sin(alpha); ...
                                  z_sign*cos(alpha)];
    end

    % Compute charge vector: Q = -(N_c / (mu_0 * R_a)) * K_I * I_vec
    % Sign convention: APDL coil winding drives flux INTO the excited pole
    % (excited pole = flux sink), so q_i has opposite sign from K_I * F_mmf.
    Q = -(c.N_c / (c.mu_0 * R_a)) * (K_I * I_vec);  % 6x1, units [A*m]

    % Accumulate B-field from each charge (vectorized over N nodes)
    for i = 1:6
        % Vector from charge i to each field point
        dx = p_wp(:,1) - charge_pos(1,i);
        dy = p_wp(:,2) - charge_pos(2,i);
        dz = p_wp(:,3) - charge_pos(3,i);

        r3 = (dx.^2 + dy.^2 + dz.^2).^(3/2);

        Bx = Bx + Q(i) * dx ./ r3;
        By = By + Q(i) * dy ./ r3;
        Bz = Bz + Q(i) * dz ./ r3;
    end

    % Scale by k_m
    Bx = c.k_m * Bx;
    By = c.k_m * By;
    Bz = c.k_m * Bz;
end
