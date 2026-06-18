%% gen_Vout_Vin_1p572.m
%  ------------------------------------------------------------------
%  Compute V_out/V_in 6x6 matrix using sensor placement at 1.572 mm
%  (sensor variant; baseline minus 3.0 mm).
%
%  Chain (per Bs_derivation.pdf):
%      B_bar    = (S_hall * k_A / I_in) * B_FEM       [dimensionless]
%      B_S      = (S_hall / I_in)        * B_FEM       [V/A]
%      V_out/V_in = k_A * B_S = (S_hall * k_A / I_in) * B_FEM
%
%  Sensor placement (changed: 4.572 -> 1.572 mm; baseline - 3.0 mm):
%      LOWER (P1/P3/P6): pole-local (1.572, 0.41, 0) mm relative to tip
%      UPPER (P2/P4/P5): pole-local decomposition of
%                          1.572 along cone slant + 0.41 outward normal
%                          axial      = 1.572*cos(beta) - 0.41*sin(beta) = 1.4609
%                          transverse = 1.572*sin(beta) + 0.41*cos(beta) = 0.7103
%
%  FEM source: magnetic_sim/ANSYS/main/ANSYS_data/long2016_hexapole_halfcut/coilN
%      coil1 -> coil1_pre_fine_mesh (matching half-cut run).
%
%  Outputs:
%    magnetic_sim/ANSYS/main/MATLAB_data/long2016_hexapole_halfcut/bs_matrix/Vout_Vin_1p572.mat
%    magnetic_sim/ANSYS/main/doc/Solve_B_matrix/long2016_hexapole_halfcut/scripts/Vout_Vin_1p572.tex
%  ------------------------------------------------------------------

clear; clc; close all;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');

cnst = mt_constants();

apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];

res_base = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';

%% --- Hall sensor geometry ---
sensor_diam_mm = 0.3;
sensor_radius  = sensor_diam_mm/2 * 1e-3;
n_rings        = 15;
ring_dr        = sensor_radius / n_rings;
disc_local     = [0, 0];
for k = 1:n_rings
    r_k = ring_dr * k;
    n_k = 6*k;
    phi = (0:n_k-1) * 2*pi / n_k;
    disc_local = [disc_local; r_k*cos(phi(:)), r_k*sin(phi(:))]; %#ok<AGROW>
end
N_disc_points = size(disc_local,1);
assert(N_disc_points == 721, 'disc point count mismatch');

%% --- Geometry ---
beta_deg = atan2d(3.0, 15.0);
beta     = deg2rad(beta_deg);

%% --- Sensor placement: 1.572 mm (CHANGED: baseline - 3.0 mm) ---
S_ALONG  = 1.572e-3;
S_OFFSET = 0.41e-3;

sensor_pl_lower = [S_ALONG; S_OFFSET; 0];
n_pl_lower      = [0; 1; 0];

sensor_pl_upper = [S_ALONG*cos(beta) - S_OFFSET*sin(beta); ...
                   S_ALONG*sin(beta) + S_OFFSET*cos(beta); ...
                   0];
n_pl_upper      = [-sin(beta); cos(beta); 0];

fprintf('Sensor placement (pole-local, mm):\n');
fprintf('  Lower: axial=%.4f  transverse=%.4f  (offset 0.41)\n', ...
        sensor_pl_lower(1)*1e3, sensor_pl_lower(2)*1e3);
fprintf('  Upper: axial=%.4f  transverse=%.4f  (= %.4f along slant + %.4f normal)\n', ...
        sensor_pl_upper(1)*1e3, sensor_pl_upper(2)*1e3, S_ALONG*1e3, S_OFFSET*1e3);

%% --- Build sensor positions and normals in ANSYS global coords ---
sensor_pos   = zeros(3, 6);
sensor_n     = zeros(3, 6);
disc_basis_u = zeros(3, 6);
disc_basis_v = zeros(3, 6);

for i = 1:6
    theta    = cnst.pole_angles(i) * pi/180;
    is_lower = cnst.pole_is_lower(i);

    if is_lower
        pole_axis = [cos(theta); sin(theta); 0];
        up_hat    = [0; 0; 1];
        tip_z     = -cnst.R_norm_z;
        s_pl      = sensor_pl_lower;
        n_pl      = n_pl_lower;
    else
        inc       = cnst.upper_incline;
        pole_axis = [cos(inc)*cos(theta); cos(inc)*sin(theta); sin(inc)];
        up_unnorm = [0; 0; 1] - sin(inc) * pole_axis;
        up_hat    = up_unnorm / norm(up_unnorm);
        tip_z     = +cnst.R_norm_z;
        s_pl      = sensor_pl_upper;
        n_pl      = n_pl_upper;
    end

    tip = [cnst.R_norm_xy*cos(theta); cnst.R_norm_xy*sin(theta); tip_z];

    side_hat = cross(pole_axis, up_hat);
    side_hat = side_hat / norm(side_hat);

    R = [pole_axis, up_hat, side_hat];

    sensor_pos(:, i) = tip + R * s_pl;
    sensor_n(:, i)   = R * n_pl;

    u = up_hat - dot(up_hat, sensor_n(:,i)) * sensor_n(:,i);
    if norm(u) < 1e-9
        u = pole_axis - dot(pole_axis, sensor_n(:,i)) * sensor_n(:,i);
    end
    u = u / norm(u);
    v = cross(sensor_n(:,i), u);
    disc_basis_u(:,i) = u;
    disc_basis_v(:,i) = v;
end

%% --- Load FEM for 6 coils, sample disc ---
fprintf('\nSampling 721-pt disc at each sensor for 6 coils...\n');

B_matrix_at_Iin = zeros(6, 6);

for k_apdl = 1:6
    coilname = sprintf('coil%d', k_apdl);
    if k_apdl == 1
        coil_dir = fullfile(res_base, 'coil1_pre_fine_mesh');
    else
        coil_dir = fullfile(res_base, coilname, 'standard');
    end
    fprintf('  %s (%s): ', coilname, coilname);

    d = import_ansys_data(coil_dir, 'all', coilname);
    z_wp = d.z - cnst.SPH_OFST;
    Fx = scatteredInterpolant(d.x, d.y, z_wp, d.bx, 'linear', 'nearest');
    Fy = scatteredInterpolant(d.x, d.y, z_wp, d.by, 'linear', 'nearest');
    Fz = scatteredInterpolant(d.x, d.y, z_wp, d.bz, 'linear', 'nearest');

    paper_j = apdl_to_paper_idx(k_apdl);
    for i = 1:6
        u = disc_basis_u(:,i);
        v = disc_basis_v(:,i);
        Bn_samples = zeros(N_disc_points, 1);
        for p = 1:N_disc_points
            xp = sensor_pos(:,i) + disc_local(p,1)*u + disc_local(p,2)*v;
            bxp = Fx(xp(1), xp(2), xp(3));
            byp = Fy(xp(1), xp(2), xp(3));
            bzp = Fz(xp(1), xp(2), xp(3));
            Bn_samples(p) = [bxp, byp, bzp] * sensor_n(:,i);
        end
        B_matrix_at_Iin(i, paper_j) = mean(Bn_samples);
    end
    fprintf('done (paper P%d)\n', paper_j);
end

%% --- Apply transfer chain ---
I_in_A = 0.6;
k_A    = 0.3614;
S_hall = 130;     % V/T (= 0.13 V/mT per EQ-730L datasheet) — FIXED unit bug

B_S   = (S_hall / I_in_A) * B_matrix_at_Iin;
V_VV  = k_A * B_S;
B_bar = V_VV;

fprintf('\n=== Diagonal |V_out/V_in| (mV/V): new (1.572) vs old (4.572) ===\n');
fprintf('  Pole   new (1.572)     old (4.572)         delta %% (abs)\n');
old_diag_mVV = [1.0978, 1.5763, 1.0600, 1.6551, 1.5723, 1.2005];
for i = 1:6
    new_mVV = V_VV(i,i);   % V/V (was *1e3 for old broken-unit mV/V display)
    delta   = (abs(new_mVV) - abs(old_diag_mVV(i))) / abs(old_diag_mVV(i)) * 100;
    fprintf('  P%d     %+8.4f       %+8.4f             %+6.2f %%\n', ...
            i, new_mVV, old_diag_mVV(i), delta);
end

%% --- Save .mat ---
data_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\bs_matrix';
if ~exist(data_dir, 'dir'); mkdir(data_dir); end
mat_path = fullfile(data_dir, 'Vout_Vin_1p572.mat');
sensor_along_mm = S_ALONG * 1e3;
sensor_offset_mm = S_OFFSET * 1e3;
save(mat_path, ...
    'V_VV', 'B_S', 'B_bar', 'B_matrix_at_Iin', ...
    'I_in_A', 'k_A', 'S_hall', ...
    'sensor_along_mm', 'sensor_offset_mm', ...
    'sensor_pos', 'sensor_n', 'disc_basis_u', 'disc_basis_v', ...
    'sensor_diam_mm', 'N_disc_points', ...
    'beta_deg', 'sensor_pl_lower', 'sensor_pl_upper');
fprintf('\nSaved .mat:   %s\n', mat_path);

%% --- LaTeX (bmatrix, mV/V, no x10^-3 factor) ---
V_VV_disp = V_VV;        % entries already V/V (after S_hall unit fix)
for i = 1:6
    for j = 1:6
        if i == j
            V_VV_disp(i,j) =  abs(V_VV_disp(i,j));
        else
            V_VV_disp(i,j) = -abs(V_VV_disp(i,j));
        end
    end
end

tex_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\doc\Solve_B_matrix\long2016_hexapole_halfcut\scripts';
if ~exist(tex_dir, 'dir'); mkdir(tex_dir); end
tex_file = fullfile(tex_dir, 'Vout_Vin_1p572.tex');

fid = fopen(tex_file, 'w');
oc  = onCleanup(@() fclose(fid));

fprintf(fid, '%% Auto-generated by magnetic_sim/ANSYS/main/matlab/long2016_hexapole_halfcut/fit/gen_Vout_Vin_1p572.m\n');
fprintf(fid, '%% Sensor placement: 1.572 mm along surface, 0.41 mm normal\n');
fprintf(fid, '\\documentclass[11pt]{article}\n');
fprintf(fid, '\\usepackage[margin=1in]{geometry}\n');
fprintf(fid, '\\usepackage{amsmath}\n');
fprintf(fid, '\\begin{document}\n\n');
fprintf(fid, '\\[\n\\frac{V_{out}}{V_{in}} = \\begin{bmatrix}\n');
for i = 1:6
    for j = 1:6
        fprintf(fid, '%+9.4f', V_VV_disp(i,j));
        if j < 6, fprintf(fid, ' & '); end
    end
    fprintf(fid, ' \\\\\n');
end
fprintf(fid, '\\end{bmatrix}\n\\]\n\n');
fprintf(fid, ['\\noindent $V_{out}/V_{in} = k_A\\,B_S = (k_A\\,S_{hall}/I_{in})\\,B_{\\mathrm{FEM}}$ ' ...
              '(entries in V/V, dimensionless).\\\\\n']);
fprintf(fid, ['Parameters: $k_A = %.4f$~A/V, $S_{hall} = %.0f$~V/T (= 0.13 V/mT, EQ-730L), $I_{in} = %.1f$~A.\\\\\n'], ...
              k_A, S_hall, I_in_A);
fprintf(fid, ['Sensor placement: %.3f mm along pole surface from tip + %.2f mm perpendicular ' ...
              '(0.3 mm $\\varnothing$ disc, 721-pt sampling).\\\\\n'], ...
              sensor_along_mm, sensor_offset_mm);
fprintf(fid, ['Row $i$ = $B_{\\mathrm{surface}}$ sensor on pole $i$;\\ ' ...
              'column $j$ = excited pole $j$;\\ order P1--P6.\\\\\n']);
fprintf(fid, 'Sign convention: diagonal $+|\\cdot|$, off-diagonal $-|\\cdot|$.\n');
fprintf(fid, '\n\\end{document}\n');

fprintf('Saved LaTeX: %s\n', tex_file);
