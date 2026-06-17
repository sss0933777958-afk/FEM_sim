%% gen_Vout_Vin_4p572.m
%  ------------------------------------------------------------------
%  Compute V_out/V_in 6x6 matrix at baseline sensor 4.572 mm,
%  with CORRECTED unit chain (was: S_hall=0.130 V/T was 1000x too small).
%
%  Correct unit chain (per EQ-730L datasheet + ANSYS MKS):
%      ANSYS B field is in TESLA (T) since /UNITS, MKS is used
%      EQ-730L sensitivity = 0.13 V/mT  (= 130 V/T)
%      Therefore:
%          B[mT]    = B_FEM[T] * 1000              (T -> mT)
%          v_H[V]   = 0.13 V/mT * B[mT]            (sensor voltage)
%          B_S[V/A] = v_H / I_in                   (per Ampere DC gain)
%          V_VV     = k_A * B_S                    (V_out/V_in dimensionless)
%
%  Old vs new (factor 1000 difference):
%      Old: S_hall=130e-3 V/T directly with B[T]  -> V_VV in mV/V (mislabeled)
%      New: S_hall=0.13 V/mT with B converted to mT -> V_VV in true V/V
%
%  Sensor placement: 4.572 mm baseline (the original LongFei spec)
%      LOWER (P1/P3/P6): pole-local (4.572, 0.41, 0) mm
%      UPPER (P2/P4/P5): pole-local decomposition of
%                          4.572 along cone slant + 0.41 outward normal
%                          axial      = 4.572*cos(beta) - 0.41*sin(beta) = 4.4051
%                          transverse = 4.572*sin(beta) + 0.41*cos(beta) = 1.2994
%
%  Outputs:
%    magnetic_sim/ANSYS/main/MATLAB_data/long2016_hexapole_halfcut/bs_matrix/Vout_Vin_4p572.mat
%    magnetic_sim/ANSYS/main/doc/Solve_B_matrix/long2016_hexapole_halfcut/scripts/Vout_Vin_4p572.tex
%  ------------------------------------------------------------------

clear; clc; close all;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');

cnst = mt_constants();

apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];

res_base = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut';

%% --- Hall sensor geometry (disc) ---
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

%% --- Sensor placement: 4.572 mm baseline ---
S_ALONG  = 4.572e-3;
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

%% --- Build sensor positions and normals ---
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

%% --- Load FEM (ANSYS B field in Tesla) ---
fprintf('\nSampling 721-pt disc at each sensor for 6 coils...\n');

B_matrix_at_Iin_T = zeros(6, 6);   % in TESLA (ANSYS MKS native)

for k_apdl = 1:6
    coilname = sprintf('coil%d', k_apdl);
    coil_dir = fullfile(res_base, coilname, 'standard');                  % [V4] Long2016 verbatim baseline, May 29
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
        B_matrix_at_Iin_T(i, paper_j) = mean(Bn_samples);
    end
    fprintf('done (paper P%d)\n', paper_j);
end

%% --- Per-pole sign correction (diag → all positive) ---
%  For each pole j, if B_matrix(j,j) < 0 in raw FEM, multiply column j by -1.
%  This redefines "+I" per pole such that diagonal = "B aligned with sensor n+"
%  is always positive. Physically equivalent to: "+I_j produces flux flowing
%  yoke → tip_j → WP" for all 6 poles.
B_matrix_at_Iin_T_raw = B_matrix_at_Iin_T;     % keep raw
col_sign              = sign(diag(B_matrix_at_Iin_T))';
B_matrix_at_Iin_T     = B_matrix_at_Iin_T .* col_sign;

fprintf('\nRaw FEM diag signs: %s\n', mat2str(col_sign));
fprintf('Columns sign-flipped: %s\n', mat2str(find(col_sign < 0)));

%% --- UNIT CHAIN (S_hall = 130 V/T directly on ANSYS raw B[T]) ---
% Simplified: skip T->mT conversion, multiply raw B[T] by 130 V/T
% (= equivalent to 0.13 V/mT * B[mT], but cleaner code path)
S_hall = 130;                    % V/T (= 0.13 V/mT, EQ-730L datasheet)
v_H    = S_hall * B_matrix_at_Iin_T;   % V (Hall sensor voltage, sign-corrected)

% Also keep mT view for console clarity
B_matrix_mT = B_matrix_at_Iin_T * 1000;

% Transfer chain
I_in_A = 1.0;                    % [V4] coil current (Long2016 verbatim baseline: CURR_ARRAY=1.0)
k_A    = 0.3614;                 % power amplifier gain (A/V)

B_S    = v_H / I_in_A;           % V/A (sensor voltage per coil current)
V_VV   = k_A * B_S;              % V/V dimensionless (= V_out/V_in DC gain)

fprintf('\n=== Diagonal: FEM (corrected units) vs paper eq(7) ===\n');
fprintf('  Pole    B_FEM (mT)   v_H (V)    B_S (V/A)   V/V        paper B (V/A)\n');
paper_diag = [0.2816, 0.2659, 0.2795, 0.2348, 0.2537, 0.2302];
for i = 1:6
    fprintf('  P%d      %+7.3f      %+7.4f    %+7.4f    %+7.4f    %.4f\n', ...
            i, B_matrix_mT(i,i), v_H(i,i), B_S(i,i), V_VV(i,i), paper_diag(i));
end
fprintf('\nFEM mean |diag(V/A)| = %.4f   Paper mean |diag(V/A)| = %.4f   ratio FEM/Paper = %.2fx\n', ...
        mean(abs(diag(B_S))), mean(paper_diag), mean(abs(diag(B_S)))/mean(paper_diag));

%% --- Save .mat ---
data_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\bs_matrix';
if ~exist(data_dir, 'dir'); mkdir(data_dir); end
mat_path = fullfile(data_dir, 'Vout_Vin_4p572.mat');
sensor_along_mm  = S_ALONG * 1e3;
sensor_offset_mm = S_OFFSET * 1e3;
save(mat_path, ...
    'V_VV', 'B_S', 'v_H', 'B_matrix_mT', 'B_matrix_at_Iin_T', ...
    'B_matrix_at_Iin_T_raw', 'col_sign', ...
    'I_in_A', 'k_A', 'S_hall', ...
    'sensor_along_mm', 'sensor_offset_mm', ...
    'sensor_pos', 'sensor_n', 'disc_basis_u', 'disc_basis_v', ...
    'sensor_diam_mm', 'N_disc_points', ...
    'beta_deg', 'sensor_pl_lower', 'sensor_pl_upper');
fprintf('\nSaved .mat:   %s\n', mat_path);

%% --- LaTeX (bmatrix, V/V units, sign-corrected physical sign) ---
%  V_VV already has per-pole sign correction (col_sign applied above).
%  Display preserves the physical sign — no cosmetic override.
V_VV_disp = V_VV;

tex_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\doc\Solve_B_matrix\long2016_hexapole_halfcut\scripts';
if ~exist(tex_dir, 'dir'); mkdir(tex_dir); end
tex_file = fullfile(tex_dir, 'Vout_Vin_4p572.tex');

fid = fopen(tex_file, 'w');
oc  = onCleanup(@() fclose(fid));

fprintf(fid, '%% Auto-generated by magnetic_sim/ANSYS/main/matlab/long2016_hexapole_halfcut/fit/gen_Vout_Vin_4p572.m\n');
fprintf(fid, '%% Sensor placement: 4.572 mm along surface, 0.41 mm normal (baseline)\n');
fprintf(fid, '%% UNIT CHAIN (corrected): ANSYS B[T] -> B[mT] (x1000) -> v_H = 0.13 V/mT * B[mT]\n');
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
fprintf(fid, ['\\noindent $V_{out}/V_{in} = k_A\\,B_S = (k_A\\,S_{hall}/I_{in})\\,B_{\\mathrm{FEM,mT}}$ ' ...
              '(entries in V/V, dimensionless).\\\\\n']);
fprintf(fid, ['Parameters: $k_A = %.4f$~A/V, $S_{hall} = %g$~V/T (EQ-730L: $0.13$~V/mT $= 130$~V/T), $I_{in} = %.1f$~A.\\\\\n'], ...
              k_A, S_hall, I_in_A);
fprintf(fid, ['Unit chain: $v_H = S_{hall}\\cdot B_{\\mathrm{ANSYS}}[\\mathrm{T}]$, ' ...
              'then $V/V = k_A\\,v_H/I_{in}$.\\\\\n']);
fprintf(fid, ['Sensor placement: %.3f mm along pole surface from tip + %.2f mm perpendicular ' ...
              '(0.3 mm $\\varnothing$ disc, 721-pt sampling).\\\\\n'], ...
              sensor_along_mm, sensor_offset_mm);
fprintf(fid, ['Row $i$ = $B_{\\mathrm{surface}}$ sensor on pole $i$;\\ ' ...
              'column $j$ = excited pole $j$;\\ order P1--P6.\\\\\n']);
fprintf(fid, ['Sign convention: per-pole sign correction — column $j$ multiplied by ' ...
              '$\\mathrm{sign}(\\bar{B}_S(j,j))$ so all diagonals positive.\\\\\n']);
fprintf(fid, 'Columns sign-flipped vs raw FEM: %s.\n', mat2str(find(col_sign < 0)));
fprintf(fid, '\n\\end{document}\n');

fprintf('Saved LaTeX: %s\n', tex_file);
