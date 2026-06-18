%% compute_B_bar_matrix.m
%  ------------------------------------------------------------------
%  Compute the 6x6 B_bar transfer matrix per `dh_derivation.pdf`:
%      B_bar(i,j) = (130e-3 / I) * B(i,j)
%  where
%      I    = I_in / k_A           (control input, V)
%      I_in = 0.6 A                (coil current — user spec)
%      k_A  = 0.3614 A/V           (Long Fei amplifier gain)
%      B(i,j) = magnitude of B along sensor i normal when coil j is
%               energised at I_in, AVERAGED over the 0.3 mm diameter
%               Hall sensor disc.
%
%  Sensor placement (design intent confirmed 2026-05-15):
%      LOWER (P1/P3/P6): pole-local (4.572, 0.41, 0) mm relative to tip
%                        — 4.572 along milled flat (= along axis),
%                          0.41 perpendicular up from flat (= +z in WP).
%      UPPER (P2/P4/P5): pole-local (4.403, 1.299, 0) mm relative to tip
%                        — 4.572 along the natural cone slant,
%                          0.41 perpendicular outward from cone face.
%
%  FEM source: magnetic_sim/ANSYS/main/ANSYS_data/long2016_hexapole_halfcut/coilN  (N = 1..6)
%  FEM current: TURNS=70, CURR_ARRAY(N)=0.6  (=> I_in_FEM = 0.6 A,
%  matches the experimental excitation directly — no linear scaling).
%
%  Outputs:
%    magnetic_sim/ANSYS/main/MATLAB_data/long2016_hexapole_halfcut/bs_matrix/B_bar.mat
%    magnetic_sim/ANSYS/main/doc/Solve B_matrix/scripts/B_bar_matrix_0p6A.tex
%  ------------------------------------------------------------------

clear; clc; close all;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');

cnst = mt_constants();

% APDL coil index k -> paper pole P-index (verbatim from compute_D_matrix.m)
apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];

res_base = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';

%% --- Hall sensor geometry ---
sensor_diam_mm   = 0.3;     % Hall sensor active area diameter (mm)
sensor_radius    = sensor_diam_mm/2 * 1e-3;    % 0.15e-3 m
% Build disc sample points in the sensor-local plane (z = 0 disc normal).
% Concentric-ring grid: centre (1 pt) + 15 rings at 0.01,0.02,...,0.15 mm
% (radial spacing 0.01 mm); ring k has 6k points (proportional to radius)
% so a plain mean approximates the disc area average.
n_rings    = 15;
ring_dr    = sensor_radius / n_rings;       % 0.01 mm radial spacing
disc_local = [0, 0];                        % centre
for k = 1:n_rings
    r_k = ring_dr * k;
    n_k = 6*k;                              % points proportional to radius
    phi = (0:n_k-1) * 2*pi / n_k;
    disc_local = [disc_local; r_k*cos(phi(:)), r_k*sin(phi(:))]; %#ok<AGROW>
end
N_disc_points = size(disc_local,1);         % 1 + 6*(1+...+15) = 721
assert(N_disc_points == 721, 'disc point count mismatch');

%% --- Geometry constants (CAD half-angle = 11.31 deg) ---
beta_deg = atan2d(3.0, 15.0);            % 11.31 deg  (CAD POLE_R / POLE_CONE_LEN)
beta     = deg2rad(beta_deg);

% Pole-local sensor offsets (relative to apex/tip), in metres
% LOWER pole: (4.572 along axis, 0.41 perpendicular)
sensor_pl_lower  = [4.572e-3; 0.41e-3; 0];
n_pl_lower       = [0; 1; 0];                            % normal = +y (perpendicular to milled flat)

% UPPER pole: (4.572 along cone slant, 0.41 perpendicular to slant outward)
%   foot  = 4.572 * (cos beta, sin beta, 0)
%   sensor = foot + 0.41 * (-sin beta, cos beta, 0)
%   normal = (-sin beta, cos beta, 0)
sensor_pl_upper  = [4.572e-3 * cos(beta) - 0.41e-3 * sin(beta); ...
                    4.572e-3 * sin(beta) + 0.41e-3 * cos(beta); ...
                    0];
n_pl_upper       = [-sin(beta); cos(beta); 0];

%% --- Build sensor positions and normals in ANSYS global coords ---
%  For each of the 6 poles, transform pole-local (sensor_pl, n_pl) to global.
%  Tip position from mt_constants. Pole axis direction and "up_hat" follow
%  the same convention as compute_D_matrix.m.

sensor_pos = zeros(3, 6);  % global ANSYS coords (m), tied to WP frame
sensor_n   = zeros(3, 6);

% Disc tangent basis storage (for later sampling)
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

    % Out-of-plane direction (perpendicular to pole_axis and up_hat)
    side_hat = cross(pole_axis, up_hat);
    side_hat = side_hat / norm(side_hat);

    % Rotation matrix from pole-local (x=axis, y=up_hat, z=side_hat) to global
    R = [pole_axis, up_hat, side_hat];     % 3x3

    sensor_pos(:, i) = tip + R * s_pl;
    sensor_n(:, i)   = R * n_pl;

    % Disc basis: 2 orthogonal unit vectors in the disc plane (perpendicular to n_hat)
    %   pick u along up_hat projected to perp(n); v = cross(n, u)
    u = up_hat - dot(up_hat, sensor_n(:,i)) * sensor_n(:,i);
    if norm(u) < 1e-9       % degenerate (n ~ up_hat); use pole_axis instead
        u = pole_axis - dot(pole_axis, sensor_n(:,i)) * sensor_n(:,i);
    end
    u = u / norm(u);
    v = cross(sensor_n(:,i), u);
    disc_basis_u(:,i) = u;
    disc_basis_v(:,i) = v;
end

fprintf('Sensor placement (WP frame, mm):\n');
fprintf('%-6s %3s %8s %8s %8s   %8s %8s %8s\n', 'Pole','lo','x','y','z','n_x','n_y','n_z');
for i = 1:6
    fprintf('P%d    %2d  %+8.4f %+8.4f %+8.4f   %+7.4f %+7.4f %+7.4f\n', i, ...
        cnst.pole_is_lower(i), sensor_pos(:,i)*1e3, sensor_n(:,i));
end

%% --- Load FEM for 6 coils, build interpolants, sample 0.3 mm disc ---
fprintf('\nLoading FEM data and sampling 0.3 mm disc at each sensor (%d points each)...\n', N_disc_points);

B_matrix_at_Iin = zeros(6, 6); % rows = sensor i (paper P), cols = excited coil paper-index. FEM was at I_in = 0.6 A directly.

for k_apdl = 1:6
    coilname = sprintf('coil%d', k_apdl);
    if k_apdl == 1
        % [FIXED 2026-05-21] the 'coil1' folder is an UNCUT-geometry run
        % (gives P1 self-act 0.0235 T, inconsistent with sibling lower poles
        % P3/P6 ~0.014). 'coil1_pre_fine_mesh' is the matching half-cut coil1
        % (same 05-15 batch as coil2-6; gives P1 = 0.01409 T). Verified by
        % diag_coil1_variants.m. File prefix inside that folder is 'coil1'.
        coil_dir = fullfile(res_base, 'coil1_pre_fine_mesh');
    else
        coil_dir = fullfile(res_base, coilname, 'standard');
    end
    fprintf('  %s (%s): ', coilname, coil_dir);

    d = import_ansys_data(coil_dir, 'all', coilname);
    z_wp = d.z - cnst.SPH_OFST;
    Fx = scatteredInterpolant(d.x, d.y, z_wp, d.bx, 'linear', 'nearest');
    Fy = scatteredInterpolant(d.x, d.y, z_wp, d.by, 'linear', 'nearest');
    Fz = scatteredInterpolant(d.x, d.y, z_wp, d.bz, 'linear', 'nearest');

    paper_j = apdl_to_paper_idx(k_apdl);
    for i = 1:6
        % Build 37 disc points in global coords
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

%% --- Compute B_bar at I_in = 0.6 A (FEM was run at this same current) ---
I_in_A    = 0.6;
k_A       = 0.3614;          % Long Fei amplifier gain (A/V)
S_hall    = 130e-3;          % Hall sensitivity (V/T)

% B_bar per derivation: B_bar = (130e-3 / I) * B,  where I = I_in / k_A
%   => B_bar = (130e-3 * k_A / I_in) * B(I_in)
B_bar = (S_hall * k_A / I_in_A) * B_matrix_at_Iin;

fprintf('\nB_matrix at I_in=%.2f A (Tesla, from FEM directly):\n', I_in_A);
print_matrix(B_matrix_at_Iin);
fprintf('\nB_bar (dimensionless) at I_in=%.2f A:\n', I_in_A);
print_matrix(B_bar);
fprintf('  |diag(B_bar)| mean = %.4e   (off-diag mean = %.4e)\n', ...
    mean(abs(diag(B_bar))), mean(abs(B_bar - diag(diag(B_bar))), 'all'));

%% --- Save .mat ---
data_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\bs_matrix';
if ~exist(data_dir, 'dir'); mkdir(data_dir); end
save(fullfile(data_dir, 'B_bar.mat'), ...
    'B_matrix_at_Iin', 'B_bar', ...
    'I_in_A', 'k_A', 'S_hall', ...
    'sensor_pos', 'sensor_n', 'disc_basis_u', 'disc_basis_v', ...
    'sensor_diam_mm', 'N_disc_points', ...
    'beta_deg', 'sensor_pl_lower', 'sensor_pl_upper');
fprintf('\nSaved: %s\n', fullfile(data_dir, 'B_bar.mat'));

%% --- Write LaTeX (just the matrix, booktabs + diagonal highlight) ---
tex_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\doc\Solve B_matrix\scripts';
tex_file = fullfile(tex_dir, 'B_bar_matrix_0p6A.tex');

write_bbar_tex(tex_file, B_bar, I_in_A, k_A);
fprintf('Saved: %s\n', tex_file);

%% --- Helpers ---
function print_matrix(M)
    fprintf('%-6s', '');
    for j = 1:size(M,2), fprintf(' P%-10d', j); end
    fprintf('\n');
    for i = 1:size(M,1)
        fprintf('P%d   ', i);
        for j = 1:size(M,2)
            fprintf('%+11.4e ', M(i, j));
        end
        fprintf('\n');
    end
end

function write_bbar_tex(filename, M, I_in_A, k_A)
    fid = fopen(filename, 'w');
    cleaner = onCleanup(@() fclose(fid));

    fprintf(fid, '%% Auto-generated by magnetic_sim/ANSYS/main/matlab/long2016_hexapole_halfcut/fit/compute_B_bar_matrix.m\n');
    fprintf(fid, '\\documentclass[11pt,a4paper]{article}\n');
    fprintf(fid, '\\usepackage[margin=1in]{geometry}\n');
    fprintf(fid, '\\usepackage{booktabs}\n');
    fprintf(fid, '\\usepackage{xcolor,colortbl}\n');
    fprintf(fid, '\\usepackage{amsmath}\n\n');
    fprintf(fid, '\\begin{document}\n\n');
    fprintf(fid, '\\section*{$\\bar{\\mathbf{B}}$ matrix — Long Fei half-cut hexapole \\\\\n');
    fprintf(fid, '($I_{in} = %.1f$~A, $k_A = %.4f$~A/V, sensor disc $\\varnothing 0.3$~mm)}\n\n', I_in_A, k_A);
    fprintf(fid, '\\noindent All entries below are multiplied by $10^{%d}$.\n\n', round(median(log10(abs(M(M~=0))))));
    fprintf(fid, '\\begin{table}[htbp]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\small\n');
    fprintf(fid, '\\begin{tabular}{c rrrrrr}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, ' & P1 & P2 & P3 & P4 & P5 & P6 \\\\\n');
    fprintf(fid, '\\midrule\n');
    % Pick a common exponent for the whole matrix so the table reads as
    % multiples of 10^p (mantissas only in cells).  Use the median |element|
    % order of magnitude to choose p; round mantissas to 4 sig figs.
    nz = M(M ~= 0);
    p  = round(median(log10(abs(nz))));
    scale = 10^p;
    for i = 1:6
        fprintf(fid, 'P%d', i);
        for j = 1:6
            mant = M(i,j) / scale;
            if i == j
                fprintf(fid, ' & \\cellcolor{yellow!20}$%+.4f$', mant);
            else
                fprintf(fid, ' & $%+.4f$', mant);
            end
        end
        fprintf(fid, ' \\\\\n');
    end
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n\n');
    fprintf(fid, '\\end{document}\n');
end
