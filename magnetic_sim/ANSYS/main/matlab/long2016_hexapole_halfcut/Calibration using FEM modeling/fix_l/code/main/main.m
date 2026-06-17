%% main.m -- hexapole point-charge model (fix-l) calibration driver
%  Document : "Magnetic hexapole model.pdf"
%  Model    : B(p) = gB * sum_i ( sum_j Khat_ij I_j ) * (p/ell - dhat_i)/||p/ell - dhat_i||^3
%             charges fixed on the pole axes at pc_i = ell*dhat_i (no bias).
%  Fit vars : Khat_I^FEM (6x6, Khat(1,1) fixed = 5/6), ell, gB   (lsqnonlin)
%  Current  : I = 1 A = FEM excitation current (per fit-current-matches-sim rule).
%
%  MODE switch:
%    'single' -> one fit at R = R_single_um
%    'sweep'  -> fits at R = R_start_um : R_step_um : R_end_um
%  Output: results-only LaTeX scripts (one per R) to results\fix_1\,
%          named  fit_<shape>_R<R>um_<I>A.tex   e.g. fit_ball_R150um_1A.tex
%  All model math lives in code\function\ ; this file is just the driver.
%
%  Pipeline (top-to-bottom call order):
%    1) load_coils      -- load the 6-coil FEM field once          (load data)
%    2) select_ball     -- keep nodes inside the sampling ball R   (pick region)
%    3) fit_KI_fixl      -- lsqnonlin fit {Khat, ell, gB}           (fit)
%    4) region_field_err -- relative RMS field error over region    (accuracy)
%    5) write_KI_tex     -- emit results-only LaTeX script           (output)

clear; clc;

%% ---- config ----------------------------------------------------------------
MODE        = 'single';        % 'single' | 'sweep'
R_single_um = 150;             % single-mode sampling-ball radius [um]
R_start_um  = 50;              % sweep start [um]
R_step_um   = 5;               % sweep step  [um]
R_end_um    = 500;             % sweep end   [um]
I_actual    = 1;               % drive current [A] = FEM excitation (1 A)
SHAPE       = 'ball';          % sampling region: ball ||p|| <= R about WP

%% ---- paths -----------------------------------------------------------------
TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration using FEM modeling\fix_l'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\hexapole-long2016\analysis'); % mt_constants/import_ansys_data/filter_iron_nodes
addpath(fullfile(TREE,'code','function'));                                       % model helpers (added last -> takes precedence)
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut';
tex_dir      = fullfile(TREE,'results','fix_l');
if ~exist(tex_dir,'dir'); mkdir(tex_dir); end

%% ---- constants + pole-tip directions ---------------------------------------
cnst = mt_constants();
apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];                 % coil k excites this paper pole
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);                             % 3x6 unit directions

%% ---- load FEM once (all 6 coils) -------------------------------------------
fprintf('loading 6 coils ...\n');
C = load_coils(results_root, cnst, apdl_to_paper_idx);

%% ---- radius list by mode ---------------------------------------------------
switch lower(MODE)
    case 'single', R_um_list = R_single_um;
    case 'sweep',  R_um_list = R_start_um:R_step_um:R_end_um;
    otherwise,     error('MODE must be ''single'' or ''sweep''.');
end

%% ---- fit + emit LaTeX result per R -----------------------------------------
for R_um = R_um_list
    [coil, nmin]       = select_ball(C, R_um*1e-6);
    [ell, gB, Khat, J] = fit_KI_fixl(coil, dhat, I_actual);
    errpct             = region_field_err(coil, J);
    fname = fullfile(tex_dir, sprintf('fit_%s_R%03dum_%gA.tex', SHAPE, R_um, I_actual));
    write_KI_tex(fname, SHAPE, R_um, I_actual, Khat, ell, gB, errpct);
    fprintf('R=%3d um | nmin/coil=%6d | ell=%.3f mm | gB=%.4e | err=%.2f%%\n', ...
            R_um, nmin, ell*1e3, gB, errpct);
end
fprintf('done (%s mode): %d .tex file(s) in %s\n', MODE, numel(R_um_list), tex_dir);
