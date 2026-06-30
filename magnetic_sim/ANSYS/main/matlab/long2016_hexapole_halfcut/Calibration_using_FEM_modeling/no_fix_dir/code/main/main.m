%% main.m -- hexapole point-charge model (no_fix_l, 18-param bias) calibration driver
%  Document : "no_fix_l.pdf"
%  Model    : charges at pc = ell*(Pc_base + E(e_hat)), actuator frame; 18 params
%             = characteristic length ell + 1x17 bias e (e6z constrained).
%             Per simulation the 6 charge magnitudes g_j are profiled out (LS);
%             gauge with k-bar_I(1,1)=5/6 -> gB = (6/5)h11, K-bar_I = (5/(6 g11))H.
%  Current  : I = 1 A = FEM excitation current (per fit-current-matches-sim rule).
%
%  MODE switch:
%    'single' -> one fit at R = R_single_um
%    'sweep'  -> fits at R = R_start_um : R_step_um : R_end_um
%  Output: results-only LaTeX scripts (one per R) to results\no_fix_l\,
%          named  fit_<shape>_R<R>um_<I>A.tex   e.g. fit_ball_R150um_1A.tex
%  All model math lives in code\function\ ; this file is just the driver.
%
%  Pipeline (top-to-bottom call order):
%    1) load_coils_actuator -- load 6-coil FEM, rotate to actuator frame  (load data)
%    2) select_ball         -- keep nodes inside the sampling ball R       (pick region)
%    3) fit_bias            -- lsqnonlin fit {ell, e_hat(17)}              (fit)
%    4) gauge_KI            -- profile g_j, gauge -> gB, K-bar_I           (gauge)
%    5) region_field_err    -- relative RMS field error over region        (accuracy)
%    6) write_KbarI_tex     -- emit results-only LaTeX script               (output)

clear; clc;

%% ---- config ----------------------------------------------------------------
MODE        = 'single';        % 'single' | 'sweep'
R_single_um = 150;             % single-mode sampling-ball radius [um]
R_start_um  = 50;              % sweep start [um]
R_step_um   = 5;               % sweep step  [um]
R_end_um    = 500;             % sweep end   [um]
I_actual    = 1;               % drive current [A] = FEM excitation (1 A)
SHAPE       = 'ball';          % sampling region: ball ||p|| <= R about WP
ell0        = 0.5e-3;          % ell initial guess [m] (= ell_design)
dataset     = 'all';           % standard-mesh dataset
VARIANT     = 'gap200um_mueq'; % [MODIFIED] FEM 變體子夾（'standard' = baseline；'gap200um_mueq' = gap200 2 段式 μ_eff）

%% ---- paths -----------------------------------------------------------------
TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\no_fix_dir'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');         % ansys_path
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');   % mt_constants/import_ansys_data/filter_iron_nodes
addpath(fullfile(TREE,'code','function'));                                          % model helpers (last -> precedence)
model   = 'long2016_hexapole_halfcut';
tex_dir = fullfile(TREE,'results');
if ~exist(tex_dir,'dir'); mkdir(tex_dir); end

%% ---- constants + conventions -----------------------------------------------
cnst = mt_constants();
apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];                 % coil k excites this paper pole
coil_sign = [1 -1 1 -1 -1 1];                           % all-source display flip (upper P2/P4/P5)

%% ---- load FEM once (all 6 coils, actuator frame) ---------------------------
fprintf('loading 6 coils (variant ''%s'', dataset ''%s'') ...\n', VARIANT, dataset);   % [MODIFIED]
D = load_coils_actuator(model, cnst, apdl_to_paper_idx, dataset, VARIANT);              % [MODIFIED] variant
vtag = '';                                                                              % [ADDED] output suffix
if ~strcmp(VARIANT,'standard'), vtag = ['_' VARIANT]; end

%% ---- radius list by mode ---------------------------------------------------
switch lower(MODE)
    case 'single', R_um_list = R_single_um;
    case 'sweep',  R_um_list = R_start_um:R_step_um:R_end_um;
    otherwise,     error('MODE must be ''single'' or ''sweep''.');
end

%% ---- fit + emit LaTeX result per R -----------------------------------------
for R_um = R_um_list
    [P, Bstack, npts] = select_ball(D, R_um*1e-6);
    [ell, e_hat, J]   = fit_bias(P, Bstack, D.Pc_base, ell0);
    Pc                = make_Pc(e_hat, D.Pc_base);
    [KbarI, gB]       = gauge_KI(ell, Pc, P, Bstack, D.F);
    errpct            = region_field_err(Bstack, J);
    fname = fullfile(tex_dir, sprintf('fit_%s_R%03dum_%gA%s.tex', SHAPE, R_um, I_actual, vtag));  % [MODIFIED] vtag
    write_KbarI_tex(fname, SHAPE, R_um, I_actual, KbarI, ell, gB, e_hat, coil_sign, errpct);
    compile_tex_pdf(fname);                                          % [ADDED] 編成可看的 PDF（清 aux/log）
    fprintf('R=%3d um | npts=%6d | ell=%.3f mm | gB=%.4e | err=%.2f%%\n', ...
            R_um, npts, ell*1e3, gB, errpct);
end
fprintf('done (%s mode, variant=%s): %d result PDF(s) in %s\n', MODE, VARIANT, numel(R_um_list), tex_dir);

%% ---- local: 編 standalone .tex -> PDF（同夾，清中間檔；results/ 留 .pdf + .tex）----
function compile_tex_pdf(texfile)
    xelatex = 'C:\Users\Kuo\AppData\Local\Programs\MiKTeX\miktex\bin\x64\xelatex.exe';
    [d,b]   = fileparts(texfile);
    old = cd(d);
    [st,out] = system(sprintf('"%s" -interaction=nonstopmode -halt-on-error "%s"', xelatex, texfile));
    cd(old);
    if st ~= 0 || ~exist(fullfile(d,[b '.pdf']),'file')
        fprintf('%s\n', out); error('xelatex 編譯失敗：%s', texfile);
    end
    for ext = {'.tex','.aux','.log','.out'}    % results 只留 .pdf（連 .tex 一起清）
        f = fullfile(d,[b ext{1}]); if exist(f,'file'); delete(f); end
    end
end
