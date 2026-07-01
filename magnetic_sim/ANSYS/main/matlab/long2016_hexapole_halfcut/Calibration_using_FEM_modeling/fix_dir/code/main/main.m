%% main.m -- hexapole point-charge model (fix-l) calibration driver
%  Document : "Magnetic hexapole model.pdf"
%  Model    : B(p) = ^Bg_I * sum_i ( sum_j K_bar_ij I_j ) * (p/ell - dhat_i)/||p/ell - dhat_i||^3
%             charges fixed on the pole axes at pc_i = ell*dhat_i (no bias).
%  Fit vars : K_bar (6x6, K_bar(1,1) fixed = 5/6), ell, ^Bg_I (ghat_I_B)   (lsqnonlin)
%  Notation : paper «Lumped-Parameter…» — ^Bg_I=ghat_I_B (current gain, T/A), K_bar (1,1)=5/6.
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
%    3) fit_KI_fixl      -- lsqnonlin fit {K_bar, ell, ^Bg_I}       (fit)
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
VARIANT     = 'gap200um_mueq'; % [MODIFIED] FEM 變體子夾（'standard' = baseline；'gap200um_mueq' = gap200 2 段式 μ_eff）

%% ---- paths -----------------------------------------------------------------
TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\fix_dir'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis'); % mt_constants/import_ansys_data/filter_iron_nodes
addpath(fullfile(TREE,'code','function'));                                       % model helpers (added last -> takes precedence)
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
tex_dir      = fullfile(TREE,'results');
if ~exist(tex_dir,'dir'); mkdir(tex_dir); end
cal_dir      = fullfile(TREE,'data');   % 規則#2：.mat 放本組 data/（fit_fixl 解；Hall_sensor_base_fix_dir/decouple 由此載 ℓ̂）
if ~exist(cal_dir,'dir'); mkdir(cal_dir); end

%% ---- constants + pole-tip directions ---------------------------------------
cnst = mt_constants();
apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];                 % coil k excites this paper pole
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);                             % 3x6 unit directions

%% ---- load FEM once (all 6 coils) -------------------------------------------
fprintf('loading 6 coils (variant ''%s'') ...\n', VARIANT);     % [MODIFIED]
C = load_coils(results_root, cnst, apdl_to_paper_idx, VARIANT);  % [MODIFIED] variant subfolder
vtag = '';                                                       % [ADDED] output suffix (non-standard variant only)
if ~strcmp(VARIANT,'standard'), vtag = ['_' VARIANT]; end

%% ---- radius list by mode ---------------------------------------------------
switch lower(MODE)
    case 'single', R_um_list = R_single_um;
    case 'sweep',  R_um_list = R_start_um:R_step_um:R_end_um;
    otherwise,     error('MODE must be ''single'' or ''sweep''.');
end

%% ---- fit + emit LaTeX result per R -----------------------------------------
for R_um = R_um_list
    [coil, nmin]              = select_ball(C, R_um*1e-6);            % coil.p [m]、bfem [mT]
    [ell, ghat_I_B, K_bar, J] = fit_KI_fixl(coil, dhat, I_actual);   % 在 SI 公尺擬合；^Bg_I [mT/A]、ell [m]
    ell                       = ell * 1e6;                           % m → µm（此後 write/save 用 µm）
    errpct                    = region_field_err(coil, J);
    fname = fullfile(tex_dir, sprintf('fit_%s_R%03dum_%gA%s.tex', SHAPE, R_um, I_actual, vtag));  % [MODIFIED] vtag
    write_KI_tex(fname, SHAPE, R_um, I_actual, K_bar, ell, ghat_I_B, errpct);
    compile_tex_pdf(fname);                                          % [ADDED] 編成可看的 PDF（清 aux/log）
    % 存 fit_KI_fixl 解成 .mat（供 Hall_sensor_base_fix_dir 載入 ℓ̂；ell 為 [µm]）
    gB = ghat_I_B;  Khat = K_bar;   % alias：維持 .mat field 名 'gB'/'Khat' 與下游 loader 相容（範圍：不改 .mat field 名）
    save(fullfile(cal_dir, sprintf('fit_fixl_R%03dum%s.mat', R_um, vtag)), ...     % [MODIFIED] vtag
         'ell','gB','Khat','J','errpct','R_um','I_actual','SHAPE','VARIANT');
    fprintf('R=%3d um | nmin/coil=%6d | ell=%.2f µm | ^Bg_I=%.4e mT/A | err=%.2f%%  -> fit_fixl_R%03dum%s.mat\n', ...
            R_um, nmin, ell, ghat_I_B, errpct, R_um, vtag);
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
