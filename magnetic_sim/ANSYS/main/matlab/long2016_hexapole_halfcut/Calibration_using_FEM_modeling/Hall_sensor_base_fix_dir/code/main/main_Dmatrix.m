%% main_Dmatrix.m -- Hall-sensor-based hexapole flux model (single-parameter, on-axis) -- recover 6x6 D matrix
% =========================================================================
%  本包 = Hall_sensor_base_fix_dir = **single-parameter（電荷在軸 d̂、無 bias）** 模型。
%  （這支 D 矩陣 driver 過去誤用 no_fix bias 模型 → ℓ̂ 跟 fix_dir 不一致；已改回在軸。）
%
%  Model (single-parameter, on-axis): point charges at p_ck = ℓ̂·Pc_base（actuator 框，無 ê bias）。
%     b_ij = S(p_i/ℓ̂)·Dv_j,  Dv_j profiled per excitation: Dv_j=(ΣSᵀS)⁻¹(ΣSᵀb_ij)。
%     只有 1 個自由幾何參數 ℓ̂：minimize J(ℓ̂)=‖A(ℓ̂)·Dv − Bstack‖²（Dv 已 profile 掉）。
%     此 ℓ̂ 與 fix_dir/fit_KI_fixl、與本包 main.m/sensor_cost_lhat **同一條 cost** → ℓ̂ 必一致。
%     Recover D (Hall-sensor): D^v=[Dv_1..Dv_NI], V=Vmat → D=D^v Vᵀ(VVᵀ)⁻¹=(6·d11/5)·D̄。
%
%  REUSES existing functions (no duplicated helpers):
%     no_fix_dir/code/function : load_coils_actuator(+variant), select_ball, build_A
%     this pkg  /code/function : build_sensor_geometry, extract_Vmat_interp
%  ℓ̂ 用在軸 profiled field cost 一維搜尋（fminbnd），不再 fit_bias/make_Pc（那是離軸 bias）。
%
%  DATA (user spec): variant 'gap200um_mueq', coil1..6 (all-source) = GRADED mesh + support-base mu_eff.
%  Sensor voltage  : extract_Vmat_interp, 100-point cylinder average (GRADED sensor-local mesh CSV
%                    from MT_Export_SensorLocalMesh_graded.txt; node-IDs match gap200 .dat).
%  Ordering        : all matrices in paper P1..P6 (excitation cols reordered
%                    APDL coil1..6 -> P1..P6 so diagonal = self-excitation).
%  Output          : xelatex matrix PDF (D_bar, D, V; P1..P6 labelled) ->
%                    results/D_<variant>.pdf ; .mat -> this pkg's data/ (rule #2).
%  Current         : I = 1 A = FEM excitation (per fit-current-matches-sim rule).
% =========================================================================

clear; clc;

%% ---- config ----------------------------------------------------------------
VARIANT    = 'gap200um_mueq';  % FEM 變體：coil1..6 場來源（gap200 2 段式 μ_eff）
N_I        = 6;                % FEM 模擬次數 = 6 個單線圈解
R_select   = 150e-6;           % 取點半徑 [m]：擬合只用此球內 air 節點
I_actual   = 1;                % 驅動電流 [A] = FEM 激發電流（1 A）
S_hall     = 130;              % Hall 靈敏度 [V/T]（EQ-730L）
ELL_LO     = 0.2e-3;           % ℓ̂ 一維搜尋下界 [m]
ELL_HI     = 3.0e-3;           % ℓ̂ 一維搜尋上界 [m]
n_uniform  = 10000;            % sensor 圓柱內均勻取樣點數（內插平均；2026-06-30 定案 10k）
dataset    = 'all';            % 取全 air 節點，再用 R_select 選球

%% ---- paths（沿用既有函式，不另寫）-----------------------------------------
CAL  = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');     % mt_constants/import_ansys_data/filter_iron_nodes
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');  % ansys_path
addpath(fullfile(CAL,'no_fix_dir','code','function'));          % 沿用：load_coils_actuator/select_ball/build_A
addpath(fullfile(CAL,'Hall_sensor_base_fix_dir','code','function'));  % 沿用：build_sensor_geometry/extract_Vmat_interp（最後 addpath→優先）
model        = 'long2016_hexapole_halfcut';
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
mesh_csv_dir = fullfile(results_root,'mesh','graded','csv');   % gap200=graded → graded sensor-local tet CSV（csv/ 子夾）
out_dir      = fullfile(CAL,'Hall_sensor_base_fix_dir','results');
if ~exist(out_dir,'dir'); mkdir(out_dir); end

%% ---- 常數 + 慣例 -----------------------------------------------------------
cnst = mt_constants();
apdl_to_paper_idx = [1,3,6,5,2,4];

%% ---- 步驟1：載 6-coil FEM（gap200, actuator 框, all-source）→ 選 R 球 ------
fprintf('載 6 coils（variant ''%s'', dataset ''%s''）...\n', VARIANT, dataset);
D = load_coils_actuator(model, cnst, apdl_to_paper_idx, dataset, VARIANT);   % [沿用，+variant]
[P, Bstack, npts] = select_ball(D, R_select);                                % [沿用]
fprintf('選出 N_p = %d 個 air 節點（R ≤ %g µm 球內）\n', npts, R_select*1e6);

% ---- all-source（只翻 sink）：對 raw FEM B 直接「只翻下極 sink(P1/P3/P6)」、上極 source 不動 ----
%  load_coils_actuator 回傳的是「全域 −B_FEM」→ 先還原 raw（×−1），再套 s_sink（下極=−1、上極=+1）。
%  與 extract_Vmat_interp 的 exc_sign 同一條 literal flip-sink 邏輯（V 與 D^v 慣例一致）。
%  ℓ̂/field error 不受符號影響，只讓 D^v/D 全 source（對角全正）。
s_sink = ones(1,6);
for j = 1:6, if ismember(apdl_to_paper_idx(j), [1 3 6]), s_sink(j) = -1; end; end   % 只翻下極 sink
Bstack = (-Bstack) .* s_sink;            % 還原 raw（−(−B_FEM)=B_FEM）後只翻下極 → all-source

%% ---- 步驟2：擬合 ℓ̂（single-parameter，在軸 Pc_base，無 bias）--------------
%  在軸 profiled field cost J(ℓ̂)=‖A(ℓ̂)·Dv − Bstack‖²（Dv=(AᵀA)⁻¹AᵀBstack 已 profile 掉）；
%  fminbnd 一維搜尋。此 cost 與 fix_dir/fit_KI_fixl 完全同一條 → ℓ̂ 應 = fix_dir。
Pc            = D.Pc_base;                                                    % 在軸電荷位置（actuator 框理想格）
opt           = optimset('TolX',1e-9,'Display','off');
[ell_hat, J]  = fminbnd(@(l) onaxis_cost(l, P, Bstack, Pc), ELL_LO, ELL_HI, opt);
fprintf('擬合：ℓ̂ = %.4f mm（single-parameter, on-axis）| J = %.6e\n', ell_hat*1e3, J);
if abs(ell_hat-ELL_LO)<1e-9 || abs(ell_hat-ELL_HI)<1e-9
    warning('ℓ̂ 落在搜尋邊界 → 放寬 ELL_LO/ELL_HI 再確認。');
end

%% ---- 步驟3：profile 每次激發電荷 Dv_j → D^v（6×N_I）-----------------------
A   = build_A(ell_hat, Pc, P);          % [沿用] 3Np×6 stacked kernel（在軸 Pc）
M   = A.' * A;                           % = Σ_i Sᵀ S
Dv  = M \ (A.' * Bstack);                % Dv_j=(ΣSᵀS)⁻¹(ΣSᵀb_ij)，逐欄 j 一次解 (6×N_I)

%% ---- 步驟4：Hall sensor 電壓 Vmat（gap200、內插 100 點、all-source）---------
[sensor_pos, sensor_n] = build_sensor_geometry(cnst);                         % [沿用]
[Vmat, exc_sign] = extract_Vmat_interp(results_root, cnst, apdl_to_paper_idx, ... % [沿用，graded CSV]
                       sensor_pos, sensor_n, S_hall, mesh_csv_dir, n_uniform, [], [], VARIANT);

%% ---- 步驟5：解 D = D^v Vᵀ(VVᵀ)⁻¹、gauge D̄ --------------------------------
Dmat  = (Dv * Vmat.') / (Vmat * Vmat.');     % D = D^v Vᵀ(VVᵀ)⁻¹
D_bar = Dmat * (5 / (6 * Dmat(1,1)));        % gauge：D̄(1,1)=5/6

%% ---- 區域場誤差 + 重建檢查 -------------------------------------------------
errpct = 100 * sqrt(J / sum(Bstack(:).^2));
recon  = norm(Dmat*Vmat - Dv,'fro') / norm(Dv,'fro');

%% ---- 激發欄重排成 paper P1..P6（sensor 列、電荷列本就 P1..P6）---------------
[~, paper_to_apdl] = sort(apdl_to_paper_idx);   % = [1 5 2 6 4 3]：欄 k = paper Pk 對應的 APDL coil
Vmat_p = Vmat(:, paper_to_apdl);                % 列=sensor P1..P6；欄重排→激發 P1..P6（對角=自激發）
Dv_p   = Dv(:,   paper_to_apdl);                % 列=電荷 P1..P6；欄重排→激發 P1..P6

%% ---- 印結果（全 P1..P6）----------------------------------------------------
fprintf('\n========= Hall-sensor single-parameter（variant=%s, ℓ̂=%.3f mm）=========\n', VARIANT, ell_hat*1e3);
fprintf('  N_p=%d | J=%.4e | region err=%.3f%% | recon ‖D·V−D^v‖/‖D^v‖=%.2e\n', npts, J, errpct, recon);
fprintf('  D^v = [Dv_1..Dv_6]（列=電荷極 P1..P6，欄=激發 P1..P6；profiled 每激發電荷，單位 T）：\n');
for i=1:6, fprintf('   % .4e % .4e % .4e % .4e % .4e % .4e\n', Dv_p(i,:)); end
fprintf('  Vmat [V]（列=sensor P1..P6，欄=激發 P1..P6，對角=自激發）：\n');
for i=1:6, fprintf('   % .3e % .3e % .3e % .3e % .3e % .3e\n', Vmat_p(i,:)); end
fprintf('  D 矩陣（列=電荷極 P1..P6，欄=sensor P1..P6；Dv = D·V）：\n');
for i=1:6, fprintf('   % .4e % .4e % .4e % .4e % .4e % .4e\n', Dmat(i,:)); end
fprintf('  D̄ 矩陣（gauge D̄(1,1)=5/6；列=電荷 P1..P6，欄=sensor P1..P6）：\n');
for i=1:6, fprintf('   % .4f % .4f % .4f % .4f % .4f % .4f\n', D_bar(i,:)); end
fprintf('===================================================================\n');

%% ---- 存 .mat 到本組 data/（results/ 只留 PDF）------------------
mat_dir = fullfile(CAL,'Hall_sensor_base_fix_dir','data');   % 規則#2：.mat 放本組 data/
if ~exist(mat_dir,'dir'); mkdir(mat_dir); end
mat_out = fullfile(mat_dir, sprintf('calib_D_%s.mat', VARIANT));
save(mat_out, 'Dmat','D_bar','Dv_p','Vmat_p','exc_sign','ell_hat','Pc', ...
              'J','errpct','recon','S_hall','R_select','npts','VARIANT', ...
              'apdl_to_paper_idx','paper_to_apdl','sensor_pos','sensor_n','n_uniform');
fprintf('已存 %s\n', mat_out);

%% ---- 產生 xelatex 矩陣 PDF（D̄、D、V，皆 P1..P6 標列/欄）-------------------
pole     = {'P1','P2','P3','P4','P5','P6'};
tex_path = fullfile(out_dir, sprintf('D_%s.tex', VARIANT));
pdf_path = fullfile(out_dir, sprintf('D_%s.pdf', VARIANT));
fid = fopen(tex_path,'w');
fprintf(fid,'%% Auto-generated by main_Dmatrix.m — Hall-sensor single-parameter D matrix (P1..P6)\n');
fprintf(fid,'\\documentclass[11pt]{article}\n');
fprintf(fid,'\\usepackage[margin=1in]{geometry}\n');
fprintf(fid,'\\usepackage{amsmath}\n');
fprintf(fid,'\\begin{document}\n');
fprintf(fid,'\\begin{center}\\large Hall-sensor single-parameter (on-axis) calibration\\quad(variant: %s)\\end{center}\n\n', ...
        strrep(VARIANT,'_','\_'));
emit_labeled_matrix(fid, 'D^{v}', Dv_p, pole, pole, 'auto', ...
    'row $i$ = charge pole P$i$;\ column $j$ = excited pole P$j$;\ profiled per-excitation charges $D^{v}=[Dv_1\ldots Dv_6]$,\ $Dv_j=(\Sigma S^{\top}S)^{-1}(\Sigma S^{\top}b_{ij})$;\ all-source (flip-sink: lower P1/P3/P6);\ $D^{v}=D\,V$;\ units T.');
emit_labeled_matrix(fid, '\bar{D}', D_bar, pole, pole, '', ...
    'row $i$ = charge pole P$i$;\ column $j$ = sensor pole P$j$;\ gauge $\bar{D}_{11}=5/6$.');
emit_labeled_matrix(fid, 'D', Dmat, pole, pole, 'auto', ...
    'row $i$ = charge pole P$i$;\ column $j$ = sensor pole P$j$;\ $D=D^{v}V^{\top}(VV^{\top})^{-1}$,\ $D^{v}=D\,V$.');
emit_labeled_matrix(fid, 'V', Vmat_p, pole, pole, 'auto', ...
    'row $i$ = sensor pole P$i$;\ column $j$ = excited pole P$j$ (diagonal = self-excitation);\ units V.');
fprintf(fid,['\\noindent\\small Single-parameter (on-axis) model;\\ order P1--P6;\\ $I=1$\\,A (= FEM excitation);\\ ' ...
             'all-source (flip-sink: lower P1/P3/P6; $V$ same convention);\\ %d-pt interp;\\ ' ...
             '$\\hat{\\ell}=%.4f$\\,mm,\\ region err $=%.3f\\%%$,\\ recon $=%.2e$.\n'], ...
        n_uniform, ell_hat*1e3, errpct, recon);
fprintf(fid,'\\end{document}\n');
fclose(fid);

xelatex = 'C:\Users\Kuo\AppData\Local\Programs\MiKTeX\miktex\bin\x64\xelatex.exe';
old = cd(out_dir);
[st, sysout] = system(sprintf('"%s" -interaction=nonstopmode -halt-on-error "%s"', xelatex, tex_path));
cd(old);
if st ~= 0 || ~exist(pdf_path,'file')
    fprintf('%s\n', sysout);
    error('xelatex 編譯失敗（見上方輸出）。');
end
% 清中間檔（results/ 只留 PDF）
for ext = {'.tex','.aux','.log','.out'}
    f = fullfile(out_dir, sprintf('D_%s%s', VARIANT, ext{1}));
    if exist(f,'file'); delete(f); end
end
fprintf('已存 %s\n', pdf_path);

%% ---- 移除舊「bias 版」產物（D_v2_* / calib_D_v2_*；本包不該有 bias 結果）---
for f = { fullfile(out_dir, sprintf('D_v2_%s.pdf', VARIANT)), ...
          fullfile(out_dir, sprintf('D_v2_%s.tex', VARIANT)), ...
          fullfile(mat_dir, sprintf('calib_D_v2_%s.mat', VARIANT)) }
    if exist(f{1},'file'); delete(f{1}); fprintf('已移除舊 bias 版 %s\n', f{1}); end
end

%% ---- local：在軸 single-parameter profiled field cost（Dv profile 掉）------
function J = onaxis_cost(l, P, Bstack, Pc)
    A  = build_A(l, Pc, P);                 % 3Np×6 在軸 kernel
    Dv = (A.'*A) \ (A.'*Bstack);            % 逐激發 profile 電荷
    R  = A*Dv - Bstack;                     % 殘差
    J  = sum(R(:).^2);                      % ‖A·Dv − Bstack‖²
end

%% ---- local function：渲染帶 P1..P6 表頭的矩陣（含自動 10^n 因子）----------
function emit_labeled_matrix(fid, name_tex, M, rowlab, collab, factor_mode, caption)
    if strcmp(factor_mode,'auto')
        mx = max(abs(M(:)));
        if mx > 0, e = floor(log10(mx)); else, e = 0; end
        if e ~= 0, Ms = M / 10^e; fac = sprintf('10^{%d}\\,', e);
        else,      Ms = M;        fac = '';
        end
    else
        Ms = M; fac = '';
    end
    fprintf(fid,'\\[\n%s = %s\\begin{array}{c|cccccc}\n', name_tex, fac);
    fprintf(fid,' ');
    for j = 1:6, fprintf(fid,'& %s ', collab{j}); end
    fprintf(fid,'\\\\\\hline\n');
    for i = 1:6
        fprintf(fid,'%s ', rowlab{i});
        for j = 1:6, fprintf(fid,'& %+9.4f ', Ms(i,j)); end
        fprintf(fid,'\\\\\n');
    end
    fprintf(fid,'\\end{array}\n\\]\n');
    fprintf(fid,'\\noindent\\small %s\n\n', caption);
end
