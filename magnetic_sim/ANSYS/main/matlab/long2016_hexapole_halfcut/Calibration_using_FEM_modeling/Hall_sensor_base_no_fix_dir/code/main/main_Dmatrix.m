%% main_Dmatrix.m -- Hall-sensor-based hexapole flux model (18-param bias, off-axis) -- recover 6x6 D matrix
% =========================================================================
%  本包 = Hall_sensor_base_no_fix_dir = **18-param bias（電荷離軸 ℓ̂·(Pc_base+E(ê))）** 模型。
%  與 Hall_sensor_base_fix_dir/main_Dmatrix.m（single-parameter 在軸）並列：唯一差別 = ℓ̂-fit。
%
%  Model (no-fix-l, 18-param bias): point charges at p_ck = ℓ̂·(Pc_base + E(ê))，actuator 框。
%     b_ij = S(p_i/ℓ̂, ê)·Dv_j,  Dv_j profiled per excitation: Dv_j=(ΣSᵀS)⁻¹(ΣSᵀb_ij)。
%     擬合 {ℓ̂, ê(17)}：minimize J(ℓ̂,ê)=‖A·Dv − Bstack‖²（Dv 已 profile；fit_bias/lsqnonlin）。
%     Recover D (Hall-sensor): D^v=[Dv_1..Dv_NI], V=Vmat → D=D^v Vᵀ(VVᵀ)⁻¹=(6·d11/5)·D̄。
%
%  REUSES existing functions (no duplicated helpers):
%     no_fix_dir/code/function : load_coils_actuator(+variant), select_ball, fit_bias, make_Pc, build_A
%     Hall_sensor_base_fix_dir/code/function : build_sensor_geometry, extract_Vmat_interp
%
%  DATA : variant 'gap200um_mueq'（GRADED mesh + support-base μ_eff），coil1..6（all-source）。
%  Sensor voltage : extract_Vmat_interp，10000-point cylinder average（GRADED sensor-local mesh CSV，
%                   修正後 sensor 位置；node-ID 對齊 gap200 .dat）。
%  Sign : all-source（literal flip-sink：只翻下極 sink P1/P3/P6；上極 source 不翻）→ D^v 對角全正。
%  Ordering : 所有矩陣 paper P1..P6（激發欄重排 APDL coil1..6 → P1..P6，對角=自激發）。
%  Output : xelatex PDF（D^v, D̄, D, V, ê）-> results/D_<variant>.pdf ; .mat -> 本包 data/。
%  Current : I = 1 A = FEM 激發（per fit-current-matches-sim）。
% =========================================================================

clear; clc;

%% ---- config ----------------------------------------------------------------
VARIANT    = 'gap200um_mueq';  % FEM 變體：coil1..6 場來源（gap200 2 段式 μ_eff）
N_I        = 6;                % FEM 模擬次數 = 6 個單線圈解
R_select   = 150e-6;           % 取點半徑 [m]：擬合只用此球內 air 節點
I_actual   = 1;                % 驅動電流 [A] = FEM 激發電流（1 A）
S_hall     = 130;              % Hall 靈敏度 [V/T]（EQ-730L）
ell0       = 0.5e-3;           % ℓ̂ 初值 [m]（= ℓ_design；fit_bias 用）
n_uniform  = 10000;            % sensor 圓柱內均勻取樣點數（內插平均；2026-06-30 定案 10k）
dataset    = 'all';            % 取全 air 節點，再用 R_select 選球

%% ---- paths（沿用既有函式，不另寫）-----------------------------------------
CAL  = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');     % mt_constants/import_ansys_data/filter_iron_nodes
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');  % ansys_path
addpath(fullfile(CAL,'no_fix_dir','code','function'));          % 沿用：load_coils_actuator/select_ball/fit_bias/make_Pc/build_A
addpath(fullfile(CAL,'Hall_sensor_base_fix_dir','code','function'));  % 沿用：build_sensor_geometry/extract_Vmat_interp（最後 addpath→優先）
model        = 'long2016_hexapole_halfcut';
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
mesh_csv_dir = fullfile(results_root,'mesh','graded','csv');   % gap200=graded → graded sensor-local tet CSV（csv/ 子夾）
TREE         = fullfile(CAL,'Hall_sensor_base_no_fix_dir');    % 本包（bias）
out_dir      = fullfile(TREE,'results');   if ~exist(out_dir,'dir'); mkdir(out_dir); end
mat_dir      = fullfile(TREE,'data');      if ~exist(mat_dir,'dir'); mkdir(mat_dir); end

%% ---- 常數 + 慣例 -----------------------------------------------------------
cnst = mt_constants();
apdl_to_paper_idx = [1,3,6,5,2,4];

%% ---- 步驟1：載 6-coil FEM（gap200, actuator 框）→ 選 R 球 ------------------
fprintf('載 6 coils（variant ''%s'', dataset ''%s''）...\n', VARIANT, dataset);
D = load_coils_actuator(model, cnst, apdl_to_paper_idx, dataset, VARIANT);   % [沿用，+variant]
[P, Bstack, npts] = select_ball(D, R_select);                                % [沿用]
fprintf('選出 N_p = %d 個 air 節點（R ≤ %g µm 球內）\n', npts, R_select*1e6);

% ---- all-source（只翻 sink）：對 raw FEM B 直接「只翻下極 sink(P1/P3/P6)」、上極 source 不動 ----
%  load_coils_actuator 回傳「全域 −B_FEM」→ 先還原 raw（×−1），再套 s_sink（下極=−1、上極=+1）。
%  與 extract_Vmat_interp 的 exc_sign 同一條 literal flip-sink 邏輯（V 與 D^v 慣例一致）。
%  ℓ̂/ê/field error 不受符號影響，只讓 D^v/D 全 source（對角全正）。
s_sink = ones(1,6);
for j = 1:6, if ismember(apdl_to_paper_idx(j), [1 3 6]), s_sink(j) = -1; end; end   % 只翻下極 sink
Bstack = (-Bstack) .* s_sink;            % 還原 raw（−(−B_FEM)=B_FEM）後只翻下極 → all-source

%% ---- 步驟2：擬合 ℓ̂、ê（18-param bias；profile Dv；lsqnonlin）---------------
[ell_hat, e_hat, J] = fit_bias(P, Bstack, D.Pc_base, ell0);                   % [沿用 no_fix_dir]
Pc = make_Pc(e_hat, D.Pc_base);                                              % [沿用] 離軸電荷位置
fprintf('擬合：ℓ̂ = %.4f mm（18-param bias）| J = %.6e | ‖ê‖ = %.4e\n', ell_hat*1e3, J, norm(e_hat));

%% ---- 步驟3：profile 每次激發電荷 Dv_j → D^v（6×N_I）-----------------------
A   = build_A(ell_hat, Pc, P);          % [沿用] 3Np×6 stacked kernel（離軸 Pc）
M   = A.' * A;                           % = Σ_i Sᵀ S
Dv  = M \ (A.' * Bstack);                % Dv_j=(ΣSᵀS)⁻¹(ΣSᵀb_ij)，逐欄 j 一次解 (6×N_I)

%% ---- 步驟4：Hall sensor 電壓 Vmat（gap200、內插 10000 點、all-source）---------
[sensor_pos, sensor_n] = build_sensor_geometry(cnst);                         % [沿用 fix_dir，修正位置]
[Vmat, exc_sign] = extract_Vmat_interp(results_root, cnst, apdl_to_paper_idx, ... % [沿用，graded CSV]
                       sensor_pos, sensor_n, S_hall, mesh_csv_dir, n_uniform, [], [], VARIANT);

%% ---- 步驟5：解 D = D^v Vᵀ(VVᵀ)⁻¹、gauge D̄ --------------------------------
Dmat  = (Dv * Vmat.') / (Vmat * Vmat.');     % D = D^v Vᵀ(VVᵀ)⁻¹
D_bar = Dmat * (5 / (6 * Dmat(1,1)));        % gauge：D̄(1,1)=5/6

%% ---- 區域場誤差 + 重建檢查 -------------------------------------------------
errpct = 100 * sqrt(J / sum(Bstack(:).^2));
recon  = norm(Dmat*Vmat - Dv,'fro') / norm(Dv,'fro');

%% ---- 激發欄重排成 paper P1..P6 + 重建 ê(3×6) -------------------------------
[~, paper_to_apdl] = sort(apdl_to_paper_idx);   % = [1 5 2 6 4 3]
Vmat_p = Vmat(:, paper_to_apdl);                % 列=sensor P1..P6；欄=激發 P1..P6
Dv_p   = Dv(:,   paper_to_apdl);                % 列=電荷 P1..P6；欄=激發 P1..P6
E36 = zeros(3,6);                               % 17-vec ê → 3×6（含約束 e6z），欄=P1..P6
E36(:,1)=e_hat(1:3); E36(:,2)=e_hat(4:6); E36(:,3)=e_hat(7:9); E36(:,4)=e_hat(10:12); E36(:,5)=e_hat(13:15);
E36(1,6)=e_hat(16);  E36(2,6)=e_hat(17);  E36(3,6)=e_hat(1)-e_hat(4)+e_hat(8)-e_hat(11)+e_hat(15);

%% ---- 印結果（全 P1..P6）----------------------------------------------------
fprintf('\n========= Hall-sensor 18-param bias（variant=%s, ℓ̂=%.3f mm）=========\n', VARIANT, ell_hat*1e3);
fprintf('  N_p=%d | J=%.4e | region err=%.3f%% | recon ‖D·V−D^v‖/‖D^v‖=%.2e\n', npts, J, errpct, recon);
fprintf('  D^v = [Dv_1..Dv_6]（列=電荷極 P1..P6，欄=激發 P1..P6；單位 T）：\n');
for i=1:6, fprintf('   % .4e % .4e % .4e % .4e % .4e % .4e\n', Dv_p(i,:)); end
fprintf('  Vmat [V]（列=sensor P1..P6，欄=激發 P1..P6，對角=自激發）：\n');
for i=1:6, fprintf('   % .3e % .3e % .3e % .3e % .3e % .3e\n', Vmat_p(i,:)); end
fprintf('  D̄ 矩陣（gauge D̄(1,1)=5/6；列=電荷 P1..P6，欄=sensor P1..P6）：\n');
for i=1:6, fprintf('   % .4f % .4f % .4f % .4f % .4f % .4f\n', D_bar(i,:)); end
fprintf('  ê(3×6, 列=x/y/z, 欄=P1..P6, 單位=ℓ̂ 比例)：\n');
for r=1:3, fprintf('   % .4f % .4f % .4f % .4f % .4f % .4f\n', E36(r,:)); end
fprintf('===================================================================\n');

%% ---- 存 .mat 到本包 data/（規則#2）----------------------------------------
mat_out = fullfile(mat_dir, sprintf('calib_D_%s.mat', VARIANT));
save(mat_out, 'Dmat','D_bar','Dv_p','Vmat_p','exc_sign','ell_hat','e_hat','E36','Pc', ...
              'J','errpct','recon','S_hall','R_select','npts','VARIANT', ...
              'apdl_to_paper_idx','paper_to_apdl','sensor_pos','sensor_n','n_uniform');
fprintf('已存 %s\n', mat_out);

%% ---- 產生 xelatex PDF（D^v、D̄、D、V、ê）---------------------------------
pole     = {'P1','P2','P3','P4','P5','P6'};
tex_path = fullfile(out_dir, sprintf('D_%s.tex', VARIANT));
pdf_path = fullfile(out_dir, sprintf('D_%s.pdf', VARIANT));
fid = fopen(tex_path,'w');
fprintf(fid,'%% Auto-generated by main_Dmatrix.m — Hall-sensor 18-param bias D matrix (P1..P6)\n');
fprintf(fid,'\\documentclass[11pt]{article}\n\\usepackage[margin=1in]{geometry}\n\\usepackage{amsmath}\n');
fprintf(fid,'\\begin{document}\n');
fprintf(fid,'\\begin{center}\\large Hall-sensor 18-param bias (off-axis) calibration\\quad(variant: %s)\\end{center}\n\n', ...
        strrep(VARIANT,'_','\_'));
emit_labeled_matrix(fid, 'D^{v}', Dv_p, pole, pole, 'auto', ...
    'row $i$ = charge pole P$i$;\ column $j$ = excited pole P$j$;\ profiled per-excitation charges $D^{v}=[Dv_1\ldots Dv_6]$;\ all-source (flip-sink: lower P1/P3/P6);\ $D^{v}=D\,V$;\ units T.');
emit_labeled_matrix(fid, '\bar{D}', D_bar, pole, pole, '', ...
    'row $i$ = charge pole P$i$;\ column $j$ = sensor pole P$j$;\ gauge $\bar{D}_{11}=5/6$.');
emit_labeled_matrix(fid, 'D', Dmat, pole, pole, 'auto', ...
    'row $i$ = charge pole P$i$;\ column $j$ = sensor pole P$j$;\ $D=D^{v}V^{\top}(VV^{\top})^{-1}$,\ $D^{v}=D\,V$.');
emit_labeled_matrix(fid, 'V', Vmat_p, pole, pole, 'auto', ...
    'row $i$ = sensor pole P$i$;\ column $j$ = excited pole P$j$ (diagonal = self-excitation);\ units V.');
% --- ê bias（3×6）---
fprintf(fid,'\\[\n\\hat{\\mathbf{e}} = \\begin{array}{c|cccccc}\n');
fprintf(fid,' '); for j=1:6, fprintf(fid,'& %s ', pole{j}); end; fprintf(fid,'\\\\\\hline\n');
rl = {'e_x','e_y','e_z'};
for r=1:3
    fprintf(fid,'%s ', rl{r});
    for j=1:6, fprintf(fid,'& %+8.4f ', E36(r,j)); end
    fprintf(fid,'\\\\\n');
end
fprintf(fid,'\\end{array}\n\\]\n');
fprintf(fid,'\\noindent\\small bias offset $\\hat{e}$ (charge positions $p_{ck}=\\hat{\\ell}\\,(P_{c,\\mathrm{base}}+E(\\hat{e}))$;\\ units = fraction of $\\hat{\\ell}$;\\ $e_{6z}$ from constraint $e_{1x}-e_{2x}+e_{3y}-e_{4y}+e_{5z}$).\n\n');
fprintf(fid,['\\noindent\\small 18-param bias (off-axis) model;\\ order P1--P6;\\ $I=1$\\,A (= FEM excitation);\\ ' ...
             'all-source (flip-sink: lower P1/P3/P6; $V$ same convention);\\ %d-pt interp;\\ ' ...
             '$\\hat{\\ell}=%.4f$\\,mm,\\ region err $=%.3f\\%%$,\\ recon $=%.2e$.\n'], n_uniform, ell_hat*1e3, errpct, recon);
fprintf(fid,'\\end{document}\n');
fclose(fid);

xelatex = 'C:\Users\Kuo\AppData\Local\Programs\MiKTeX\miktex\bin\x64\xelatex.exe';
old = cd(out_dir);
[st, sysout] = system(sprintf('"%s" -interaction=nonstopmode -halt-on-error "%s"', xelatex, tex_path));
cd(old);
if st ~= 0 || ~exist(pdf_path,'file'); fprintf('%s\n', sysout); error('xelatex 編譯失敗（見上方輸出）。'); end
for ext = {'.tex','.aux','.log','.out'}    % results 只留 .pdf
    f = fullfile(out_dir, sprintf('D_%s%s', VARIANT, ext{1})); if exist(f,'file'); delete(f); end
end
fprintf('已存 %s\n', pdf_path);

%% ---- local function：渲染帶 P1..P6 表頭的 6×6 矩陣（含自動 10^n 因子；沿用 fix_dir）----
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
