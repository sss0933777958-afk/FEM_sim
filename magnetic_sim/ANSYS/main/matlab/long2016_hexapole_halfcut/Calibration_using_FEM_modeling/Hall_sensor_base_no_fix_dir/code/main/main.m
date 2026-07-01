%% main.m -- Hall-sensor 統一校正 driver（18-param bias, off-axis）：輸出 V, D̄, ^Bĝ_V, Ĥ_V
% =========================================================================
%  本包 = Hall_sensor_base_no_fix_dir = **18-param bias（電荷離軸 ℓ̂·(Pc_base+E(ê))）** 模型。
%  這支是兩包統一後的「唯一 main」（取代舊 main_Dmatrix）。與 Hall_sensor_base_fix_dir/main.m
%  （single-parameter 在軸）並列：唯一差別 = ℓ̂-fit（本包多 17 個 ê 離軸自由度，fit_bias/lsqnonlin）。
%  記號統一論文《Lumped-Parameter…》step 9（見 doc/charge_model_fitting/Current&Hall_sensor）。
%
%  Model (no-fix-l, 18-param bias): point charges at p_ck = ℓ̂·(Pc_base + E(ê))，actuator 框。
%     b_ij = S(p_i/ℓ̂, ê)·G_j,  G_j profiled per excitation: G_j=(ΣSᵀS)⁻¹(ΣSᵀb_ij)（code 變數 Dv，= G = D^v）。
%     擬合 {ℓ̂, ê(17)}：minimize J(ℓ̂,ê)=‖A·Dv − Bstack‖²（Dv 已 profile）。
%     Recover（論文 step 9）：Ĥ_V = G·Vᵀ(VVᵀ)⁻¹（code 變數 H_V）→ ^Bĝ_V=(6/5)h₁₁、D̄=(5/(6h₁₁))Ĥ_V。
%
%  sensor（使用者定案）：
%     位置 = build_sensor_geometry 預設（下極 −β 底錐面修正版）。
%     取值 = extract_Vmat_interp，Ø0.3mm × 0.1mm 圓柱（sensor_r=0.15e-3、axial_tol=0.10e-3 = 函式預設）、
%            均勻內插 n_uniform=10000 點平均。
%
%  主輸出（V, D̄, ^Bĝ_V, Ĥ_V；兩包同一組）：
%     ^Bĝ_V (ghat_V_B) = (6/5)·Ĥ_V(1,1) = 電壓側增益（T/V；Ĥ_V = ^Bĝ_V·D̄）。
%
%  論文↔code 變數：G=Dv（profiled charges=D^v）、Ĥ_V=H_V（舊 Dmat）、^Bĝ_V=ghat_V_B（舊 g_V）、D̄=D_bar。
%  REUSES（不重寫 helper）：
%     no_fix_dir/code/function : load_coils_actuator(+variant), select_ball, fit_bias, make_Pc, build_A
%     Hall_sensor_base_fix_dir/code/function : build_sensor_geometry, extract_Vmat_interp
%
%  DATA   : variant 'gap200um_mueq'（GRADED mesh + support-base μ_eff），coil1..6（all-source）。
%  Sign   : all-source（literal flip-sink：只翻下極 sink P1/P3/P6；上極 source 不翻）→ D^v 對角全正。
%  Order  : 所有矩陣 paper P1..P6（激發欄重排 APDL coil1..6 → P1..P6，對角=自激發）。
%  Output : xelatex PDF（V, D̄, ^Bĝ_V, Ĥ_V）-> results/D_<variant>.pdf；.mat -> 本包 data/（field 名不變，含 ê）。
%  Current: I = 1 A = FEM 激發（per fit-current-matches-sim）。
% =========================================================================

clear; clc;

%% ---- config ----------------------------------------------------------------
VARIANT    = 'gap200um_mueq';  % FEM 變體：coil1..6 場來源（gap200 2 段式 μ_eff）
N_I        = 6;                % FEM 模擬次數 = 6 個單線圈解
R_select   = 150e-6;           % 取點半徑 [m]：擬合只用此球內 air 節點
I_actual   = 1;                % 驅動電流 [A] = FEM 激發電流（1 A）
S_hall     = 130;              % Hall 靈敏度 [mV/mT]（EQ-730L；數值同 V/T，配 B[mT]→V[mV]）
ell0       = 0.5e-3;           % ℓ̂ 初值 [m]（= ℓ_design；fit_bias 在 SI 擬合、well-scaled）
n_uniform  = 10000;            % sensor 圓柱內均勻取樣點數（內插平均；定案 10k）
AXIAL_TOL  = 1e-6;             % [ADDED] sensor 圓柱厚度 [m]：canonical 0.10e-3；設 1e-6 = 1µm 薄片 what-if（檔名加 _h1um、不蓋 canonical；sensor_r 維持 Ø0.3mm）
htag       = '';  if abs(AXIAL_TOL-0.10e-3)>1e-12, htag = sprintf('_h%gum', AXIAL_TOL*1e6); end  % ≠canonical → 檔名後綴、不蓋原版
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
Bstack = (-Bstack) .* s_sink;            % 還原 raw（−(−B_FEM)=B_FEM）後只翻下極 → all-source（Bstack 已 mT）

%% ---- 步驟2：擬合 ℓ̂、ê（18-param bias；profile Dv；lsqnonlin，在 SI 公尺擬合）--
[ell_hat, e_hat, J] = fit_bias(P, Bstack, D.Pc_base, ell0);                   % [沿用 no_fix_dir]（P/ell 均公尺，well-scaled）
Pc = make_Pc(e_hat, D.Pc_base);                                              % [沿用] 離軸電荷位置（ê 無因次）
fprintf('擬合：ℓ̂ = %.2f µm（18-param bias）| J = %.6e | ‖ê‖ = %.4e\n', ell_hat*1e6, J, norm(e_hat));

%% ---- 步驟3：profile 每次激發電荷 Dv_j → D^v（6×N_I）-----------------------
A   = build_A(ell_hat, Pc, P);          % [沿用] 3Np×6 stacked kernel（離軸 Pc；ell/P 均公尺）
M   = A.' * A;                           % = Σ_i Sᵀ S
Dv  = M \ (A.' * Bstack);                % Dv_j=(ΣSᵀS)⁻¹(ΣSᵀb_ij)，逐欄 j 一次解 (6×N_I)
ell_hat = ell_hat * 1e6;                 % m → µm（此後 print/save/PDF 用 µm；build_A 已用完）

%% ---- 步驟4：Hall sensor 電壓 Vmat（gap200、Ø0.3mm 圓柱 × AXIAL_TOL 厚、內插 10k、all-source）---
[sensor_pos, sensor_n] = build_sensor_geometry(cnst);                         % [沿用 fix_dir，下極 −β 修正位置]
[Vmat, exc_sign] = extract_Vmat_interp(results_root, cnst, apdl_to_paper_idx, ... % [沿用，graded CSV]
                       sensor_pos, sensor_n, S_hall, mesh_csv_dir, n_uniform, [], AXIAL_TOL, VARIANT);  % sensor_r=[]→Ø0.3mm；axial_tol=AXIAL_TOL

%% ---- 步驟5：解 Ĥ_V = G·Vᵀ(VVᵀ)⁻¹、gauge D̄、增益 ^Bĝ_V（論文 step 9）------
%  論文 notation：G = Dv（profiled charges = D^v）；Ĥ_V = G·Vᵀ(VVᵀ)⁻¹（code 變數 H_V）；
%  ^Bĝ_V = (6/5)·h₁₁、D̄ = (5/(6·h₁₁))·Ĥ_V（h₁₁ = Ĥ_V(1,1)）。
H_V      = (Dv * Vmat.') / (Vmat * Vmat.');   % Ĥ_V = G Vᵀ(VVᵀ)⁻¹（= 舊 Dmat，單位 T/V）
D_bar    = H_V * (5 / (6 * H_V(1,1)));         % gauge：D̄(1,1)=5/6
ghat_V_B = (6/5) * H_V(1,1);                   % ^Bĝ_V：電壓側增益（T/V；H_V = ^Bĝ_V·D̄）

%% ---- 區域場誤差 + 重建檢查 -------------------------------------------------
errpct = 100 * sqrt(J / sum(Bstack(:).^2));
recon  = norm(H_V*Vmat - Dv,'fro') / norm(Dv,'fro');   % ‖Ĥ_V·V − G‖/‖G‖

%% ---- 激發欄重排成 paper P1..P6 + 重建 ê(3×6)（存檔/下游用）-----------------
[~, paper_to_apdl] = sort(apdl_to_paper_idx);   % = [1 5 2 6 4 3]
Vmat_p = Vmat(:, paper_to_apdl);                % 列=sensor P1..P6；欄=激發 P1..P6
Dv_p   = Dv(:,   paper_to_apdl);                % 列=電荷 P1..P6；欄=激發 P1..P6
E36 = zeros(3,6);                               % 17-vec ê → 3×6（含約束 e6z），欄=P1..P6
E36(:,1)=e_hat(1:3); E36(:,2)=e_hat(4:6); E36(:,3)=e_hat(7:9); E36(:,4)=e_hat(10:12); E36(:,5)=e_hat(13:15);
E36(1,6)=e_hat(16);  E36(2,6)=e_hat(17);  E36(3,6)=e_hat(1)-e_hat(4)+e_hat(8)-e_hat(11)+e_hat(15);

%% ---- 印結果（主輸出 V, D̄, ^Bĝ_V, Ĥ_V；V 先印供驗證）-----------------------
fprintf('\n========= Hall-sensor 18-param bias（variant=%s, l_hat=%.2f µm）=========\n', VARIANT, ell_hat);
fprintf('  N_p=%d | J=%.4e | region err=%.3f%% | recon ||H_V*V-G||/||G||=%.2e\n', npts, J, errpct, recon);
fprintf('  V [mV]（列=sensor P1..P6，欄=激發 P1..P6，對角=自激發）：\n');
for i=1:6, fprintf('   % .3e % .3e % .3e % .3e % .3e % .3e\n', Vmat_p(i,:)); end
fprintf('  D_bar（gauge D_bar(1,1)=5/6；列=電荷 P1..P6，欄=sensor P1..P6）：\n');
for i=1:6, fprintf('   % .4f % .4f % .4f % .4f % .4f % .4f\n', D_bar(i,:)); end
fprintf('  ghat_V_B (^B g_V) = %.4f x10^-3 mT/mV（電壓側增益，H_V = ghat_V_B·D_bar）\n', ghat_V_B*1e3);
fprintf('  H_V (Hhat_V) [mT/mV] 矩陣（列=電荷極 P1..P6，欄=sensor P1..P6；G[mT] = H_V·V[mV]）：\n');
for i=1:6, fprintf('   % .4e % .4e % .4e % .4e % .4e % .4e\n', H_V(i,:)); end
fprintf('===================================================================\n');

%% ---- 存 .mat 到本包 data/（規則#2；保留 Dv_p/ê 供下游 plot）---------------
mat_out = fullfile(mat_dir, sprintf('calib_D_%s%s.mat', VARIANT, htag));   % htag=_h1um → 不蓋 canonical
Dmat = H_V;  g_V = ghat_V_B;     % alias：維持 .mat field 名 'Dmat'/'g_V' 與下游 loader 相容（範圍：不改 .mat field 名）
save(mat_out, 'Dmat','D_bar','g_V','Dv_p','Vmat_p','exc_sign','ell_hat','e_hat','E36','Pc', ...
              'J','errpct','recon','S_hall','R_select','npts','VARIANT', ...
              'apdl_to_paper_idx','paper_to_apdl','sensor_pos','sensor_n','n_uniform');
fprintf('已存 %s\n', mat_out);

%% ---- 產生 xelatex 矩陣 PDF（V、D̄、g_V、D）-------------------------------
pole     = {'P1','P2','P3','P4','P5','P6'};
tex_path = fullfile(out_dir, sprintf('D_%s%s.tex', VARIANT, htag));
pdf_path = fullfile(out_dir, sprintf('D_%s%s.pdf', VARIANT, htag));
fid = fopen(tex_path,'w');
fprintf(fid,'%% Auto-generated by main.m — Hall-sensor 18-param bias V / Dbar / ^Bg_V / Hhat_V (P1..P6)\n');
fprintf(fid,'\\documentclass[11pt]{article}\n\\usepackage[margin=1in]{geometry}\n\\usepackage{amsmath}\n\\usepackage{amssymb}\n');
fprintf(fid,'\\begin{document}\n');
fprintf(fid,'\\begin{center}\\large Hall-sensor 18-param bias (off-axis) calibration\\quad(variant: %s)\\end{center}\n\n', ...
        strrep(VARIANT,'_','\_'));
emit_labeled_matrix(fid, 'V', Vmat_p, pole, pole, 'auto', ...
    'row $i$ = sensor pole P$i$;\ column $j$ = excited pole P$j$ (diagonal = self-excitation);\ units mV.');
emit_labeled_matrix(fid, '\bar{D}', D_bar, pole, pole, '', ...
    'row $i$ = charge pole P$i$;\ column $j$ = sensor pole P$j$;\ gauge $\bar{D}_{11}=5/6$.');
fprintf(fid,'\\[\n {}^{B}\\hat{g}_{V} = %.4f\\times10^{-3}~\\mathrm{mT/mV}\n\\]\n', ghat_V_B*1e3);
fprintf(fid,'\\noindent\\small voltage-side gain ${}^{B}\\hat{g}_{V}=\\tfrac{6}{5}\\hat{H}_{V,11}$ (mT/mV; $\\hat{H}_{V}={}^{B}\\hat{g}_{V}\\bar{D}$).\n\n');
emit_labeled_matrix(fid, '\hat{H}_{V}', H_V, pole, pole, 'auto', ...
    'row $i$ = charge pole P$i$;\ column $j$ = sensor pole P$j$;\ $\hat{H}_{V}=G\,V^{\top}(VV^{\top})^{-1}$,\ $G=\hat{H}_{V}V$ ($G=D^{v}$);\ units mT/mV.');
fprintf(fid,['\\noindent\\small 18-param bias (off-axis) model;\\ order P1--P6;\\ $I=1$\\,A (= FEM excitation);\\ ' ...
             'all-source (flip-sink: lower P1/P3/P6; $V$ same convention);\\ ' ...
             'sensor $\\varnothing0.3$\\,mm$\\times%.4g$\\,mm cylinder,\\ %d-pt interp;\\ ' ...
             '$\\hat{\\ell}=%.1f$\\,$\\mu$m,\\ region err $=%.3f\\%%$,\\ recon $=%.2e$.\n'], AXIAL_TOL*1e3, n_uniform, ell_hat, errpct, recon);
fprintf(fid,'\\end{document}\n');
fclose(fid);

xelatex = 'C:\Users\Kuo\AppData\Local\Programs\MiKTeX\miktex\bin\x64\xelatex.exe';
old = cd(out_dir);
[st, sysout] = system(sprintf('"%s" -interaction=nonstopmode -halt-on-error "%s"', xelatex, tex_path));
cd(old);
if st ~= 0 || ~exist(pdf_path,'file'); fprintf('%s\n', sysout); error('xelatex 編譯失敗（見上方輸出）。'); end
for ext = {'.tex','.aux','.log','.out'}    % results 只留 .pdf
    f = fullfile(out_dir, sprintf('D_%s%s%s', VARIANT, htag, ext{1})); if exist(f,'file'); delete(f); end
end
fprintf('已存 %s\n', pdf_path);

%% ---- local function：渲染帶 P1..P6 表頭的 6×6 矩陣（含自動 10^n 因子）----------
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
