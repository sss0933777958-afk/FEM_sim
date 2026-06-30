%% compare_sensor_methods.m -- Task1+4：完美接觸下三種 sensor 讀值取法的比較
% =========================================================================
%  使用者要求（完美接觸 = 幾何完全接觸、μ_steel=280 的 baseline）：
%   比較三種「sensor 電壓矩陣 Vmat」取法，以「單點內插」為參考算誤差百分比：
%     M1 單點內插    extract_Vmat_interp_center  (standard 粗網格，sensor 底面中心單點重心內插) ← 參考
%     M2 1000 點內插 extract_Vmat_interp         (standard 粗網格，圓柱內均勻撒 1000 點重心內插平均)
%     M3 加密網格    extract_Vmat                (sensor_spheres 加密網格，圓柱內真實節點平均)
%   三者共用同一組 sensor 幾何(build_sensor_geometry) + 同 S_hall + 同 all-source 翻號。
%
%  Task4：baseline(standard 粗網格)下，單一 sensor 的取樣範圍(Ø0.3mm×0.1mm 圓柱)是否落在
%   「同一個 tet」內？→ 把 M2 的 1000 取樣點丟 pointLocation 數不重複 tet id。
%   若每顆都唯一 → 解釋為何 M1≈M2（tet 內 FEM 場線性，體積平均≈形心值≈中心點值）。
%
%  本檔只算 + 印 + 寫 results txt（圖另外定案後再存）。電流 = 1A(FEM 激發)。
% =========================================================================

clear; clc;

%% ---- paths ----------------------------------------------------------------
TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % mt_constants/import_ansys_data
addpath(fullfile(TREE,'code','function'));                                                      % 三個 extract_Vmat*
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
mesh_csv_dir = fullfile(results_root,'mesh','standard','csv');   % sensor_local_{nodes,elems}.csv（standard 拓樸，csv/ 子夾）
resdir       = fullfile(TREE,'results');
if ~exist(resdir,'dir'); mkdir(resdir); end

%% ---- config ---------------------------------------------------------------
S_hall   = 130;            % Hall 靈敏度 [V/T]（EQ-730L）
N_UNIF   = 1000;           % M2 每 sensor 取樣點數
SENSOR_R = 0.15e-3;        % sensor 圓柱半徑 [m]
AXIAL_T  = 0.10e-3;        % sensor 圓柱高 [m]
plabel   = {'P1','P2','P3','P4','P5','P6'};

%% ---- 共用：常數 + sensor 幾何 + APDL→paper 對映 --------------------------
cnst = mt_constants();
apdl_to_paper_idx = [1,3,6,5,2,4];                       % APDL coil j → paper pole（判下極翻號用）
[sensor_pos, sensor_n] = build_sensor_geometry(cnst);     % 6 顆 sensor 中心 + 法線 n+（兩法 M1/M2 共用）

%% ===== M1 單點內插（參考）==================================================
fprintf('===== M1：單點內插（standard，sensor 底面中心）=====\n');
[V1, exc_sign] = extract_Vmat_interp_center(results_root, cnst, apdl_to_paper_idx, ...
                                            sensor_pos, sensor_n, S_hall, mesh_csv_dir, 'standard');

%% ===== M2 1000 點內插 ======================================================
fprintf('\n===== M2：1000 點內插（standard，圓柱均勻取樣平均）=====\n');
[V2, ~] = extract_Vmat_interp(results_root, cnst, apdl_to_paper_idx, ...
                              sensor_pos, sensor_n, S_hall, mesh_csv_dir, N_UNIF, SENSOR_R, AXIAL_T);

%% ===== M3 加密網格 real-node ===============================================
fprintf('\n===== M3：加密網格（sensor_spheres，真實節點平均）=====\n');
[V3, ~] = extract_Vmat(results_root, cnst, apdl_to_paper_idx, ...
                       sensor_pos, sensor_n, S_hall, 'sensor_spheres', SENSOR_R, AXIAL_T);

%% ---- 誤差百分比（以 M1 為參考）-------------------------------------------
% 元素誤差 %（小分母會放大，故另給 Frobenius 相對誤差當穩健單一指標）
E2 = (V2 - V1)./V1 * 100;
E3 = (V3 - V1)./V1 * 100;
% 物理上的「self 讀值」不在矩陣對角：欄 j(APDL coil) 激發 paper pole apdl_to_paper_idx(j)，
% 故 self = (row i == apdl_to_paper_idx(j))。其餘為 cross。
selfidx = false(6,6);
for j = 1:6, selfidx(apdl_to_paper_idx(j), j) = true; end
% cross 再分「強 cross」(|V1| 夠大、% 有意義) 與「近零 cross」(對極 opposing pair，|V1|≈0、% 失真)
tolV  = 0.10 * max(abs(V1(selfidx)));        % 門檻 = self 最大值的 10%（分開「強鄰極耦合」與「近零對極耦合」）
small = abs(V1) < tolV;                       % 近零元素：% 無意義
crossStrong = ~selfidx & ~small;              % 強 cross
crossNear   = ~selfidx &  small;              % 近零 cross（opposing pair）

froErr = @(Vk) 100*norm(Vk-V1,'fro')/norm(V1,'fro');
mae    = @(Em,idx) mean(abs(Em(idx)));
mxe    = @(Em,idx) max(abs(Em(idx)));

%% ===== Task4：baseline 取樣範圍是否落在單一 tet ============================
fprintf('\n===== Task4：standard 網格下 sensor 取樣範圍的 tet 數 =====\n');
% 重建 standard 局部 triangulation（同 extract_Vmat_interp 的 tet 化）
N = readmatrix(fullfile(mesh_csv_dir,'sensor_local_nodes.csv'));   % [id x y z]
E = readmatrix(fullfile(mesh_csv_dir,'sensor_local_elems.csv'));   % [eid s1..s8]
nid = N(:,1); Pn = N(:,2:4); mxid = max(nid); g2l = zeros(mxid,1); g2l(nid) = 1:numel(nid);
Sl = E(:,2:9); tets = zeros(size(Sl,1),4); kk = 0;
for r = 1:size(Sl,1)
    u = unique(Sl(r,:),'stable');
    if numel(u) == 4, kk = kk+1; tets(kk,:) = g2l(u); end
end
tets = tets(1:kk,:);
v1 = Pn(tets(:,2),:)-Pn(tets(:,1),:); v2 = Pn(tets(:,3),:)-Pn(tets(:,1),:); v3 = Pn(tets(:,4),:)-Pn(tets(:,1),:);
vol = dot(v1,cross(v2,v3,2),2); bad = vol < 0; tets(bad,[3 4]) = tets(bad,[4 3]);
TR = triangulation(tets, Pn);

% 重現 M2 的取樣點（rng(0) + 同圓柱幾何、同 rand 呼叫順序 → 與 extract_Vmat_interp 一致）
rng(0); samp = cell(1,6);
for i = 1:6
    ci = sensor_pos(:,i) + [0;0;cnst.SPH_OFST]; ni = sensor_n(:,i);
    t1 = [-ni(2); ni(1); 0]; if norm(t1) < 1e-9, t1 = [1;0;0]; end
    t1 = t1/norm(t1); t2 = cross(ni, t1);
    a  = AXIAL_T  * rand(N_UNIF,1);
    rr = SENSOR_R * sqrt(rand(N_UNIF,1));
    th = 2*pi     * rand(N_UNIF,1);
    samp{i} = ci.' + a.*ni.' + (rr.*cos(th)).*t1.' + (rr.*sin(th)).*t2.';
end

n_tet   = zeros(1,6);   % 每 sensor 取樣點覆蓋的不重複 tet 數
n_in    = zeros(1,6);   % 落在網格內的取樣點數
ctr_tet = zeros(1,6);   % sensor 中心點(M1 評估點)所在 tet
same_as_ctr = false(1,6);
for i = 1:6
    ti = pointLocation(TR, samp{i}); good = ~isnan(ti);
    n_in(i)  = nnz(good);
    utet     = unique(ti(good));
    n_tet(i) = numel(utet);
    pc       = (sensor_pos(:,i) + [0;0;cnst.SPH_OFST]).';
    ctr_tet(i) = pointLocation(TR, pc);
    same_as_ctr(i) = ismember(ctr_tet(i), utet);
    fprintf('  sensor %s: 取樣點 %d/%d 在網格內、覆蓋 %d 個 tet（中心 tet=%d，%s）\n', ...
            plabel{i}, n_in(i), N_UNIF, n_tet(i), ctr_tet(i), ...
            ternary(n_tet(i)==1,'單一 tet','跨多 tet'));
end
all_single = all(n_tet==1);

%% ---- console 摘要 ---------------------------------------------------------
fprintf('\n========================= 摘要 =========================\n');
fprintf('Frobenius 相對誤差(對 M1)：       M2 = %.2f %%   M3 = %.2f %%\n', froErr(V2), froErr(V3));
fprintf('self 讀值 |err%%| mean / max ：    M2 = %.2f / %.2f %%   M3 = %.2f / %.2f %%\n', ...
        mae(E2,selfidx), mxe(E2,selfidx), mae(E3,selfidx), mxe(E3,selfidx));
fprintf('強 cross |err%%| mean / max ：     M2 = %.2f / %.2f %%   M3 = %.2f / %.2f %%\n', ...
        mae(E2,crossStrong), mxe(E2,crossStrong), mae(E3,crossStrong), mxe(E3,crossStrong));
fprintf('近零 cross(opposing) |err%%| max ：M2 = %.1f %%   M3 = %.1f %%  (|V1|<self 10%%，%% 失真、僅供參考)\n', ...
        mxe(E2,crossNear), mxe(E3,crossNear));
fprintf('Task4：每顆 sensor 取樣範圍是否單一 tet？ → %s\n', ternary(all_single,'是（全 6 顆單一 tet）','否（跨多 tet）'));
fprintf('=========================================================\n');

%% ---- 寫 results txt -------------------------------------------------------
txt = fullfile(resdir,'sensor_methods_compare.txt');
fid = fopen(txt,'w');
fprintf(fid,'Sensor 讀值取法比較（完美接觸 baseline，I=1A，S_hall=%g V/T）\n', S_hall);
fprintf(fid,'參考 = M1 單點內插(extract_Vmat_interp_center, standard)\n');
fprintf(fid,'M2 = 1000 點內插(extract_Vmat_interp, standard, n=%d)\n', N_UNIF);
fprintf(fid,'M3 = 加密網格 real-node(extract_Vmat, sensor_spheres)\n');
fprintf(fid,'圓柱：半徑 %.3f mm、高 %.3f mm；all-source 翻下極(P1/P3/P6)。\n\n', SENSOR_R*1e3, AXIAL_T*1e3);

writemat(fid,'V1 = M1 單點內插 [V]（列=sensor P1..P6，欄=激發 coil1..6）', V1, '% .4e');
writemat(fid,'V2 = M2 1000 點內插 [V]', V2, '% .4e');
writemat(fid,'V3 = M3 加密網格 [V]', V3, '% .4e');
writemat(fid,'E2 = (V2-V1)/V1 ×100 [%]', E2, '% .2f');
writemat(fid,'E3 = (V3-V1)/V1 ×100 [%]', E3, '% .2f');

fprintf(fid,'--- 摘要（self = 被激發極自己的 sensor 讀值，不在矩陣對角；欄 j 激發 paper P%s）---\n', '');
fprintf(fid,'self 讀值位置：');
for j=1:6, fprintf(fid,'(row %d,col %d) ', apdl_to_paper_idx(j), j); end
fprintf(fid,'\n');
fprintf(fid,'Frobenius 相對誤差(對 M1)：M2 = %.2f %% , M3 = %.2f %%\n', froErr(V2), froErr(V3));
fprintf(fid,'self 讀值 |err%%| mean/max：M2 = %.2f/%.2f %% , M3 = %.2f/%.2f %%\n', ...
        mae(E2,selfidx),mxe(E2,selfidx),mae(E3,selfidx),mxe(E3,selfidx));
fprintf(fid,'強 cross |err%%| mean/max：M2 = %.2f/%.2f %% , M3 = %.2f/%.2f %%\n', ...
        mae(E2,crossStrong),mxe(E2,crossStrong),mae(E3,crossStrong),mxe(E3,crossStrong));
fprintf(fid,'近零 cross(opposing pair) |err%%| max：M2 = %.1f %% , M3 = %.1f %%（|V1|<self 10%%，%% 失真僅參考）\n\n', ...
        mxe(E2,crossNear),mxe(E3,crossNear));

fprintf(fid,'--- Task4：standard 網格 sensor 取樣範圍 tet 數 ---\n');
for i = 1:6
    fprintf(fid,'sensor %s：取樣點 %d/%d 在網格內、覆蓋 %d 個 tet（中心 tet=%d，%s）\n', ...
            plabel{i}, n_in(i), N_UNIF, n_tet(i), ctr_tet(i), ternary(n_tet(i)==1,'單一 tet','跨多 tet'));
end
fprintf(fid,'結論：%s\n', ternary(all_single, ...
    '6 顆 sensor 的取樣範圍都落在單一 tet → M1(中心點) ≈ M2(體積平均)，因 tet 內 FEM 場線性。', ...
    '部分 sensor 取樣跨多 tet → M1 與 M2 會有差異。'));
fclose(fid);
fprintf('已寫 %s\n', txt);

% 暫存供畫圖定案用（不存最終圖；待與使用者確認 style 後再出 PNG）
save(fullfile(resdir,'sensor_methods_compare_data.mat'), ...
     'V1','V2','V3','E2','E3','n_tet','n_in','ctr_tet','same_as_ctr','plabel','S_hall','N_UNIF');

%% ---- local helpers --------------------------------------------------------
function s = ternary(c,a,b); if c, s=a; else, s=b; end; end
function writemat(fid,ttl,M,fmt)
    fprintf(fid,'%s：\n',ttl);
    for i=1:size(M,1)
        for j=1:size(M,2); fprintf(fid,[fmt ' '],M(i,j)); end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end
