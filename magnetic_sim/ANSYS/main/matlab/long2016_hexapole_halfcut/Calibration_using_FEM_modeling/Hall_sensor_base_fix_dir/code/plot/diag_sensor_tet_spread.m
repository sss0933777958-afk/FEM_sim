%% diag_sensor_tet_spread.m -- sensor 圓柱「跨幾顆 tet + tet 四角值散布」內插失真診斷
% =========================================================================
% 回答兩個問題（針對 extract_Vmat_interp 的 100 點圓柱重心內插這條路）：
%   1. 每顆 sensor 的取樣圓柱是「包在一顆 tet 內」還是「跨幾顆 tet」？
%   2. 落入的 tet，四個角節點的 B·n+ 差異大不大？→ 100 點平均後有沒有失真？
%
% 與實際 Vmat 抽取「逐點一致」：內聯複製 extract_Vmat_interp 的「建 TR」段與
%   「圓柱撒點」段（同 rng(0)、同 sensor_r/axial_tol/n_uniform、同修正 sensor 位置）。
% 參數對齊 main_Vmat.m：VARIANT='gap200um_mueq'、S_hall=130、n_uniform=100、graded csv、
%   all-source（翻下極 P1/P3/P6）。純 console 印表，不出圖 / 不寫 .mat。
%
% A. 幾何（與 coil 無關，6 顆各一次）：命中數/100、相異 tet 數 K、單 tet 最大佔比、
%    鄰近 tet 邊長 vs 圓柱尺寸（Ø0.30mm × 高 0.10mm）。
% B. tet 四角散布（self 激發）：每顆被命中 tet 4 角 B·n+ 的全距%/CoV%，跨 tet median/max。
% C. 失真量化：100 點內插平均 vs (a)最近角節點平均 (b)形心單點內插 的相對差%。
% =========================================================================
clear; clc;

%% ---- config（對齊 main_Vmat.m）------------------------------------------
VARIANT   = 'gap200um_mueq';
S_hall    = 130;            % [V/T]
n_uniform = 100;           % 圓柱內均勻取樣點數
SENSOR_R  = 0.15e-3;       % 圓柱半徑 [m]（同 extract_Vmat_interp 預設）
AXIAL_T   = 0.10e-3;       % 圓柱高 [m]
RUN_ALL_COILS = false;     % true = B/C 對 6 coil 全掃（幾何 A 不變）；預設只 self 激發

%% ---- paths（沿用，不另寫）-----------------------------------------------
CAL  = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % mt_constants/import_ansys_data
addpath(fullfile(CAL,'Hall_sensor_base_fix_dir','code','function'));                            % build_sensor_geometry
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
mesh_csv_dir = fullfile(results_root,'mesh','graded','csv');

%% ---- 常數 + 慣例 + 修正 sensor 位置 -------------------------------------
cnst = mt_constants();
apdl_to_paper_idx = [1,3,6,5,2,4];
plabel = {'P1','P2','P3','P4','P5','P6'};
[sensor_pos, sensor_n] = build_sensor_geometry(cnst);   % 3×6 WP 框 [m]（下極 −β 修正版）

% 每顆 sensor(paper Pi) 的「self 激發」coil（apdl index）
self_coil = zeros(1,6);
for i = 1:6, self_coil(i) = find(apdl_to_paper_idx == i); end

%% ---- 建 TR（內聯複製 extract_Vmat_interp）-------------------------------
N = readmatrix(fullfile(mesh_csv_dir,'sensor_local_nodes.csv'));   % [id x y z]
E = readmatrix(fullfile(mesh_csv_dir,'sensor_local_elems.csv'));   % [eid s1..s8]
nid = N(:,1); P = N(:,2:4); mxid = max(nid); g2l = zeros(mxid,1); g2l(nid) = 1:numel(nid);
Sconn = E(:,2:9); tets = zeros(size(Sconn,1),4); kk = 0;
for r = 1:size(Sconn,1)
    u = unique(Sconn(r,:),'stable');                              % SOLID96 8 槽 → tet 4 相異
    if numel(u) == 4, kk = kk+1; tets(kk,:) = g2l(u); end
end
tets = tets(1:kk,:);
v1 = P(tets(:,2),:)-P(tets(:,1),:); v2 = P(tets(:,3),:)-P(tets(:,1),:); v3 = P(tets(:,4),:)-P(tets(:,1),:);
vol = dot(v1,cross(v2,v3,2),2); bad = vol < 0; tets(bad,[3 4]) = tets(bad,[4 3]);
TR = triangulation(tets, P);
fprintf('graded sensor-local 網格：%d 節點 / %d tet\n', numel(nid), size(tets,1));

%% ---- 圓柱撒點（內聯複製 extract_Vmat_interp，同 rng(0)）-----------------
rng(0);
samp = cell(1,6);
for i = 1:6
    ci = sensor_pos(:,i) + [0;0;cnst.SPH_OFST];     % sensor 中心(ANSYS 框)
    ni = sensor_n(:,i);
    t1 = [-ni(2); ni(1); 0]; if norm(t1) < 1e-9, t1 = [1;0;0]; end
    t1 = t1/norm(t1); t2 = cross(ni, t1);
    a  = AXIAL_T  * rand(n_uniform,1);
    r  = SENSOR_R * sqrt(rand(n_uniform,1));
    th = 2*pi     * rand(n_uniform,1);
    samp{i} = ci.' + a.*ni.' + (r.*cos(th)).*t1.' + (r.*sin(th)).*t2.';   % n_uniform×3
end

%% ---- A. 幾何：跨幾顆 tet（與 coil 無關）---------------------------------
% 預先算每顆 tet 的「特徵邊長」（6 邊均值），給「圓柱 ⊂ 1 tet?」做尺度對照
edge_pairs = [1 2;1 3;1 4;2 3;2 4;3 4];
tet_edge_mean = zeros(size(tets,1),1);
for t = 1:size(tets,1)
    pp = P(tets(t,:),:);
    el = zeros(6,1);
    for e = 1:6, el(e) = norm(pp(edge_pairs(e,1),:) - pp(edge_pairs(e,2),:)); end
    tet_edge_mean(t) = mean(el);
end

ti_all   = cell(1,6);   % 每顆 sensor 各取樣點所屬 tet index（NaN=未命中）
good_all = cell(1,6);
fprintf('\n=============== A. 幾何：sensor 圓柱跨幾顆 tet（與 coil 無關）===============\n');
fprintf('圓柱：Ø%.2f mm × 高 %.2f mm；n=%d 點\n', 2*SENSOR_R*1e3, AXIAL_T*1e3, n_uniform);
fprintf('%-4s %8s %7s %10s %12s %12s\n','P','命中/n','K tet','最大佔比','tet邊長均[mm]','vs 圓柱直徑');
for i = 1:6
    pts = samp{i};
    ti  = pointLocation(TR, pts); good = ~isnan(ti);
    ti_all{i} = ti; good_all{i} = good;
    tih = ti(good);
    uK  = unique(tih);
    K   = numel(uK);
    % 每顆被命中 tet 的點數 → 最大單 tet 佔命中點比例
    if ~isempty(tih)
        cnts = histc(tih, uK); maxshare = max(cnts)/numel(tih);
        emean = mean(tet_edge_mean(uK));                  % 被命中 tet 的平均邊長
    else
        maxshare = NaN; emean = NaN;
    end
    fprintf('%-4s %5d/%-3d %7d %9.0f%% %11.4f %11.2fx\n', ...
            plabel{i}, nnz(good), n_uniform, K, 100*maxshare, emean*1e3, emean/(2*SENSOR_R));
end
fprintf('（K tet=圓柱跨的相異四面體數，K=1 即「包在一顆 tet 內」；最大佔比=最多點落同一顆 tet 的比例）\n');

%% ---- B+C：四角散布 + 失真（逐 coil；預設只 self 激發）-------------------
if RUN_ALL_COILS, coil_list = 1:6; else, coil_list = unique(self_coil); end

for kc_idx = 1:numel(coil_list)
    kc = coil_list(kc_idx);
    cn = sprintf('coil%d', kc);
    ds = import_ansys_data(fullfile(results_root, cn, VARIANT),'all',cn);
    m2 = max(max(nid), max(ds.node_id)); id2 = zeros(m2,1); id2(ds.node_id) = 1:numel(ds.node_id);
    li = zeros(numel(nid),1); inb = nid <= m2; li(inb) = id2(nid(inb));
    if any(li == 0), error('局部節點 ID 對不上 %s 的 .dat（網格不一致）。', VARIANT); end
    Bnode = [ds.bx(li), ds.by(li), ds.bz(li)];          % 每局部節點 B（對齊局部節點順序）

    % 哪些 sensor 要在這顆 coil 報？（self 模式：只報 self_coil==kc 的那顆）
    if RUN_ALL_COILS, sensors_here = 1:6; else, sensors_here = find(self_coil == kc); end

    fprintf('\n=============== B+C. coil%d（paper %s 激發）：四角散布 + 失真 ===============\n', ...
            kc, plabel{apdl_to_paper_idx(kc)});
    fprintf('%-4s %11s | 四角全距%% (med/max) | 四角CoV%% (med/max) | %11s %11s\n', ...
            'P(sensor)','interp[V]','失真%vs最近','vs形心');
    for ii = 1:numel(sensors_here)
        i  = sensors_here(ii);
        pts = samp{i}; ni = sensor_n(:,i);
        ti = ti_all{i}; good = good_all{i};
        conn = TR.ConnectivityList(ti(good),:);          % nGood×4 角節點(局部 index)
        bc   = cartesianToBarycentric(TR, ti(good), pts(good,:));   % nGood×4 重心權重

        % --- 內插 B·n+（= 實際 Vmat 抽法）---
        Bp = zeros(nnz(good),3);
        for c = 1:4, Bp = Bp + bc(:,c).*Bnode(conn(:,c),:); end
        bn_interp = Bp * ni;                             % 各點內插 B·n+
        V_interp  = S_hall * mean(bn_interp);            % = Vmat(i,kc) 的量（未套 all-source 符號）

        % --- B. 每顆被命中 tet 的 4 角 B·n+ 散布 ---
        corner_bn = (Bnode * ni);                        % 每局部節點的 B·n+
        uK = unique(ti(good));
        rng_pct = zeros(numel(uK),1); cov_pct = zeros(numel(uK),1);
        for t = 1:numel(uK)
            vals = corner_bn(tets(uK(t),:));             % 該 tet 4 角 B·n+
            mu = mean(vals);
            rng_pct(t) = 100*(max(vals)-min(vals))/abs(mu);
            cov_pct(t) = 100*std(vals)/abs(mu);
        end

        % --- C. 失真：interp 平均 vs (a)最近角節點 (b)形心單點內插 ---
        % (a) 每點取重心權重最大的角 → 該角 B·n+，再平均
        [~, amax] = max(bc, [], 2);
        lin = sub2ind(size(conn), (1:size(conn,1)).', amax);
        bn_near = corner_bn(conn(lin));
        V_near  = S_hall * mean(bn_near);
        % (b) 圓柱形心單點內插
        ci  = sensor_pos(:,i) + [0;0;cnst.SPH_OFST];
        cen = (ci + AXIAL_T/2*ni).';                     % 形心 = 中心 + 半高沿 n+
        tic2 = pointLocation(TR, cen);
        if ~isnan(tic2)
            bcc = cartesianToBarycentric(TR, tic2, cen);
            cc  = TR.ConnectivityList(tic2,:);
            Bc  = zeros(1,3); for c=1:4, Bc = Bc + bcc(c)*Bnode(cc(c),:); end
            V_cen = S_hall * (Bc*ni);
        else
            V_cen = NaN;
        end
        d_near = 100*(V_near - V_interp)/abs(V_interp);
        d_cen  = 100*(V_cen  - V_interp)/abs(V_interp);

        fprintf('%-4s % .4e | %7.1f / %7.1f | %7.1f / %7.1f | %+9.1f%% %+9.1f%%\n', ...
                plabel{i}, V_interp, median(rng_pct), max(rng_pct), ...
                median(cov_pct), max(cov_pct), d_near, d_cen);
    end
end
fprintf(['\n判讀：四角全距%%/CoV%% 小 ⇒ tet 內場近均勻；失真%% 小 ⇒ 內插平均≈真值、不失真。\n' ...
         '（全距%%=(max-min)/|mean|；CoV%%=std/|mean|；失真%%>~5%% 建議改用加密 sensor_spheres 重抽。）\n']);
