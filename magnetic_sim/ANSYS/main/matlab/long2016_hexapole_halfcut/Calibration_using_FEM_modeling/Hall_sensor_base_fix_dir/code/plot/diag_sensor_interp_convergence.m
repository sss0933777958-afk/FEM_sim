%% diag_sensor_interp_convergence.m -- 提高取樣點數「會更無偏嗎」收斂實證
% =========================================================================
% 概念：重心內插 = 精確重建 FEM 分片線性場（節點值線性內插，無內插誤差）。
%   N 點均勻撒 = 對「圓柱體積平均」做 Monte-Carlo 估計。
%   MC 估計器對任意 N 都「無偏」（E[est]=體積平均）；提高 N 只降「變異 SE∝1/√N」，
%   不降偏差。體積平均 vs 形心值之差（場曲率）是固定性質，N 再大也不消失。
%
% 實證（每顆 self 激發 sensor）：
%   (1) ground truth 體積平均 V_GT = 巨量 N_GT 單跑（≈真值）。
%   (2) N ∈ {100,1000,10000}：跑 M 個不同 seed → 估計值的「平均」(偏差 vs V_GT)
%       與「標準差」(實測 SE)；對照理論 SE = 每點 std/√N。
%   (3) 形心單點內插 V_cen（確定值）→ (V_GT−V_cen)/V_GT = 曲率偏移（與 N 無關）。
% 參數/網格/位置/符號全對齊 main_Vmat.m / extract_Vmat_interp。純印表。
% =========================================================================
clear; clc;

%% ---- config ----------------------------------------------------------------
VARIANT  = 'gap200um_mueq';
S_hall   = 130;
SENSOR_R = 0.15e-3; AXIAL_T = 0.10e-3;
N_GT     = 5e5;                 % ground-truth 體積平均的取樣點數
N_LEVELS = [100, 1000, 10000];  % 要比較的取樣點數
M        = 30;                  % 每個 N 跑幾個 seed（估 SE）

%% ---- paths -----------------------------------------------------------------
CAL  = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
addpath(fullfile(CAL,'Hall_sensor_base_fix_dir','code','function'));
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
mesh_csv_dir = fullfile(results_root,'mesh','graded','csv');

cnst = mt_constants();
apdl_to_paper_idx = [1,3,6,5,2,4];
plabel = {'P1','P2','P3','P4','P5','P6'};
[sensor_pos, sensor_n] = build_sensor_geometry(cnst);
self_coil = zeros(1,6); for i=1:6, self_coil(i)=find(apdl_to_paper_idx==i); end

%% ---- 建 TR（同 extract_Vmat_interp）---------------------------------------
N = readmatrix(fullfile(mesh_csv_dir,'sensor_local_nodes.csv'));
E = readmatrix(fullfile(mesh_csv_dir,'sensor_local_elems.csv'));
nid = N(:,1); P = N(:,2:4); mxid=max(nid); g2l=zeros(mxid,1); g2l(nid)=1:numel(nid);
Sconn = E(:,2:9); tets=zeros(size(Sconn,1),4); kk=0;
for r=1:size(Sconn,1)
    u=unique(Sconn(r,:),'stable'); if numel(u)==4, kk=kk+1; tets(kk,:)=g2l(u); end
end
tets=tets(1:kk,:);
v1=P(tets(:,2),:)-P(tets(:,1),:); v2=P(tets(:,3),:)-P(tets(:,1),:); v3=P(tets(:,4),:)-P(tets(:,1),:);
vol=dot(v1,cross(v2,v3,2),2); bad=vol<0; tets(bad,[3 4])=tets(bad,[4 3]);
TR=triangulation(tets,P);

%% ---- 逐 self 激發 sensor 掃 N ---------------------------------------------
fprintf('提高取樣點數「會更無偏嗎」收斂實證（self 激發、%s）\n', VARIANT);
fprintf('每顆：V_GT=體積平均(N=%g)；對各 N 跑 %d seed；理論 SE=每點std/√N\n\n', N_GT, M);

for i = 1:6
    kc = self_coil(i);
    cn = sprintf('coil%d', kc);
    ds = import_ansys_data(fullfile(results_root, cn, VARIANT),'all',cn);
    m2=max(max(nid),max(ds.node_id)); id2=zeros(m2,1); id2(ds.node_id)=1:numel(ds.node_id);
    li=zeros(numel(nid),1); inb=nid<=m2; li(inb)=id2(nid(inb));
    fnode = [ds.bx(li), ds.by(li), ds.bz(li)] * sensor_n(:,i);   % 每節點 B·n+（純量場）

    ci = sensor_pos(:,i) + [0;0;cnst.SPH_OFST]; ni = sensor_n(:,i);

    % ground truth 體積平均
    rng(999);
    V_GT = S_hall * est_mean(TR, fnode, sample_cyl(ci,ni,N_GT,SENSOR_R,AXIAL_T));
    % 每點 std（預測 SE）：用 GT 那批點的散布近似
    rng(999);
    fp_gt = interp_pts(TR, fnode, sample_cyl(ci,ni,N_GT,SENSOR_R,AXIAL_T));
    perpt_std = std(fp_gt)*S_hall;
    % 形心單點
    cen = (ci + AXIAL_T/2*ni).'; tic2 = pointLocation(TR, cen);
    bcc = cartesianToBarycentric(TR, tic2, cen); cc = TR.ConnectivityList(tic2,:);
    V_cen = S_hall * sum(bcc(:).*fnode(cc(:)));

    fprintf('--- sensor %s（coil%d 自激）：V_GT=% .5e V，形心偏移=%+.2f%% ---\n', ...
            plabel{i}, kc, V_GT, 100*(V_GT-V_cen)/abs(V_GT));
    fprintf('  %7s %12s %12s %12s %12s\n','N','估計均[V]','偏差%vsGT','實測SE%','理論SE%(std/√N)');
    for nl = 1:numel(N_LEVELS)
        Nn = N_LEVELS(nl); ests = zeros(M,1);
        for s = 1:M
            rng(s);   % 不同 seed
            ests(s) = S_hall * est_mean(TR, fnode, sample_cyl(ci,ni,Nn,SENSOR_R,AXIAL_T));
        end
        bias = 100*(mean(ests)-V_GT)/abs(V_GT);
        se   = 100*std(ests)/abs(V_GT);
        se_t = 100*(perpt_std/sqrt(Nn))/abs(V_GT);
        fprintf('  %7d % .5e %+11.2f %11.2f %11.2f\n', Nn, mean(ests), bias, se, se_t);
    end
    fprintf('\n');
end
fprintf(['判讀：偏差%% 對所有 N 都≈0（在 SE/√%d 內）⇒ 估計器本就無偏，提高 N「不會更無偏」；\n' ...
         '      實測 SE%% 隨 N 每 ×10 約 /3.16（∝1/√N）⇒ 提高 N 降的是「雜訊/變異」，使估計更穩；\n' ...
         '      形心偏移%% 與 N 無關（場曲率固定）⇒ 那是體積平均 vs 點值的本質差，撒再多點也不消失。\n'], M);

%% ---- locals ----------------------------------------------------------------
function pts = sample_cyl(ci, ni, Nn, r0, h0)
    t1=[-ni(2);ni(1);0]; if norm(t1)<1e-9, t1=[1;0;0]; end
    t1=t1/norm(t1); t2=cross(ni,t1);
    a=h0*rand(Nn,1); r=r0*sqrt(rand(Nn,1)); th=2*pi*rand(Nn,1);
    pts = ci.' + a.*ni.' + (r.*cos(th)).*t1.' + (r.*sin(th)).*t2.';
end
function fp = interp_pts(TR, fnode, pts)
    ti=pointLocation(TR,pts); good=~isnan(ti);
    bc=cartesianToBarycentric(TR,ti(good),pts(good,:));
    conn=TR.ConnectivityList(ti(good),:);
    fp = sum(bc.*fnode(conn),2);
end
function v = est_mean(TR, fnode, pts)
    v = mean(interp_pts(TR, fnode, pts));
end
