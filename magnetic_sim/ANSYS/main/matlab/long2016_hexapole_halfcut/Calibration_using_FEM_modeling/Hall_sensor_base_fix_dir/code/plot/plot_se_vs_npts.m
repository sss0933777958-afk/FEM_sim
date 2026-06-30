%% plot_se_vs_npts.m -- 6 顆 sensor 的 MC 標準誤差 SE = σ/√N vs 內插取樣點數
% =========================================================================
% 橫軸：圓柱內均勻撒點數 N（0~60000）；縱軸：SE = σ_i/√N [V]。
% σ_i = 各 sensor 自激發下，圓柱內內插 B·n+ 的每點標準差 × S_hall（巨量樣本估）。
% 6 顆 sensor 各一條不同顏色線疊在一起。風格 ①粗體框圖。
% 純理論曲線（σ/√N），對齊 main_Vmat / extract_Vmat_interp 的網格/位置/符號。
% 輸出實檔 → figures/se_vs_npts.png（覆蓋迭代）。
% =========================================================================
clear; clc;

%% ---- config ----------------------------------------------------------------
VARIANT  = 'gap200um_mueq';
S_hall   = 130;            % [V/T]
SENSOR_R = 0.15e-3; AXIAL_T = 0.10e-3;
N_SIGMA  = 2e5;            % 估 σ 用的巨量取樣點數
N_MAX    = 60000;         % 橫軸上限

%% ---- paths -----------------------------------------------------------------
CAL  = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
addpath(fullfile(CAL,'Hall_sensor_base_fix_dir','code','function'));
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
mesh_csv_dir = fullfile(results_root,'mesh','graded','csv');
TREE = fullfile(CAL,'Hall_sensor_base_fix_dir');
fig_dir = fullfile(TREE,'figures'); if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

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

%% ---- 每顆 sensor 的 σ（自激發、巨量樣本）----------------------------------
sigma = zeros(1,6);
for i = 1:6
    kc = self_coil(i); cn = sprintf('coil%d', kc);
    ds = import_ansys_data(fullfile(results_root, cn, VARIANT),'all',cn);
    m2=max(max(nid),max(ds.node_id)); id2=zeros(m2,1); id2(ds.node_id)=1:numel(ds.node_id);
    li=zeros(numel(nid),1); inb=nid<=m2; li(inb)=id2(nid(inb));
    fnode = [ds.bx(li), ds.by(li), ds.bz(li)] * sensor_n(:,i);   % 每節點 B·n+
    ci = sensor_pos(:,i) + [0;0;cnst.SPH_OFST]; ni = sensor_n(:,i);
    rng(7);
    pts = sample_cyl(ci,ni,N_SIGMA,SENSOR_R,AXIAL_T);
    ti=pointLocation(TR,pts); good=~isnan(ti);
    bc=cartesianToBarycentric(TR,ti(good),pts(good,:));
    conn=TR.ConnectivityList(ti(good),:);
    fp = sum(bc.*fnode(conn),2);                                 % 每點內插 B·n+
    sigma(i) = std(fp) * S_hall;                                 % σ_i [V]
    fprintf('  %s: σ=% .4e V（%d/%d 命中）\n', plabel{i}, sigma(i), nnz(good), N_SIGMA);
end

%% ---- SE = σ/√N 曲線 --------------------------------------------------------
Nv = unique(round(logspace(0, log10(N_MAX), 1000)));   % log-spaced → 線性軸上平滑渲染近 0 的陡降
SE = sigma(:) ./ sqrt(Nv);                              % 6×numel(Nv)，列=sensor

%% ---- 畫圖（風格 ①粗體框圖）------------------------------------------------
cols = [0.00 0.20 0.60;    % P1 深藍
        0.85 0.10 0.10;    % P2 紅
        0.00 0.55 0.20;    % P3 綠
        0.95 0.55 0.00;    % P4 橘
        0.50 0.10 0.65;    % P5 紫
        0.00 0.60 0.72];   % P6 青
fig = figure('Color','w','Units','inches','Position',[1 1 7.8 5.2]);
ax = axes(fig); hold(ax,'on');
h = gobjects(1,6);
for i = 1:6
    h(i) = plot(ax, Nv, SE(i,:), '-', 'Color', cols(i,:), 'LineWidth', 2.5);
end
hold(ax,'off');
xlim(ax,[0 N_MAX]);
set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
box(ax,'on'); grid(ax,'off');
xt=get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));
yt=get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
xlabel(ax,'interpolation points  N','FontWeight','bold');
ylabel(ax,'SE = \sigma/\surdN  (V)','FontWeight','bold');
lg = legend(h, plabel, 'Location','northeast'); lg.FontSize=14; lg.FontWeight='bold';

out = fullfile(fig_dir,'se_vs_npts.png');
exportgraphics(fig, out, 'Resolution', 150);
fprintf('已輸出 %s\n', out);

%% ---- local -----------------------------------------------------------------
function pts = sample_cyl(ci, ni, Nn, r0, h0)
    t1=[-ni(2);ni(1);0]; if norm(t1)<1e-9, t1=[1;0;0]; end
    t1=t1/norm(t1); t2=cross(ni,t1);
    a=h0*rand(Nn,1); r=r0*sqrt(rand(Nn,1)); th=2*pi*rand(Nn,1);
    pts = ci.' + a.*ni.' + (r.*cos(th)).*t1.' + (r.*sin(th)).*t2.';
end
