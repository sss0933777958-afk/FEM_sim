function diag_P2P1_single(variant)
% DIAG_P2P1_SINGLE  單一 coil1（P1 激發）下 P2 sensor 的 all-source 單點內插 B·n+ + sign。
%   只需 coil1/<variant> 的 .dat（不需 6 顆全跑）。連接性恆用 baseline 標準網格 CSV
%   （gap μ_r 等效＝同網格同 node_id）。用來看「改 μ_eff 後，P2-under-P1 翻負了沒」。
%   variant 例：'standard' / 'gap200um_mueq' / 'gap20um_mueq'。
    if nargin < 1 || isempty(variant), variant = 'standard'; end
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath(['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
             'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir\code\function']);
    cnst = mt_constants();
    [sp, sn] = build_sensor_geometry(cnst);
    rr  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
    mcd = fullfile(rr,'mesh','standard','csv');

    % ---- 標準網格 → triangulation ----
    N=readmatrix(fullfile(mcd,'sensor_local_nodes.csv'));  E=readmatrix(fullfile(mcd,'sensor_local_elems.csv'));
    nid=N(:,1); P=N(:,2:4); mxid=max(nid); g2l=zeros(mxid,1); g2l(nid)=1:numel(nid);
    Sl=E(:,2:9); tets=zeros(size(Sl,1),4); kk=0;
    for r=1:size(Sl,1), u=unique(Sl(r,:),'stable'); if numel(u)==4, kk=kk+1; tets(kk,:)=g2l(u); end, end
    tets=tets(1:kk,:);
    v1=P(tets(:,2),:)-P(tets(:,1),:); v2=P(tets(:,3),:)-P(tets(:,1),:); v3=P(tets(:,4),:)-P(tets(:,1),:);
    vol=dot(v1,cross(v2,v3,2),2); bad=vol<0; tets(bad,[3 4])=tets(bad,[4 3]);
    TR=triangulation(tets,P);

    % ---- 載 coil1/variant ----
    ds=import_ansys_data(fullfile(rr,'coil1',variant),'all','coil1');
    fprintf('[%s] coil1：%d 節點，|B|max=%.4f T\n', variant, numel(ds.x), max(vecnorm([ds.bx ds.by ds.bz],2,2)));
    m2=max(max(nid),max(ds.node_id)); id2=zeros(m2,1); id2(ds.node_id)=1:numel(ds.node_id);
    li=zeros(numel(nid),1); inb=nid<=m2; li(inb)=id2(nid(inb));
    if any(li==0), error('局部節點 ID 對不上 %s .dat（網格不一致）。', variant); end
    Bnode=[ds.bx(li), ds.by(li), ds.bz(li)];

    % ---- P2 sensor 底面中心單點內插 ----
    ctr=(sp(:,2)+[0;0;cnst.SPH_OFST]).';  ni=sn(:,2);
    ti=pointLocation(TR,ctr);
    if isnan(ti), error('P2 sensor 中心不在局部網格內。'); end
    bc=cartesianToBarycentric(TR,ti,ctr); conn=TR.ConnectivityList(ti,:);
    Bp=bc*Bnode(conn,:);
    raw=Bp*ni; allsrc=-raw;          % P1 為下極 → all-source ×(-1)
    tags={'正(回流)','負(匯進P2)'}; tag=tags{(allsrc<0)+1};
    fprintf('  P2-under-P1：raw B·n+=%+.4e T ｜ all-source=%+.4e T → %s\n', raw, allsrc, tag);
end
