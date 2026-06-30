%% diag_sensor_node_spread.m -- 加密網格下，sensor 圓柱內各真實節點 B·n+ 的散布
% 問題：那 ~170 個節點投影到 n+ 的 B·n+ 都一樣嗎？
% 對每顆 sensor 取「self 激發」(自己那顆極通電) 的 coil，列圓柱內節點 B·n+ 的
% 平均 / 標準差 / min / max / 變異係數 CoV / 全距%，以及沿軸向(離極遠近)的趨勢。
% I=1A、sensor_spheres 加密網格、圓柱 Ø0.3mm×0.1mm（同 extract_Vmat）。

clear; clc;
TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
addpath(fullfile(TREE,'code','function'));
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';

SENSOR_R = 0.15e-3; AXIAL_T = 0.10e-3;
cnst = mt_constants();
apdl_to_paper_idx = [1,3,6,5,2,4];
[sensor_pos, sensor_n] = build_sensor_geometry(cnst);
plabel = {'P1','P2','P3','P4','P5','P6'};

% 每顆 sensor 的「self 激發」coil（apdl_to_paper_idx(j)==i 的 j）
self_coil = zeros(1,6);
for i=1:6, self_coil(i) = find(apdl_to_paper_idx==i); end

fprintf('加密網格 sensor 圓柱內 B·n+ 散布（self 激發、I=1A）\n');
fprintf('%-4s %5s %11s %10s %9s %9s %7s %8s\n', ...
        'P','N','mean[T]','std[T]','min[T]','max[T]','CoV%%','全距%%');
for i = 1:6
    cn = sprintf('coil%d', self_coil(i));
    da = import_ansys_data(fullfile(results_root, cn, 'sensor_spheres'),'all',cn);
    X  = [da.x, da.y, da.z - cnst.SPH_OFST];
    Bn = [da.bx, da.by, da.bz];
    ni = sensor_n(:,i);
    r  = X - sensor_pos(:,i).';
    axial = r * ni;
    rho   = vecnorm(r - axial*ni.', 2, 2);
    sel = (rho <= SENSOR_R) & (axial >= 0) & (axial <= AXIAL_T);
    bdotn = Bn(sel,:) * ni;               % 各節點 B·n+ [T]
    ax    = axial(sel);                   % 各節點離 sensor 面的軸向距離 [m]
    mu = mean(bdotn); sd = std(bdotn);
    cov = 100*sd/abs(mu); rng_pct = 100*(max(bdotn)-min(bdotn))/abs(mu);
    fprintf('%-4s %5d % .4e % .3e % .3e % .3e %6.1f %7.1f\n', ...
            plabel{i}, numel(bdotn), mu, sd, min(bdotn), max(bdotn), cov, rng_pct);
    % 沿軸趨勢：近極(axial<0.03mm) vs 遠極(axial>0.07mm) 的平均，看是否單調下降
    near = ax < 0.03e-3; far = ax > 0.07e-3;
    if any(near)&&any(far)
        drop = 100*(mean(bdotn(near))-mean(bdotn(far)))/abs(mu);
        cm = corrcoef(ax, bdotn); cc = cm(1,2);
        fprintf('        近極面均 % .4e / 遠極面均 % .4e → 跨圓柱降 %.1f%%，corr(軸距,B·n)=% .2f\n', ...
                mean(bdotn(near)), mean(bdotn(far)), drop, cc);
    end
end
fprintf('\n（CoV=std/|mean|；全距%%=(max-min)/|mean|；軸距大=離極遠）\n');
