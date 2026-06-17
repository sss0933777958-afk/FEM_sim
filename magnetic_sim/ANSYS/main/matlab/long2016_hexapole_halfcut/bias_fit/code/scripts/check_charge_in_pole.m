%% check_charge_in_pole.m
%  ============================================================================
%  檢查 no_fix_l(含 bias)fit 出的 6 顆等效電荷,掃 R=50:5:500,是否都落在磁極鐵內。
%  電荷位置回報在 actuator 座標系(使用者要的);in-iron 判定用既有 filter_iron_nodes
%  (錐幾何:沿 pole_axis 投影 + 該處錐半徑;safety_r=0 純幾何)。只讀,不寫檔。
%  ============================================================================
clear; clc;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\hexapole-long2016\analysis');
model = 'long2016_hexapole_halfcut';

cnst = mt_constants();
S = load(fullfile(matlab_path(model,'charge_fit','fitting_trend'),'sweep_nofixl_vs_R.mat'));
R_um = S.R_um;  Rrot = S.R_act;  Pc_base = S.Pc_base;
plabel = {'P1','P2','P3','P4','P5','P6'};
nR = numel(R_um);

in_iron = false(6, nR);                                      % 每極每 R 是否在鐵內
for ri = 1:nR
    E = reshape(S.Esave(:,ri), 3, 6);                       % 3x6 bias(含 e6z)
    charge_act = S.ell_R(ri) * (Pc_base + E);               % actuator 框電荷位置 [m]
    charge_meas = Rrot.' * charge_act;                      % 旋回 measure/WP 框 [m]
    % 用 filter_iron_nodes 判 in-iron(傳 APDL 座標:z_apdl = z_wp + SPH_OFST)
    air = filter_iron_nodes(charge_meas(1,:).', charge_meas(2,:).', ...
                            (charge_meas(3,:)+cnst.SPH_OFST).', cnst, ...
                            struct('visualize',false,'safety_r',0));
    in_iron(:,ri) = ~air(:);
end

% ---- 每極對「自己錐」的離面餘量 margin = r_perp - r_cone(>0 = 在錐外多少 um)----
margin = nan(6, nR);                                        % um
for ri = 1:nR
    E = reshape(S.Esave(:,ri),3,6); charge_act = S.ell_R(ri)*(Pc_base+E);
    charge_meas = Rrot.' * charge_act;
    for k = 1:6
        tip = [cnst.pole_tip_x(k); cnst.pole_tip_y(k); cnst.pole_tip_z_wp(k)];
        v = charge_meas(:,k) - tip;  ax = cnst.pole_axis(:,k);
        s = v.'*ax;  rperp = sqrt(max(0, v.'*v - s^2));
        rcone = cnst.POLE_TIP_R + s*(cnst.POLE_R-cnst.POLE_TIP_R)/cnst.POLE_CONE_LEN;
        margin(k,ri) = (rperp - rcone)*1e6;                % um;>0 在錐外
    end
end
fprintf('\n--- P3 對自己錐的離面餘量(um;>0=錐外多少)在掉出去的 R ---\n');
outidx = find(margin(3,:) > 0);
for ri = outidx, fprintf('  R=%3d: +%.1f um 在錐外\n', R_um(ri), margin(3,ri)); end
fprintf('  P3 最大超出 = %.1f um\n', max(margin(3,:)));
fprintf('--- 下極在「最深在內」時離錐面還有多少(min margin, 越負越深)---\n');
for k=[1 3 6], fprintf('  %s: 最深 %.1f um、最淺 %.1f um\n', plabel{k}, min(margin(k,:)), max(margin(k,:))); end

% ---- 總結:每極在整段 50-500 內 in-iron 的比例 ----------------------------
fprintf('\n=== 含 bias 電荷在磁極鐵內?(掃 R=%d:%d:%d, %d 點)===\n', R_um(1), R_um(2)-R_um(1), R_um(end), nR);
for k = 1:6
    nin = nnz(in_iron(k,:));
    fprintf('  %s: in-iron %3d / %d 個 R  (%s)\n', plabel{k}, nin, nR, ...
        ternary(nin==nR,'全程在內', ternary(nin==0,'全程在外','部分')));
end
allin = all(in_iron(:));
fprintf('\n>> 6 顆是否「每個 R 都在磁極內」: %s\n', ternary(allin,'是','否'));
fprintf('\n--- 每極「掉到磁極外」的 R [um] ---\n');
for k = 1:6
    outR = R_um(~in_iron(k,:));
    if isempty(outR), fprintf('  %s: (none)\n', plabel{k});
    else, fprintf('  %s: %s\n', plabel{k}, mat2str(outR)); end
end
firstAllIn = find(all(in_iron,1), 1);                       % 第一個「6 顆全 in」的 R index
if ~isempty(firstAllIn)
    contig = all(all(in_iron(:,firstAllIn:end),1));
    fprintf('\n>> 第一個 6 顆全 in 的 R = %d um;此後到 500 是否全 in: %s\n', ...
        R_um(firstAllIn), ternary(contig,'是','否'));
end

% ---- 代表性表:R=50/150/300/500 的 actuator 座標(um)+ in/out ------------
showR = [50 150 300 500];
for rr = showR
    ri = find(R_um==rr,1);  if isempty(ri), continue; end
    E = reshape(S.Esave(:,ri),3,6);  charge_act = S.ell_R(ri)*(Pc_base+E)*1e6;  % um
    fprintf('\n--- R=%d um (ell_hat=%.3f mm) | actuator 座標 [um] ---\n', rr, S.ell_R(ri)*1e3);
    fprintf('   pole    x_a      y_a      z_a    |r|um   in-iron\n');
    for k = 1:6
        rmag = norm(charge_act(:,k));
        fprintf('   %-4s %7.1f %8.1f %8.1f %7.1f   %s\n', plabel{k}, ...
            charge_act(1,k), charge_act(2,k), charge_act(3,k), rmag, ternary(in_iron(k,ri),'YES','no'));
    end
end

function out = ternary(cond,a,b), if cond, out=a; else, out=b; end, end
