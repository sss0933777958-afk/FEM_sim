%% anim_charge_side.m
%  ============================================================================
%  兩個側視動畫(GIF):下極 P1 / 上極 P2 的含-bias 等效磁荷位置 vs 取樣半徑 R=50→500µm。
%  磁極輪廓 = CAD 幾何(STEP intent,同 plot_lower/upper_pole_sensor_placement.m):
%    POLE_R=3.175mm、POLE_CONE_LEN=15.876mm、half-angle 11.31°;
%    P1 = milled half-cone(磨平面 y=0 朝 WP、下半錐),軸水平;
%    P2 = full natural cone,軸傾斜 35.26°(=90−54.74)。
%  座標:apex 在原點(pole-local)。zoom 到 apex/磁荷區。
%  ----------------------------------------------------------------------------
%  ⚠ 幾何一致性註記:磁荷是擬合到 FEM 場(FEM 錐含 40µm tip、POLE_R=3.0mm);此圖用 CAD
%    幾何(sharp apex、3.175mm)。下極磁荷貼錐面,CAD(sharp)下會落在錐面外側 ~20µm,
%    FEM(含 tip)下則在內 ~13µm —— 差別來自 CAD-vs-FEM 的極尖定義。圖上照 CAD 畫。
%  ============================================================================
clear; clc; close all;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
model = 'long2016_hexapole_halfcut';

cnst = mt_constants();
S = load(fullfile('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\bias_fit\data','sweep_nofixl_vs_R.mat'));
fig_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
mm = 1e3;  nR = numel(S.R_um);  Rrot = S.R_act;  Pc_base = S.Pc_base;

% ---- 含 bias 磁荷(measure 框,mm)----
chg = zeros(3,nR,6);
for ri = 1:nR
    E = reshape(S.Esave(:,ri),3,6);
    chg(:,ri,:) = reshape(Rrot.'*(S.ell_R(ri)*(Pc_base+E)), 3,1,6);
end

% ---- CAD 幾何(mm)----
POLE_R = 3.175;  POLE_CONE_LEN = 15.876;  cyl_end = 30.0;

% ===== P1(下極,milled half-cone,軸水平)=====
k=1; tip = [cnst.pole_tip_x(k);cnst.pole_tip_y(k);cnst.pole_tip_z_wp(k)];
% pole-local:x = 沿水平極軸(measure +x)、y = 垂直(measure +z);磨平面 y=0、下半錐 y<0
xl1 = zeros(2,nR);
for ri=1:nR, v=(chg(:,ri,k)-tip)*mm; xl1(:,ri)=[v(1); v(3)]; end
outl1_x = [0, POLE_CONE_LEN, cyl_end, cyl_end, 0];
outl1_y = [0, -POLE_R,       -POLE_R, 0,       0];
make_anim(xl1, outl1_x, outl1_y, [0 cyl_end],[0 0], S.R_um, S.ell_R, ...
    'P1 (lower, milled half-cone)', [-0.3 1.7], [-0.55 0.25], ...
    fullfile(fig_dir,'anim_charge_P1.gif'), fig_dir,'P1');

% ===== P2(上極,full cone,軸傾斜 35.26°)=====
k=2; tip = [cnst.pole_tip_x(k);cnst.pole_tip_y(k);cnst.pole_tip_z_wp(k)];
tilt = deg2rad(90 - 54.7356);  u=[cos(tilt);sin(tilt)]; p=[-sin(tilt);cos(tilt)];
A=[0;0]; J=POLE_CONE_LEN*u; T=(POLE_CONE_LEN+14.124)*u;
Ju=J+POLE_R*p; Jl=J-POLE_R*p; Tu=T+POLE_R*p; Tl=T-POLE_R*p;
outl2_x = [A(1),Ju(1),Tu(1),Tl(1),Jl(1),A(1)];
outl2_y = [A(2),Ju(2),Tu(2),Tl(2),Jl(2),A(2)];
% pole-local(ref 框):x_ref = 水平徑向外(= −measure x)、y_ref = 垂直(measure z)
xl2 = zeros(2,nR);
for ri=1:nR, v=(chg(:,ri,k)-tip)*mm; xl2(:,ri)=[-v(1); v(3)]; end
make_anim(xl2, outl2_x, outl2_y, [A(1) T(1)],[A(2) T(2)], S.R_um, S.ell_R, ...
    'P2 (upper, full cone, tilt $35.26^{\circ}$)', [-0.3 1.7], [-0.25 1.35], ...
    fullfile(fig_dir,'anim_charge_P2.gif'), fig_dir,'P2');

fprintf('done.\n');

%% ---------------------------------------------------------------------------
function make_anim(xl, outl_x, outl_y, axln_x, axln_y, R_um, ell_R, ttl, xlim_, ylim_, gifpath, fig_dir, tag)
    nR = numel(R_um);  previewR = [50 275 500];
    f = figure('Color','w','Position',[100 100 900 470]);
    for ri = 1:nR
        clf; hold on; axis equal; box on;
        fill(outl_x, outl_y, [0.80 0.83 0.88],'EdgeColor','k','LineWidth',1.8);   % 磁極(CAD 輪廓)
        plot(axln_x, axln_y, '-.','Color',[0.6 0.6 0.65],'LineWidth',0.6);        % 極軸
        plot(0,0,'k.','MarkerSize',16); text(0,0,'  apex','FontSize',9,'VerticalAlignment','top');
        plot(xl(1,1:ri), xl(2,1:ri),'-','Color',[0.55 0.65 0.95],'LineWidth',1.0);% 軌跡
        plot(xl(1,ri), xl(2,ri),'o','MarkerSize',10,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
        xlim(xlim_); ylim(ylim_);
        xlabel('x [mm]  (along pole axis from apex)','FontSize',11); ylabel('y [mm]','FontSize',11);
        title(sprintf('%s   $R=%d~\\mu$m,  $\\hat{\\ell}=%.3f$ mm', ttl, R_um(ri), ell_R(ri)*1e3), ...
              'Interpreter','latex','FontSize',13);
        set(gca,'FontSize',10); drawnow;
        fr = getframe(f); [Aim,map] = rgb2ind(frame2im(fr),256);
        if ri==1, imwrite(Aim,map,gifpath,'gif','LoopCount',Inf,'DelayTime',0.08);
        else,     imwrite(Aim,map,gifpath,'gif','WriteMode','append','DelayTime',0.08); end
        if any(R_um(ri)==previewR)
            exportgraphics(f, fullfile(fig_dir, sprintf('anim_%s_frame_R%03d.png',tag,R_um(ri))),'Resolution',130);
        end
    end
    fprintf('wrote %s\n', gifpath);  close(f);
end
