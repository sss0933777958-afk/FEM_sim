function [Vmat, exc_sign] = extract_Vmat(results_root, cnst, apdl_to_paper_idx, ...
                                         sensor_pos, sensor_n, disc_u, disc_v, disc_local, Ndisc, S_hall)
% EXTRACT_VMAT  抽 v_jk:Vmat(i,j) = S_hall * disc 平均((−B)·n_i)。
%   列 i = sensor 極 P1..P6、欄 j = 激發 coil(sim)。回傳 all-source 翻號後的 Vmat
%   與翻號向量 exc_sign(下極激發 P1/P3/P6 → −1)。
%   用 coil1..coil6(v4 baseline,'all' dataset),sensor 在遠處錐面故用全域。
%   需 import_ansys_data 在 path。 (ported verbatim from calib_fem.m PAGE 2)
    Vmat = zeros(6,6);                                         % sensor 電壓矩陣 [V]
    for kc = 1:6                                               % 逐 coil(= sim j = kc)
        cn = sprintf('coil%d', kc);                           % coil 資料夾名
        da = import_ansys_data(fullfile(results_root, cn, 'standard'),'all',cn);  % 載 'all' 全域(sensor 在遠處錐面)
        zc = da.z - cnst.SPH_OFST;                            % ANSYS z → WP 框
        Fx = scatteredInterpolant(da.x,da.y,zc,da.bx,'linear','nearest');  % Bx 內插器
        Fy = scatteredInterpolant(da.x,da.y,zc,da.by,'linear','nearest');  % By
        Fz = scatteredInterpolant(da.x,da.y,zc,da.bz,'linear','nearest');  % Bz
        for i = 1:6                                           % 逐 sensor 極
            acc = 0;                                          % disc 上 (−B)·n+ 累加
            for p = 1:Ndisc                                  % 逐 disc 取樣點
                xp = sensor_pos(:,i) + disc_local(p,1)*disc_u(:,i) + disc_local(p,2)*disc_v(:,i);  % 取樣點全域座標
                Bv = [Fx(xp(1),xp(2),xp(3)); Fy(xp(1),xp(2),xp(3)); Fz(xp(1),xp(2),xp(3))];        % 該點 B
                acc = acc + Bv.' * sensor_n(:,i);            % 物理 signed B·n+(朝 n+ 正、朝 n− 負)
            end
            Vmat(i,kc) = S_hall * acc / Ndisc;               % 面積平均 × S_hall = sensor 電壓 [V]
        end
        fprintf('Page2: coil%d sensor 電壓抽取完成\n', kc);   % 進度
    end

    % ---- all-source 慣例:翻「下極激發」(P1/P3/P6)的欄,使每顆激發極都當 source ----
    %  物理:翻線圈繞線方向 ⟺ 把該 coil 的 B 整個變號(磁靜場對電流線性)→ 後處理翻號
    %  與「反向繞線重跑 FEM」bit-for-bit 等價,免重跑。翻後 self 對角全正、off-diag 幾乎全負。
    exc_sign = ones(1,6);                                    % 各 sim(coil)激發極的 source 翻號
    for j = 1:6                                              % 下極 P1/P3/P6 激發 → 翻號變 source
        if ismember(apdl_to_paper_idx(j), [1 3 6]), exc_sign(j) = -1; end
    end
    Vmat = Vmat .* exc_sign;                                 % all-source sensor 電壓(翻下極激發欄)
end
