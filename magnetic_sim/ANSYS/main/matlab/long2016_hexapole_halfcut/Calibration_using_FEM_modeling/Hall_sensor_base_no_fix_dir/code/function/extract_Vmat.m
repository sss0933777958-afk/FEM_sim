function [Vmat, exc_sign] = extract_Vmat(results_root, cnst, apdl_to_paper_idx, ...
                                         sensor_pos, sensor_n, S_hall, variant, sensor_r, axial_tol)
% EXTRACT_VMAT  從 6 顆線圈的 FEM 場抽 sensor 電壓矩陣 Vmat（真實節點，含 all-source 翻號）。
% -------------------------------------------------------------------------
% 用途（電壓轉換 V = S_hall·⟨B·n+⟩）：
%   對每次模擬 j（= 第 j 顆線圈通 1A），算每顆 sensor i 的輸出電壓：
%       Vmat(i,j) = S_hall · mean_over_selected_nodes( B·n_i )
%   B 直接取「真實 FEM 節點」的場（不做 scatteredInterpolant 內插）。
%
% 真實節點選點規則（使用者拍板：圓柱底面=sensor 面、往 n+ 長高 0.1mm，半徑 0.15mm）：
%   圓柱「底面」= sensor 中心 sensor_pos(:,i)（離鐵 0.41mm 的 air-gap 面），沿 n+ 往外(WP)長：
%     - 半徑   rho ≤ sensor_r（= Ø0.3mm 感測盤半徑 0.15mm）：in-plane 落在盤內
%     - 軸向 0 ≤ axial ≤ axial_tol（= 0.1mm 高，只往 n+）：全程離鐵 ≥0.41mm，不往鐵側下探
%   圓柱內的真實節點 → 對 (B·n_i) 平均，當作該 sensor 的面積平均近似。
%   若圓柱內「沒有」節點 → 在 rho ≤ sensor_r 且 axial ≥ 0 的節點中，取 axial 最小者
%   （n+ 側、沿 n 最近的節點）；連 n+ 側盤內都沒有才退回整體最近節點。
%   每顆 sensor 都 fprintf 命中節點數，避免靜默取到 0 點。
%
% 符號慣例：
%   - 取 B·n+：n+ 為「出鋼」方向，B·n+ 帶正負（朝 n+ 正、朝 n− 負）。raw FEM
%     場直接點乘 n+（本函式以 da.b* 原始場計算，不先取負）。
%   - all-source：把「下極激發」(P1/P3/P6)那幾欄整欄變號(exc_sign=−1)，使每顆被
%     激發的極都視為 source。物理上等價於「反向繞線重跑 FEM」(磁靜場對電流線性，
%     翻電流 ⟺ 翻整個 B)，故可後處理翻號、免重跑。翻後 self 對角全正、off-diag 幾乎全負。
%
% 輸入：
%   results_root        FEM 結果根目錄（其下有 coil1..coil6/standard）
%   cnst                幾何常數（用 SPH_OFST 把 ANSYS z 轉成 WP 框）
%   apdl_to_paper_idx   1×6，APDL coil j → paper pole 索引（判是否下極）
%   sensor_pos          3×6，sensor 中心（build_sensor_geometry）
%   sensor_n            3×6，sensor 法線 n+
%   S_hall              霍爾靈敏度 [V/T]（EQ-730L = 130）
%   variant    (選填)   FEM 結果變體子夾名，預設 'standard'（baseline）；sensor 加密用 'sensref'
%   sensor_r   (選填)   感測盤半徑 [m]，預設 0.15e-3（Ø0.3mm）
%   axial_tol  (選填)   圓柱高度 [m]（只往 n+，底面=sensor 面），預設 0.1e-3
% 輸出：
%   Vmat      6×6 sensor 電壓 [V]（列=sensor 極 P1..P6，欄=激發 coil1..6，已翻號）
%   exc_sign  1×6 各欄(激發極)的 all-source 翻號（下極激發 = −1，上極 = +1）
%
%   需 import_ansys_data 在 path。
% -------------------------------------------------------------------------
    if nargin < 7 || isempty(variant),   variant   = 'standard'; end % 讀哪個 FEM 變體子夾（預設 baseline）
    if nargin < 8 || isempty(sensor_r),  sensor_r  = 0.15e-3; end   % 感測盤半徑預設 0.15mm
    if nargin < 9 || isempty(axial_tol), axial_tol = 0.1e-3;   end  % 圓柱高度預設 0.1mm（只往 n+，底面=sensor 面，全程離鐵≥0.41mm）

    Vmat = zeros(6,6);                                         % sensor 電壓矩陣 [V]（先填未翻號的）
    for kc = 1:6                                               % 逐線圈 kc（= 第 j=kc 次模擬）
        cn = sprintf('coil%d', kc);                           % 該線圈結果資料夾名 coil<kc>
        da = import_ansys_data(fullfile(results_root, cn, variant),'all',cn);  % 載 'all' 全域真實節點（variant 子夾）
        X  = [da.x, da.y, da.z - cnst.SPH_OFST];             % 節點座標(WP 框)：ANSYS z 扣球心偏移 → Nn×3
        Bn = [da.bx, da.by, da.bz];                          % 節點場 B（raw FEM）→ Nn×3
        for i = 1:6                                           % 逐 sensor 極 i
            ni = sensor_n(:,i);                              % 該 sensor 法線 n+（3×1）
            r  = X - sensor_pos(:,i).';                      % 各節點相對 sensor 中心的位移（Nn×3，廣播相減）
            axial = r * ni;                                  % 沿 n+ 的軸向分量（Nn×1；>0 = n+ 側、離鐵更遠）
            rho   = vecnorm(r - axial*ni.', 2, 2);          % in-plane 徑向距離 ‖r − axial·n‖（Nn×1）
            sel = (rho <= sensor_r) & (axial >= 0) & (axial <= axial_tol);  % 圓柱：盤內 + 從 sensor 面往 n+ 0~axial_tol
            if ~any(sel)                                     % 圓柱內無節點 → 放寬軸向（仍只 n+ 側）
                cand = find((rho <= sensor_r) & (axial >= 0)); %   盤內且在 n+ 側的節點
                if isempty(cand)                            %   連 n+ 側盤內都沒有 → 取整體最近節點兜底
                    [~, j0] = min(rho); sel = j0;
                else
                    [~, k0] = min(axial(cand)); sel = cand(k0);  % 盤內、n+ 側沿 n 最近的那顆
                end
            end
            Bdotn = Bn(sel,:) * ni;                         % 選中節點的 B·n+（朝出鋼正負；列向量）
            Vmat(i,kc) = S_hall * mean(Bdotn);              % 真實節點平均 × 靈敏度 = sensor 電壓 [V]
            fprintf('  coil%d sensor P%d: 命中 %d 個真實節點\n', kc, i, numel(Bdotn));  % 透明回報命中數
        end
        fprintf('Page2: coil%d sensor 電壓抽取完成\n', kc);   % 進度訊息
    end

    % ---- all-source 慣例：翻「下極激發」(P1/P3/P6)那幾欄，使每顆激發極都當 source ----
    %  物理：翻線圈繞線方向 ⟺ 把該 coil 的 B 整個變號（磁靜場對電流線性）→ 等價於後處理翻號，
    %  與「反向繞線重跑 FEM」bit-for-bit 相同、免重跑。翻後 self 對角全正、off-diag 幾乎全負。
    exc_sign = ones(1,6);                                    % 各模擬(欄)激發極的 source 翻號，預設 +1
    for j = 1:6                                              % 若第 j 欄激發的是下極 P1/P3/P6 → 翻成 −1
        if ismember(apdl_to_paper_idx(j), [1 3 6]), exc_sign(j) = -1; end
    end
    Vmat = Vmat .* exc_sign;                                 % 套用翻號：得 all-source 慣例的 sensor 電壓矩陣
end
