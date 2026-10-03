# 座標系

## 場的擬合一律在 actuator frame

電荷模型、K̄_I、ℓ̂、ĝ 的載入 → 擬合 → 呈現，一律在 actuator（磁極軸）座標系，不可停在 Maxwell / ANSYS 原始座標。

- **六極**：平移 −SPH_OFST 後用 `R_act` 旋轉位置與場：
  ```matlab
  dhat  = tip ./ vecnorm(tip);                    % 3x6 極尖單位向量
  R_act = [dhat(:,1), dhat(:,3), dhat(:,5)].';    % P1、P3、P5 當列
  P_act = (R_act * P_meas.').';
  B_act = (R_act * B_meas.').';
  ```
  `R_act` 必為正交旋轉：`det(R_act)` ≈ 1，不成立就停下查極軸。
- **單極**：模型直接把極軸建在 +x、原點在中心（`SPH_OFST = 0`），不旋轉。
- 讀別人的 `.mat` 或場之前，確認對方也是 actuator frame。

## 例外：電壓路徑留在 Maxwell 全域座標

算電壓（V 矩陣）的整條路徑在 Maxwell 全域座標做，不做任何平移：sensor 中心與法線、中心面取樣網格點都是 Maxwell 座標，`.fld` 直接內插、直接投影 `b·n̂`（V 是純量，與框無關）。

## 畫圖

只有畫圖時才平移 z（把中心移到原點）。

座標系與號誌（source / sink）是兩回事，不可混為一談。
