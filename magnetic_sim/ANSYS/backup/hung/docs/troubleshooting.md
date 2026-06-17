# Hung Hexapole — Troubleshooting Log

每次執行前必須閱讀此文件，避免重蹈覆轍。

---

## KMODIF 移動單一 tip KP → 表面扭曲

**問題：** 用 KMODIF 只移動 pole tip 的 KP，導致連接的 cone 表面變成 non-planar / twisted。

**原因：** tip KP 連接多條 line 和 area（cone surface、fillet arc 等）。只移動一個 KP 不會自動重建 surface of revolution。

**解法：** 不能只移動 tip KP。需要移動整根 pole 的所有 KP（整體平移），或重建 pole 幾何。

---

## VSBV 連續切削 → 第二刀退化

**問題：** 對同一 volume 做連續 VSBV，第二刀報 topological degeneracy。

**原因：** 第一刀產生的切割面和第二刀的 cutting block 邊界在 hex center 附近重疊。

**解法：** 用 NUMSTR 強制 cutting block 編號，每次 VSBV 前 ALLSEL 但只 VSEL 目標 pole volumes + cutting block。

---

## BTOL=1e-6 存在 db 中

**問題：** BTOL 設定會存在 db 裡，RESUME 後仍然生效。導致後續 VADD/VOVLAP 失敗。

**解法：** 做完需要 BTOL=1e-6 的操作後，立即 `BTOL, DEFA` 重設。

---

## *get NUM MAX 在 VSBV 後不可靠

**問題：** VSBV 刪除 volume 後留下空位，新建 volume 填入空位。`*get, VOLU, NUM, MAX` 回傳的不是新建的 volume。

**解法：** 用 `NUMSTR, VOLU, 500+i` 強制指定新 volume 編號。

---

## VOVLAP 37 volumes 極慢

**問題：** VOVLAP,ALL 對 35 steel + 2 air = 37 volumes 需要 >10 分鐘，可能不會完成。

**解法：** 先 VADD 所有 steel 成 1 volume，再 VOVLAP 少量 volumes（Long2016 方法）。

---

## Pole flat-cut Z=0 → topological degeneracy

**問題：** cutting block Z=0（切削面在 pole 軸線上）和 VROTAT quarter boundary 重合。

**解法：** Z 偏移 +0.2mm（5 根 pole）或 +0.5mm（GP6）。

---

## Pole flat-cut Z=-0.1mm → tip 被切鈍

**問題：** 負 Z 偏移在 tip 端（cone R < 0.1mm）會切掉整個截面。

**解法：** 只用正 Z 偏移。GP5 需 +0.1mm，GP6 需 +0.5mm。

---

## _RETURN 被 CSYS 覆蓋

**問題：** BLOCK 指令的 _RETURN 在下一個 CSYS 指令後被覆蓋為 0。

**解法：** 不依賴 _RETURN，用 NUMSTR 預設編號。

---

## KMODIF 整體平移 pole → 共享 KP 扭曲 block

**問題：** 即使移動整根 pole 的所有 KP（不只 tip），pole 和 block/yoke 的共享 KP 也會被移動，導致 block 表面扭曲。

**原因：** Pole 的 block-end KPs 和 T-block volume 共享。KMODIF 移動共享 KP 同時改變兩個 volume。

**解法：** 不能用 KMODIF。必須**刪除舊 pole volumes → 用調整後的參數重建**。

---

## Pole 重建時參數公式產生錯誤的 tip 位置

**問題：** 用 POLE_TILT + end_r + POLE_TOTAL_LEN 公式重建 pole 時，tip 位置在 z≈-27mm，但原始 db 的 tip 在 z≈-18mm，差 9mm。

**原因：** 原始 db 不是用 full15.out 的參數公式建的（或參數已改變）。公式和實際 db 幾何不匹配。

**解法：** 重建時直接指定 tip 座標（從最佳化結果），不要用參數公式計算。end 座標從原始 db 的 KP 讀取。

---

## 變數名衝突（TIP_Z 等）

**問題：** db 中已有 *DIM 定義的陣列（如 TIP_Z），用同名 scalar 會報 "needs subscripts"。

**解法：** 所有 flat-cut 變數加 `fc_` 前綴避免衝突。

---

## COIL_H vs COIL_DZ 混淆

**問題：** `COIL_H = 14mm` 是 block 頂到 yoke 底的間距（建模用），`COIL_DZ = 15mm` 是 SOURC36 coil 截面高度（Real constant RVAL2）。兩者名稱相似但意義完全不同。

**解法：** 程式碼中明確區分：
- `COIL_H` → 僅用於 yoke Z 位置計算（`fc_yoke_zbot = fc_up_blk_top + COIL_H`）
- `COIL_DZ` → 用於 Real constants 和 coil Z 位置（`block_top + COIL_DZ/2`）

---

## Upper block BLK_T/sin(TILT_UP) 對淺角度失效

**問題：** 原始 upper block 定位用 `fc_dt = BLK_T/sin(TILT_UP)` 計算 block 中心 XY。當 TILT_UP < 28° 時，fc_dt 過大，block 中心超出 yoke 外徑（62.5mm），pole 無法接上 block。

**原因：** 淺角度時 pole 幾乎水平穿過 block，需要很長的水平距離才能上升 BLK_T 的高度。

**解法：** 改用 pole end 居中定位（跟 lower block 一樣）：
```apdl
fc_zbot = fc_endz - BLK_T/2
fc_ztop = fc_endz + BLK_T/2
fc_cx = fc_endx
fc_cy = fc_endy
```
此方法不依賴 sin(TILT_UP)，任何角度都能正確定位。

---

## IGES unit flag 與 SolidWorks STEP 轉換

**問題：** ANSYS 匯出 IGES 時 unit flag = 6（mm），但實際數值為英寸（MM=1/25.4）。SolidWorks 開啟時忽略 flag 直接讀取為英寸 ×25.4 = 正確 mm。但另存為 STEP 時可能因 flag 不一致導致尺寸錯誤。

**解法：** 匯出後修正 unit flag 為 1（inches）：
```bash
sed -i "s/,1.0,6,,/,1.0,1,,/" file.iges
```

---

## Pole tip fillet：方案 A vs 方案 B

**問題：** 為 pole tip 加 40 µm 直徑 fillet 時，幾何上不可能同時滿足三個條件：
1. Fillet 在原 tip 位置
2. Cone 半角維持 11.31°（原始 STEP 量值）
3. Cone-cyl 接點維持在 X=15.875 mm

必須有一個讓步。

**方案 A（已採用）：smooth tangent fillet**
- Cone 半角不變（11.31°）
- Junction X 從 15.875 → 15.793 mm（內縮 0.082 mm）
- Fillet 與 cone 平滑相切（無銳角）
- Cylinder 段補回 0.082 mm，總長 43 mm 不變
- **這是 Long2016 的做法，幾何上正確**

**方案 B（不採用）：sharp corner**
- Cone 半角變成 12.96°（變陡）
- Junction X = 13.74 mm
- Cone 從 fillet foremost (0,0) 直接出發，與 fillet 形成銳角
- 改變了 pole 的物理形狀

**Long2016 公式（複用）：**
```
POLE_TIP_R = 20e-6                     ! 半徑（40 µm 直徑）
ANG2 = atan(POLE_R/POLE_CONE_LEN)      ! 11.31°
L_FILLET = POLE_CONE_LEN*POLE_TIP_R/POLE_R
POLE_TIP_CENTER = L_FILLET/cos(ANG2)
FILLET_DROP = POLE_TIP_CENTER - POLE_TIP_R = 0.082 mm
JUNC_X = POLE_CONE_LEN - FILLET_DROP = 15.793 mm
```

**實作：** `apdl/geom/MT_Hung_Assembly_Dfillet.txt`、`apdl/geom/export_pole_filleted.txt`

---

## Post-processing 從獨立 session 執行時找不到結果檔

**問題：** 用獨立 ANSYS session 跑 post_extract_wp.txt 時報 "No results file"。

**原因：** 預設 jobname 為 `post1`，但結果檔是 `coil1.rmg`。RESUME 找的是 `post1.db` 而非 `sim1.db`。

**解法：** 在 post script 開頭加 `/FILNAME, coil1`，確保 jobname 與結果檔一致。
```apdl
/FILNAME, coil1
RESUME, 'sim1', 'db'
/POST1
SET, LAST
```
