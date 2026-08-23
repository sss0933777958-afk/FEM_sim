# model_check/long_fei/ — Long Fei 六極半切 檢查用幾何（IGES + STEP）

**用途**：`../../IGES/long2016_hexapole_halfcut/` 經單位/旗標轉換後、可被 ANSYS `IGESIN` 正確讀入的 Long Fei 6 極下極半切 hexapole 幾何。**新交付一律出 STEP**（per deliver-step-for-check）。

**gap 變體 STEP（交付檢查用；皆由 `apdl/.../geom/export/MT_Geom_BaseGap.txt` 的 `GAP_CFG` preset 產、`../scripts/make_gap_step.py <name>` 轉 mm STEP，中繼 IGES 走 scratch `ANSYS_data/.../db/geom_hexvariants/`）**：
- `long2016_hexapole_gap_300um.step`（2026-07-12）：6 極 base 磁隙都往 tip 移、厚 **300µm**、截面 conformal（下極 T 形 45mm³/極、上極 y±11 矩形 66mm³/極）。7 solids（鋼 + 6 base slab）。
- **`long2016_hexapole_gapdist.step`（2026-07-13，分布式）**：把「均勻 400µm base 磁阻」拆成 **base 100µm ＋ 導柱 protrusion 兩端各一 gap**（下極 z=−7/0 厚 78.6µm、上極 z=2/9 厚 86.87µm、圓 r=5mm）。**19 solids**（鋼 32838 + 6 base slab[下15/上22mm³] + 12 導柱圓盤 conformal r=5＝**截面 78.54mm²**[下6.17/上6.82mm³]）。OCC 驗：導柱盤 @radial47.5 z=−7/0/2/9、面積 78.54mm²、base @x42.5/40.95。
- **`long2016_hexapole_gap_400um.step`（2026-07-13，均勻對照）**：只 base 磁隙、厚 **400µm**、無導柱 gap。7 solids（鋼 + 3 下極 T 60mm³ + 3 上極 88mm³）。
- gapdist 與 gap_400um 總磁阻等價（使用者串聯磁阻設計）；μ_r=1 屬 sim 階段（STEP 不帶材料）。

**極尖倒圓半徑變體 IGES（交付檢查用；研究「加工尖端半徑 → lumped 模型預測能力」；scenario A＝錐體固定、尖端隨半徑後退）**：
- **`long2016_hexapolehalfcut_tip400um_cnc.iges`**：極尖倒圓 **400µm**（一般 CNC 車削，鈍尖）。尖端後退 → tip-to-WP=1.90mm、tip-to-tip≈3.8mm。
- **`long2016_hexapolehalfcut_tip20um_edm.iges`**：極尖倒圓 **20µm**（微 EDM，尖）。tip-to-WP≈0.44mm、尖端略前進。⚠ 原要 15µm，但 15µm 倒圓弧僅 19µm（模型 1.4 萬分之一）→ metre-deck boolean 建不出、IGES 匯出退化（SolidWorks 看破碎）；**放寬到 20µm** 兼顧「代表 EDM 尖尖」+ 乾淨可檢查可 mesh。
- 由 `apdl/.../geom/export/MT_Geom_Export_mm.txt` 改 `POLE_TIP_R`（0.4 / 0.020）重建（scenario A：flank/base 用 `POLE_TIP_R_REF=40µm` 凍結、尖端 kp 後退 `(R−REF)×4.099`；40µm 重現 baseline）。flag 6→2 patch（同 no_gap.iges，SolidWorks 讀）。`.db`（FEM-ready）在 `ANSYS_data/.../db/geom/tip_variants/tip{400,20}um.db`（+ tip40um baseline，皆只留 .db）。
- ⚠ **不出 STEP**：ANSYS IGESOUT → OCC 對此模型（大量 WPLANE/WPROTA）會**非均勻變形**（xy×1.411、z 另比例）；現有 gap `*.step`（make_gap_step 產）**同樣 xy=149.58mm 錯**（只驗過 z/體積、沒量 xy）。ANSYS .db 幾何本身已驗證正確（scenario A：40µm tip@0.5mm、400µm 後退到 1.90mm、15µm 0.42mm）。要 mm STEP 須走 OCC-primitive 重建（deliver-step 規則），不可用 ANSYS IGES→OCC。

**Sensor 加密球位置檢查 IGES（2026-08-05）**：
- **`sensor_refine_p1p2.iges`**：鐵件總成 ＋ **3 顆 R0.3mm 純空氣標記球**，用來在建網格前目視確認加密球位置。
  球心（mm，ANSYS 全域）＝貼附點 + 0.46·n̂（＝取樣圓柱幾何中心；0.41 氣隙 + 半個柱高 0.05），**SOFF = 4.5mm**：
  **P1 削平面** `[4.9082, 0, −12.5400]`、**P1 底錐面** `[4.7306, 0, −14.3336]`、**P2 上錐面** `[−3.0839, 0, −8.7754]`。
  三顆皆離鋼面 0.46mm、R0.3 → **清鋼 0.16mm**；球**不參與布林**（保持獨立 solid 才看得出有無吃鋼）。
- 由 `apdl/.../geom/export/MT_Geom_Export_mm_SensorRefine.txt`（＝ `MT_Geom_Export_mm.txt` verbatim + 3 球 + 球心診斷）產。
  APDL 端已驗：3 顆質心與設計值完全相符、單顆體積 0.113096 mm³、0 error / 0 warning。
- **inch 量級（`VLSCALE 1/25.4`）+ flag 6→1 patch**（＝ `MT_Geom_Export_mm.txt` 原生做法，
  該 deck 註明 2026-07-06 更正：SW 忽略 flag 當 inch 讀 → inch 座標 + INCH flag = 正確 mm；
  使用者 2026-08-05 確認採此路線）。同 `gap_200um` / `tip400um_cnc`。
  ⚠ 本夾旗標**歷史上不統一**：`no_gap`(May)、`vp_coil1`(Aug 4) 是 mm/flag 2；`gap_200um`(Jul 8)、
  `tip400um_cnc`(Jul 22) 與**本檔**是 inch/flag 1。**要跟產生它的 deck 一致**，不可一律套同一個數；
  `MT_Geom_Export_mm.txt` 的註解與本檔 §命名/慣例那句「MKS 用 flag 2」尚未統一（待日後定案）。
- ⚠ 同上不出 STEP（OCC 對此模型非均勻變形）。開檔第一件事：量 yoke 外徑應 ≈ **106mm**。
- 後續：確認球位後才改 `apdl/.../mesh/MT_Mesh_SensorRefine.txt` 的 `SC_X/Y/Z`（公尺）跑 lv1/2/3 收斂。

**內容**：`long2016_hexapole_no_gap.iges`（主）、`..._Geom_HollowProt.iges`、`..._Geom_HollowProt_Plain.iges`、`..._Geom_gap200um.iges`、`..._Geom_sphtip.iges`、`..._Geom_WPsphere.iges`（鐵件總成 ＋ WP 7mm 空氣球殼，raw 重疊、CAD 檢視用；mm/flag-2），另含 `Geom_WithCoil.iges`（6 coil rings）、`Geom_hp_split.iges`、`SinglePoleFilled.iges`（**單一下極：削平半錐填回完整圓錐** ＋ 4 塊支撐座 ＋ 1 根 protrusion 鐵柱；無 yoke / 無上極 / 無 coil 實體；mm/flag-2）。

**`..._Geom_WPsphere.iges` 出處**：由 `apdl/long2016_hexapole_halfcut/geom/export/MT_Geom_Export_mm_WPsphere.txt`（＝ `MT_Geom_Export_mm.txt` 鐵件 mm 建構 ＋ APDL `SPHERE` 加 WP 球，不做布林）`IGESOUT` → flag 6→2 patch。球心 z = −12.71 mm、R = 7 mm，與 FEM `V7` 同。

**`SinglePoleFilled.iges` 出處**：由 `apdl/long2016_hexapole_halfcut/geom/export/MT_Geom_Export_mm_SinglePole.txt`（＝ `MT_Geom_Export_mm.txt` 改：只建 1 下極 ＋ **略過 half-cut BLOCK+VSBV → 完整圓錐**、只建該極 1 根 protrusion、無 yoke、無上極、無 coil）`IGESOUT` → flag 6→2 patch。只出 mm 版（使用者要求）。

**`..._Geom_P1BaseGap.iges` 出處**（**檔名 legacy，現含全 6 極**）：由 `apdl/long2016_hexapole_halfcut/geom/export/MT_Geom_Export_mm_P1BaseGap.txt`（＝ `MT_Geom_Export_mm.txt` 於最終 VOVLAP 後加「6 極迴圈」：每極在旋轉 LOCAL frame `VSEL` 選該極 base ＋ 旋轉 yz `WPLANE` 兩面 ＋ `VSBW,ALL` 雙刀切）`IGESOUT` → flag 6→2 patch。**6 極 base 各 200µm 磁隙 = 貼合 base 實形的 conformal 實體 slab（之後 μ_r=1）**：下極 P1/P3/P6(ROTA 0/120/240) local x[42.5,42.7]、每極 5 vol@30mm³；上極(ROTA 60/180/300) local x[41,41.2]、每極 1 vol@44mm³(blockB)。**VSBW 分割真實鋼、非挖洞**（鋼總量不變）；slab 可日後每極 `VSEL,LOC` 指派 μ_r=1。⚠ protrusion 柱未切、仍跨接；只在 base 做隙。⚠ 磁隙用 VSBW 實體 slab、不要 VSBV void。圖示 `figures/long2016_hexapole_halfcut/geom_P1basegap/{fullmodel_top,fullmodel_iso,p1_base_gap_topview,p1_base_gap_iso}.png`。metre 版 → `IGES/long2016_hexapole_halfcut/..._Geom_P1BaseGap.iges`（flag 6）。μ_r=1 為 sim 階段（IGES 不帶材料）。

**資料來源 / 流向**：由 `../../IGES/long2016_hexapole_halfcut/`（+ `CAD_model/long_fei/STEP` 直轉的含 coil / split 版）轉換而來 → `apdl/.../geom` 用 `IGESIN` → ANSYS。

**命名 / 慣例**：⚠ 同一物理模型在 `model_check/` 用 topic 名 **`long_fei`**，但在 `IGES/`、`ANSYS_data/`、`apdl/` 仍用 **`long2016_hexapole_halfcut`**。與 `IGES/` 對應、須同步。MKS 轉換用 **flag 2 / 重 export**，不可抄 hung sed 6→1。`WithCoil`/`hp_split` 為由 STEP 直轉的額外變體（IGES/ 側無同名件）。

**相關**：`../README.md`、`../../IGES/long2016_hexapole_halfcut/`、`reference/workflows/iges-sync-quick.md`。
