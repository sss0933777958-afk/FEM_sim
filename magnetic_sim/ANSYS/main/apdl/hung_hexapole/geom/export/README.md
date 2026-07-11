# apdl/hung_hexapole/geom/export/ — 幾何 IGES 匯出腳本（APDL primitive 重建）

**用途**：hung_hexapole 用 **native APDL primitive 一個一個部件重建**（匯入 CAD 走不通：MAPDL 布林失敗/SpaceClaim 卡/Workbench 非 MAPDL 法），純建幾何 → `IGESOUT` 出 .iges 供檢查/下游 mesh。直接複用 backup 驗證過的建法。

**內容**：
- `MT_Geom_Poles.txt` — **只建 6 磁極**（D 形錐-圓角圓柱：前 28mm D 形 VROTAT180 + 後 15mm 全圓 VROTAT360、40µm 尖端圓角；上 tilt35°/下 tilt5.71°；尖端在 R0.5mm magic-angle 球面）。→ `hung_poles.iges`。使用者已驗證正確。
- `MT_Geom_FullAssembly.txt` — **整個組立**（97 vol）：6 磁極 + 6 支撐座 block（上 25×22×10/下 22×22×10 L 形，VADD）+ 6 線圈 solid（Ri6/Ro10/h14）+ yoke（環 Ri38/Ro62.5/厚2）+ 3 導柱（下 R4/h46）+ 3 上核（上 R4）。→ `hung_full_assembly.iges`。**直接複用 `backup/hung/apdl/geom/MT_Hung_Assembly_Dfillet.txt` 全文**，只改 CWD/IGESOUT。
- `MT_Geom_FullAssembly_BaseGap.txt` — 上檔 + **6 支撐座 block 各切 200µm BLOCK-ONLY 磁隙**（→ sim 設 μ_r=1）。每極旋轉 `LOCAL,15` frame、`VSBW` 兩旋轉 yz 面在 local radial R & R+0.2mm 切；**只選 block（radial x[40,60] 排除 pole 鐵棒 back-cyl 心~35、z-band 排 guide-post/core/coil/yoke）→ 不切 pole 鐵棒/導柱**。下極 azim 0/120/240 @ **radial 48mm**（**移到 pole 鐵棒外**：鐵棒沿軸延到 radial~43.5mm，48 在純 Part2 延伸段，氣隙不穿磁極）、上極 60/180/300 @ radial 44.43mm（已在上極鐵棒外）。每極 slab=1 vol/20mm³(驗證)。CWD→`db/geom_basegap`、IGESOUT→`IGES/hung_hexapole/hung_full_assembly_basegap.iges`(flag6)→ sed 6→1 → `model_check/hung_hexapole/`。（單位英寸 MM=1/25.4；座標 mm×MM；不出 figures。）

**單位慣例（SolidWorks 相容）**：內部用**英寸值**（`MM=1/25.4`）；ANSYS `IGESOUT` 寫 flag 6 → `model_check/` 版 patch 成 **flag 1(inch)** 自洽（`sed 's/,1.0,6,,/,1.0,1,,/'`），SolidWorks 開起來為正確 mm。

**流向**：raw IGES → `IGES/hung_hexapole/`；patch 後 → `model_check/hung_hexapole/`（兩邊同步）。ANSYS scratch 跑在 `ANSYS_data/hung_hexapole/db/geom_*`，`/EXIT,NOSAV` 不留 .db、跑完清空。

**⚠ 註**：此 primitive 組立為**簡化模型**（6 導柱/核 + 6 線圈 solid），與 CAD STEP 的 24 導柱/無線圈 solid 不完全一致；磁極/yoke/支撐座對得上。CAD 真實幾何參考 `model_check/hung_hexapole/hung_hexapole_full.iges`（OCC 匯入版）。

**相關**：`../README.md`、`backup/hung/docs/`（建法/陷阱）、`.claude/rules/{apdl-editing,ansys-cad-alignment}.md`。
