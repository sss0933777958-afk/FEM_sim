# apdl/long2016_hexapole_halfcut/mesh/ — 純 mesh（網格產生）腳本 + 尺寸表

**用途**：把「**mesh 產生**」跟「**solve**」分開的起點。mesh deck 在此建幾何＋劃網格＋`SAVE` 成 `.db`；
solve deck（在 `../sim/`）`RESUME` 那個 `.db` → 只加電流/邊界/材料 → 解（不重劃）。
本夾收**通用、各模型可套**的 region-graded mesh 腳本與尺寸規格（先以 long_fei hexapole halfcut 為對象）。

**內容**：
- `mesh_size_table.pdf` — **部件 → 網格 size 目標對照表**（參考目標，非 deck 硬尺寸；deck 用 smrt5+EREFINE 逼近）：
  內球 0.3／中球殼 1.5／磁極 cone 0.3（下極更細）／yoke 1.5／外圓柱 ~4mm。實際達成：內球~0.15／中球~1.2／yoke~0.6／外圓柱~3.7／cone~0.3mm。
- `MT_Mesh_Graded.txt` — **canonical region-graded mesh-only deck（已定案 2026-06-25、跑通驗證）**。一跑即產
  **~656,251 節點 / 3,925,014 元素** 的 `mesh_graded.db`（無 /SOLU、無 MAG BC，供 solve deck `RESUME`）。三步 pipeline：
  ① 純 `SMRTSIZE,5 + VMESH,ALL`（零強制尺寸）= 自適應 base ~511k（**鐵律：此幾何任何硬性 ESIZE/AESIZE/LESIZE 都 thrash，只有純 smrt 能跑**）；
  ② `EREFINE` 全極 cone +1 級（→~579k）；③ `EREFINE` 下極 P1/P3/P6 cone 再 +1 級（→~656k，下極 2 級、其餘極 1 級）。
  幾何 merged V8 照 baseline、單一 VOVLAP、**材料/分區一律 by location+size**（booleans 清 VATT、不可硬編 volume 號）。
  支撐座 + coil 鐵柱(protrusion) 維持 smrt5 粗（far-field、不影響精度；EREFINE 在超粗外圈空氣中造不出過渡、無法局部細化）。
  .db 輸出 → `ANSYS_data/long2016_hexapole_halfcut/db/mesh_graded/mesh_graded.db`（該夾**只留此定案 .db**）。

- `MT_Mesh_Baseline.txt` — **baseline smrt5 mesh-only deck（2026-06-26 抽離）**。= `sim/baseline/MT_Sim_P1.txt` 的 geom+材料+6 SOURC36 coil+VATT+`SMRTSIZE,5`+`VMESH,ALL` 逐行切出 + `SAVE`（去 BC/solve/POST1）。重現 canonical `db/mesh_baseline/mesh_baseline.db`（~494,889 節點、含 coils；**.db 歸 db/，2026 整理後與 data/ 分離**）。**baseline + gap sweep 共用此 .db**（同網格＝可比）。為 mesh/solve 分離的起點（`sim/mueq_sweep/` 與 `sim/gap100um_mueq/` RESUME 它）。

- `MT_Mesh_LowerFilled_smrt5.txt` — **下極填滿（full cone）+ 純 baseline smrt5 mesh-only deck（2026-06-26）**。= `MT_Mesh_LowerFilled.txt`（graded）**移除 2× EREFINE**（保留 mid-sphere + robust VATT-by-location）→ 純 smrt5 base。`SAVE db/lower_filled/mesh_lowerfilled_smrt5.db`（**不蓋** graded `mesh_lowerfilled.db`）。供 `sim/lower_filled/MT_Sim_LowerFilled_smrt5_coil1.txt` RESUME。

- `MT_Mesh_SinglePole.txt` — **單極模型** mesh-only deck（下極填回完整圓錐 + 支撐座 + 鐵柱 + 空氣域）。**均勻鐵件 ESIZE 0.3mm** + 空氣漸變（內球 0.4mm / 外圈 4mm）；材料 by location+size；建 SOURC36 coil（placeholder，solve deck 設 1A）→ `SAVE db/singlepole/mesh_singlepole.db`。供 `../sim/singlepole/MT_Sim_SinglePole.txt` RESUME。（注意：此為單極**均勻**網格，與上面 hexapole graded 不同。）

**慣例**：本夾只留最終 PDF（`.tex/.aux/.log` 編完即清，比照 `../../../matlab/.../results/`）；
mesh `.db` 不放這裡（很大、gitignore），跑出來一律歸 `ANSYS_data/long2016_hexapole_halfcut/db/<meshname>/`（**.db 全在 db/；`data/mesh/<variant>/` 只放 `csv/` 內的 sensor_local CSV**，2026 整理後 db↔data 分離）。

**慣例（mesh 腳本收斂）**：本夾的 `MT_Mesh_Graded.txt` 是**往後唯一 canonical mesh 產生腳本**；
舊散落的 mesh-only deck（`sim/mesh/`、`geom/mesh/`）已刪除（既有輸出 `.db/.cdb` 保留、但不再可重新產生）。

**相關**：solve 決見 `../sim/`（mesh 後的 solve / `sim/resolve/` RESUME 範例）；
分離可行性與 `.db` 角色見對話 / memory。3 種網格曾用：smrt5 baseline（baseline+gap+sweep 共用）、region-graded（本夾）、sensor_spheres（由 `sim/sensor_spheres/` deck 自建）。
