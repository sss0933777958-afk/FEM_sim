# ANSYS_data/long2016_hexapole_halfcut/db/ — ANSYS 模型 .db / sim scratch

**用途**：放 ANSYS `.db`（已 mesh / 已解的模型），給 **GUI Resume 互動檢視**、或**不重解重抽**
新物理量 / 新區域（RESUME → POST1 再匯出）。也當部分 sim 的 scratch 工作夾。

**內容**（子夾依用途）：
- `sensor_spheres/` — sensor 加密球 mesh 的 `.db`（`meshsens.db`）。
- `graded_p2/coil1/` — graded 密網格 P2-region 解的 `.db`（保留供重匯更大區域 / ESEL 抽鋼件外框）。
- `allsource/`、`coils_same_dir/` — GUI 檢視用小 `.db`（看繞向 / solid-model 鐵件；build-only 不解）。

**慣例（half-clean）**：sim 解完**立刻**把 `.dat` 複製到 `../data/coil<N>/<variant>/`，
scratch 只保 `.db`（+ 主 log），**刪掉** `.esav/.full/.rmg/.dbb/.page/.lock/.DSP*/.bat` 等中間檔
（一顆 solve 中間檔 16–18 GB，務必清，不可囤積）。探路性、已抽完的 scratch 可整夾刪
（`.dat` 已在 `../data/`）。清理 SOP：`.claude/rules/sim-cleanup.md`。

**注意**：`.db`/`.rmg` 等 ANSYS 輸出**不進 git**（見 repo `.gitignore`）。

**相關**：見 `../data/README.md`、`../RESULTS_MAP.md`、`.claude/rules/sim-cleanup.md`。
