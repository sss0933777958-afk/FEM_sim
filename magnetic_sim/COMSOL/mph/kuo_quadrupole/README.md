# magnetic_sim/COMSOL/mph/kuo_quadrupole/ — Kuo Quadrupole COMSOL 模型

**用途**：Kuo 4-pole MEMS Quadrupole 的 COMSOL `.mph` 模型（建模迭代 + baseline + 求解版本）。

**內容**（實際 `.mph`）：
- `quadrupole_new.mph` — 目前主用大模型（含解，~360 MB；mphinterp/mpheval 座標用 **mm**，見 memory）。
- `quadrupole_new_split.mph` — 拆分版（精簡）。
- `0.46mm_baseline.mph` — L_P=0.46 mm baseline。
- `tweezer.mph` — 較早的整體 tweezer 模型。
- `monolithic_T55.mph` 及 `quadrupole_*_T55.mph` 系列（imported / with_air / with_coils / with_blocks / with_physics / solved / solved_blocks）— T=55 µm 漸進建模各階段。

**相關**：上層 model 索引見 `../README.md`；LiveLink 連線見 `../../../../.claude/rules/comsol-livelink.md`；mm 座標探測坑見 memory `feedback_comsol_mm_coord_pitfall`；ANSYS 對應設計見 `../../../ANSYS/main/`。
