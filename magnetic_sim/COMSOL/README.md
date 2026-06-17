# magnetic_sim/COMSOL/ — COMSOL 求解器子層

**用途**：磁學模擬中以 COMSOL Multiphysics 求解的部分，與 `../ANSYS/` 並列為同一磁學類別下的另一求解器。

**內容**：
- `mph/` — COMSOL `.mph` 模型，依 model topic 分子夾（`kuo_quadrupole/`、`long2016_hexapole_halfcut/`）。

**相關**：上層磁學類別見 `../README.md`；COMSOL LiveLink headless 連線通用法見 `../../.claude/rules/comsol-livelink.md`（拆兩 process + `mphstart(2036)`）；各 model 模型清單見 `mph/README.md`。
