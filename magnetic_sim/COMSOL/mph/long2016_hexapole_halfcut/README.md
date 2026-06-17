# magnetic_sim/COMSOL/mph/long2016_hexapole_halfcut/ — Long Fei 半切六極 COMSOL 模型

**用途**：Long Fei 六極**下極半切** hexapole 的 COMSOL `.mph` 模型（DC / AC 渦電流 / B-H 飽和 變體與備份）。

**內容**（實際 `.mph`）：
- `long_fei_model_DC.mph` — DC 磁靜態（驗證 WP center |B|，coil1=P5 對應見 memory）。
- `long_fei_model_DC_BH.mph` — DC 加 B-H 非線性飽和（AISI1018，μ_r(B)）。
- `long_fei_model_AC.mph` — AC 頻掃 / 渦電流（導電鋼 σ=7e6）。
- `0.46mm_baseline_V2.mph` — baseline V2（大檔 ~930 MB）。
- `long_fei_model.bak.mph`、`long_fei_model_AC.bak.mph` — 備份版。

**相關**：上層 model 索引見 `../README.md`；LiveLink 連線見 `../../../../.claude/rules/comsol-livelink.md`；DC/AC/BH 結果與坑見 memory `project_long_fei_DC_comsol`、`project_long_fei_AC_comsol`、`project_long_fei_bh_saturation`、`feedback_comsol_multiturn_ac_artifact`；ANSYS 對應設計見 `../../../ANSYS/backup/hexapole-long2016/` 與 `../../../ANSYS/main/`（halfcut topic）。
