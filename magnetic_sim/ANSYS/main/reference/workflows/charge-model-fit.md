# charge-model-fit

從 FEM B 場擬合等效 monopole / dipole / multipole 強度,給力學模型用。
**必先跑 validity sweep**(memory `feedback_J_fit_validity_kuo_vs_hung`)。

## 何時用

- 新極形要建力學模型,需要等效電荷
- 比較不同極形(cylinder vs cone-fillet)有效範圍

## 輸入

| 參數 | 範例 |
|---|---|
| `{topic}` | `kuo_quadrupole` / `long2016_hexapole_full` |
| `{case_tag}` | 對應 FEM 結果 |
| `{sample_region}` | `cube` / `cylinder` / `sphere_shell` |
| `{r_max}` | 採樣半徑上限(µm) |
| `{pole_shape}` | `cylinder`(validity ~65%)/ `cone-fillet`(~15%)— **決定 r_max 安全範圍** |
| `{model_order}` | `monopole` / `dipole` / `quadrupole` / `multipole` |

## 前置

- FEM 跑完 + `apdl-postproc.md` 的 `xyz_extract` 拉出 `coord` + `bfield` `.dat`
- 已知極形分類(影響 validity)

## 樣板腳本

| 採樣形 | 既有 |
|---|---|
| cube | `magnetic_sim/ANSYS/main/analysis/kuo_quadrupole/fit/fit_J_cube70_and_tex.m`<br>`fit_J_F20_cube70.m` / `fit_J_F20_real4_cube70.m` |
| cylinder | `magnetic_sim/ANSYS/main/analysis/kuo_quadrupole/fit/fit_J_ScaleDown_V2_P1_cylinder.m` |
| sphere shell | `magnetic_sim/ANSYS/main/analysis/kuo_quadrupole/fit/fit_J_sphere_shell.py` |
| region scan | `fit_J_region_scan.m`(掃 r_max 找最佳)|
| 多極 | `magnetic_sim/ANSYS/main/analysis/kuo_quadrupole/fit/fit_J_quadrupole_R0500.{m,py}` |

## 步驟

1. **驗 FEM 資料齊**:`coord_all.dat` + `bfield_all.dat` 在 `magnetic_sim/ANSYS/main/ANSYS_data/{topic}/{case_tag}/coilN/`
2. **跑 validity sweep**(必!):掃 `r_max = 20..150 µm`,看 err vs r_max 曲線。
   - cylinder pole 預期 break ~65% pole-distance;cone-fillet 預期 ~15%
3. ⏸ **檢查點:err < 5% 才採信**;否則縮 r_max 或換 model
4. **正式 fit**:選定 r_max 跑對應 `fit_J_*.m`
5. **產出 LaTeX 表**:`magnetic_sim/ANSYS/main/analysis/{topic}/fit/gen_latex_hung_style.m` 風格
6. **存 .mat**:`magnetic_sim/ANSYS/main/ANSYS_data/{topic}/J_fit_{case_tag}.mat`

## 產物

- [ ] validity sweep 圖(err vs r_max)
- [ ] `magnetic_sim/ANSYS/main/ANSYS_data/{topic}/J_fit_*.mat`
- [ ] `magnetic_sim/ANSYS/main/figures/{topic}/{case_tag}/J_fit_*.png`
- [ ] (可選)LaTeX 表 → `magnetic_sim/ANSYS/main/pdf/{topic}/scripts/`

## 常見坑

- 跳過 validity sweep 直接 fit → cone-fillet 設計 r_max 抓太大 err 爆 30%+
- 非正交基底投影手算 dot-product → 用 `lstsq` 或矩陣 `λ = E⁻¹ v`
  (memory `feedback_nonorthogonal_basis_projection`)
- 採樣 box 跨 pole 邊界 → fit 不收斂;box 內必純氣域
- LaTeX 在 chat 顯示 → 用純文字數學式(memory `feedback_no_latex`)

## 適用 topic

所有 topic。**cone-fillet 極形必小心 r_max**;
參考 `magnetic_sim/ANSYS/backup/hung/figures/coil1/fit_J_validity_boundary_hung.png` 看 hung 邊界形狀。
