# bs-matrix-derive

從 N×N FEM 解(N=6 hexapole / 4 quadrupole)推出三層轉移矩陣:
**B̄_S(T) → B_S(V/A) → V_out/V_in(V/V)**。

## 何時用

- 新設計做完全 coil 自激/跨激 FEM,要算 sensor 輸出矩陣
- 改 sensor 規格(`S_hall`、`I_in`、`k_A`)要重算

## 輸入

| 參數 | 範例 |
|---|---|
| `{topic}` | `long2016_hexapole_halfcut` |
| `{coil_results}` | N 個 `magnetic_sim/ANSYS/main/ANSYS_data/{topic}/coil[1..N]/`(全 coil 自激解齊) |
| `{S_hall}` | 130 V/T(預設;Long Fei) |
| `{I_in}` | 0.6 A(預設) |
| `{k_A}` | 0.36 A/V |

## 既有腳本鏈

| 步驟 | 腳本(`magnetic_sim/ANSYS/main/analysis/long2016_hexapole_halfcut/fit/`) | 輸出 |
|---|---|---|
| 1 | `compute_B_bar_matrix.m`(disc 721 點 + 正確 coil1 = `coil1_pre_fine_mesh`) | `B_bar.mat`、`Bbar_S.mat` |
| 2 | `gen_Bbar_S_latex.m` | `Bbar_S_matrix.tex`(T,bmatrix) |
| 3 | `gen_Bs_latex.m` | `Bs_matrix.tex`(V/A) |
| 4 | `gen_Vout_Vin_latex.m` | `Vout_Vin_kA0p36.tex`(V/V) |
| 診斷 | `diag_coil1_variants.m` | coil1 資料夾陷阱檢測 |

## 前置

- N 個 coil 都跑完 FEM([apdl-fem-run.md](apdl-fem-run.md) × N)
- sensor placement 已定([h1h2-analysis.md](h1h2-analysis.md) 同樣 placement);
  Long Fei 用 cone surface + 0.41 mm normal

## 步驟

1. **驗 coil 資料夾正確**:`diag_coil1_variants.m` 跑一次,確認 coil1 對(P1 自激
   ≈ 0.014 T 量級,跟兄弟下極 P3/P6 一致);**hexapole halfcut 的 `coil1` 資料夾
   是 uncut 幾何,正確 = `coil1_pre_fine_mesh`**(memory `project_long_fei_B_bar`)
2. **跑 step 1**:`compute_B_bar_matrix.m` → `B_bar.mat` / `Bbar_S.mat`
3. **產 3 份 LaTeX**:依序跑 `gen_Bbar_S_latex.m` → `gen_Bs_latex.m` → `gen_Vout_Vin_latex.m`
4. ⏸ **檢查點 → [`model-check.md` §3](model-check.md#3-參數表-sanity-check)** 變體:
   對角線符號 + off-diag 結構合理?(下極對角 +、上極對角 −)
5. (可選)LaTeX 編譯:`magnetic_sim/ANSYS/main/pdf/{topic}/scripts/` → `magnetic_sim/ANSYS/main/pdf/{topic}/doc/`

## 產物

- [ ] `magnetic_sim/ANSYS/main/ANSYS_data/{topic}/B_bar.mat` + `Bbar_S.mat` + `Bs.mat` + `Vout_Vin_kA0p36.mat`
- [ ] `magnetic_sim/ANSYS/main/pdf/{topic}/scripts/{Bbar_S,Bs,Vout_Vin_kA0p36}_matrix.tex`

## 常見坑

- coil1 資料夾陷阱(上述)
- raw FEM 保留物理號;**cosmetic 對角+/off-diag−** 只在 LaTeX 輸出層套,raw `.mat`
  不要硬壓號型(memory `feedback_sensor_sign_convention_toward_wp`)
- sensor disc 採樣點數改(37 ↔ 721)結果 bit-identical(對平滑場 no-op),但**新採樣
  pattern 都要驗 ≥ P2-P6 column 不變**
- 上極 SOURC36 winding 方向跟下極反 → 自激對角 raw 為負(未定論:真實設計 vs APDL bug)

## 適用 pole 配置

hexapole(N=6)/ quadrupole(N=4)。**dipole 不適用**(N=2 直接走 [h1h2-analysis.md](h1h2-analysis.md))。
