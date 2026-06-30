# …/Hall_sensor_base_fix_dir/code/main/ — 求 d 主程式

**用途**：`main.m` — Hall-sensor 每極 `d` 的指定流程 driver。config 在頂部（R_select=150 µm、I=1 A、S_hall=130 V/T）。

**流程（四步）**：
1. **拿 ℓ̂**：載入 fix_dir fit_KI_fixl 解出的 `../../../fix_dir/data/fit_fixl_R<RRR>um.mat` 的 `ell`（規則#2；≈0.856；不再 fminbnd、不用 fit_KI_ball sweep 版）。
2. **抽電壓 V=S·B**：載 6-coil FEM（`wp` 建殘差場 B / `all` 抽 sensor 電壓）→ `build_sensor_geometry` → `extract_Vmat`（真實節點、沿 n 圓柱、all-source）。
3. **解 d / 殘差**：`solve_d`（內外層雙重加總、閉式解 d，**無 g_H**：`d = (Σ_j V_j(Σ_i SᵀS)V_j)⁻¹(Σ_j V_j Σ_i Sᵀb)`）→ `sensor_residual` 回 cost `J = Σ‖ε‖²`（模型 b=S·V·d）。
4. **存解**：`calib_sensor_d.mat`（本組 `../../data/`，規則#2；含 d/Vmat/exc_sign/ell_hat/J，**無 gH**）。

**預期數值（R=150 µm）**：ℓ̂≈0.856 mm（fix_dir fit_KI_fixl）、cost J≈0.142 T²（全域座標、上極 CAD 傾角 36.59° 錐面 + 下極磨平面）、d 為 no-gain（~1e-2）。

**`main_interp.m`（內插版）**：同流程，但 step2 改用 `extract_Vmat_interp`（standard 粗網格、sensor 圓柱內
均勻 1000 點真·FEM tet 重心內插）抽 Vmat → 解對角 d；存 `../../data/calib_sensor_d_interp.mat`（規則#2）。用來處理 standard
粗網格「sensor 圓柱內 0 真實節點」的情況（baseline 沒 sensor 加密時）。

**`calib_gap100um.m`（gap100um 校正 driver）**：與 main.m 同流程，差別：(1) 讀 `gap100um_mueq` 變體
（μ_eff=56 氣隙，P2←P1 應翻負）；(2) sensor 讀值用 `extract_Vmat_interp` **圓柱內均勻插 100 點平均、六顆統一**
（網格拓樸用 standard CSV、場用 gap100um）；(3) 解寫到本組 `../../data/`（`.mat`+`.txt`，規則#2/#3；
排版成 `results/calib_gap100um.pdf`，`results/` 只留 PDF）。sensor 位置維持現狀（下極仍在磨平面；錐面搬移待後續）。
預期：ℓ̂≈0.860 mm、cost J≈0.071 T²、V(感測 P2,激發 P1) 為**負**（−1.88e-3 V，μ56<翻負臨界 ~72）。
> PDF 排版：用一次性 xelatex（`results/` 只留最終 PDF，`.tex/.mat/.txt` 編完即清；Vmat 欄重排成 paper P1--P6、對角=self）。

**命名 / 慣例**：主程式組 → `code/main/{main.m, main_interp.m, calib_gap100um.m}`；模型數學一律在 `../function/`。

**相關**：見上層 `../README.md`、`../function/README.md`。
