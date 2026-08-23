# 單位統一慣例（Unit Reference Sheet，強制）

**使用者拍板（2026-07-06）**：本專案（hexapole / single-pole lumped-parameter 建模與校正）**所有討論、圖、結果 PDF、變數註解的單位，統一照 Unit Reference Sheet**。

**source of truth**：`magnetic_sim/ANSYS/main/reference/Unit Reference Sheet/pdf/Unit_Reference_Sheet.pdf`
（原稿 `.../reference/Unit Reference Sheet/scripts/Unit_Reference_Sheet.tex`）

相關規則：`figure-style.md`（圖 + 結果 PDF 的數值標註：10^0 / 無單位不標）、`fit-current-matches-sim.md`。

---

## 🔒 單位表（lumped-parameter model）

| 量 | 符號 | 單位 |
|---|---|---|
| 有效長度 | ℓ̂ | **µm** |
| 磁通密度（WP 場） | b(p) | **mT** |
| 電流增益 | ᴮĝ_I | **mT/A** |
| 電壓增益 | ᴮĝ_V | **mT/mV** |
| 電流 | I | **A** |
| 電壓 | V | **mV** |
| 磁阻 | R_a | **A/Wb** |
| 力增益 | ᴹĝ_B | (µm)³·A²/pN |
| ᴹĝ_B/2ℓ̂·(ᴮĝ_I)² | ᶠĝ_I | (µm·mT)²/pN |
| 力 | m(p)、F | **pN** |
| 磁常數 | k_m=µ₀/4π | 10⁻⁷ |

## 🔒 規則
1. **長度一律 µm、場一律 mT、電流 A、電壓 mV、力 pN**（討論、圖軸、結果 PDF、變數註解都照此）。
   - 例：極尖端位置寫 **500 µm**（不是 0.5 mm）、ℓ̂ 報 **µm**、WP 場報 **mT**。
   - **位置 / 偏移座標也是長度 → µm**：電荷位置、bias 橫向偏移 `e_y/e_z`（＝物理偏移 ℓ̂·e，非無因次）、sensor/取樣點座標等，一律用 **µm**（例：單極 bias PDF `ℓ̂ / e_y / e_z` 三欄都 µm）。
2. **幾何物理大小不變**：APDL geom deck 仍以 mm-magnitude 建（500 µm = 0.5 mm），只是**對外表述 / 分析單位**統一成上表。
3. 數值標註沿用 `figure-style.md`：10⁰ 因子不標、無因次不加單位標記。

## 觸發片語
- 「單位」/「用 µm / mT」/「照 Unit Reference Sheet」
- 寫結果 PDF / 圖軸標 / fit 報告數值時

## 何時不適用
- 純 APDL 幾何 magnitude（內部 mm）、純 FEM raw `.dat`（ANSYS 輸出 Tesla、SI）——那是內部計算單位，**對外表述時換算成上表**。
