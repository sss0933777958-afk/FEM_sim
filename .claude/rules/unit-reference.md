# 單位

討論、圖軸、結果 PDF、變數註解一律用下表。來源：`magnetic_sim/ANSYS/main/reference/Unit Reference Sheet/pdf/Unit_Reference_Sheet.pdf`。

| 量 | 符號 | 單位 |
|---|---|---|
| 有效長度 | ℓ̂ | µm |
| 位置、偏移、座標 | — | µm |
| 磁通密度 | b(p) | mT |
| 電流增益 | ᴮĝ_I | mT/A |
| 電壓增益 | ᴮĝ_V | mT/mV |
| 電流 | I | A |
| 電壓 | V | mV |
| 磁阻 | R_a | A/Wb |
| 粒子磁化常數 | ᴹĝ_B | A·µm²/mT |
| 力增益（電流） | ᶠĝ_I = ᴹĝ_B/(2ℓ̂)·(ᴮĝ_I)² | pN/A² |
| 力增益（電壓） | ᶠĝ_V | pN/mV² |
| 力 | F | pN |
| 磁常數 | k_m = µ₀/4π | 10⁻⁷ |

- 10⁰ 因子不標；無因次量不加單位。
- APDL 幾何內部仍以 mm 建、ANSYS 輸出 Tesla —— 那是內部單位，對外表述時換算成上表。
