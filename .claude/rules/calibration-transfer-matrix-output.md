# 校正結果 PDF 的內容

適用：
- Flux：`matlab/Flux/{APDL/Calibration_using_FEM_modeling,Maxwell}/function/emit_results.m`（兩分支要同步）
- Force：`matlab/Force/Maxwell/function/emit_force_results.m`

## Flux 印出項目（依序）

| 電流版（current） | 電壓版（voltage） |
|---|---|
| ᴮĝ_I [mT/A] | ᴮĝ_V [mT/mV] |
| K̄_I | M̄ |
| ᴮĤ_I = ᴮĝ_I · K̄_I [mT/A] | ᴮĤ_V = ᴮĝ_V · M̄ [mT/mV] |
| ℓ̂ [µm] | ℓ̂ [µm] |
| ê | ê |
| 校正誤差 RMS [mT] | V |
| | 校正誤差 RMS [mT] |

## Force 印出項目（依序）

Flux 的全部項目，再加上力增益與力轉移矩陣：

| 電流版（current） | 電壓版（voltage） |
|---|---|
| ᴮĝ_I [mT/A] | ᴮĝ_V [mT/mV] |
| K̄_I | M̄ |
| ᴮĤ_I = ᴮĝ_I · K̄_I [mT/A] | ᴮĤ_V = ᴮĝ_V · M̄ [mT/mV] |
| ᶠĝ_I = ᴹĝ_B/(2ℓ̂)·(ᴮĝ_I)² [pN/A²] | ᶠĝ_V = ᴹĝ_B/(2ℓ̂)·(ᴮĝ_V)² [pN/mV²] |
| ᶠĤ_I = √(ᶠĝ_I) · K̄_I [√pN/A] | ᶠĤ_V = √(ᶠĝ_V) · M̄ [√pN/mV] |
| ℓ̂ [µm] | ℓ̂ [µm] |
| ê | ê |
| 校正誤差 RMS [pN] | V |
| | 校正誤差 RMS [pN] |

## 共通格式

- 矩陣符號用 hat：LaTeX 標籤 `{}^{B}\hat{H}_{I}`、`{}^{B}\hat{H}_{V}`、`{}^{F}\hat{H}_{I}`、`{}^{F}\hat{H}_{V}`。
- ê 是無因次偏移（以 ℓ̂ 為長度單位），用標準 `bmatrix`、無欄列標籤、無單位；single（`USE_BIAS=false`）時 ê ≡ 0，不印。

## 校正誤差

只印校正誤差 RMS，不印 NMAE、RMSPE：

```
RMS = sqrt(J / (3 · M · N))
```

J 為校正點上的殘差平方和（擬合 cost），3 為分量數，M 為激發數（六極 = 6），N 為校正點數；即每個殘差分量的均方根。Flux 的單位是 mT，Force 是 pN。
