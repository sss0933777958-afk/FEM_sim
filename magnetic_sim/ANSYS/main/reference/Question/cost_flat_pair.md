# 對極結構單一激發、中心取樣時，cost 曲線在 ℓ̂ 大的一側變平

狀態：**待討論**（2026-10-06）。來源：組會簡報 `report/406/10.6/10.pptx` 第 3、4 頁。

## 1. 觀察

### 1.1 之前的對極三維校正（簡報第 3 頁，志鵬 Pair_pole）

- 場：`D:\Maxwell_sim\Zhi_peng\export\Pair_pole\{P1,P2}.fld`。
- 取樣：原點＝兩極尖中點，沿致動軸 xa（斜 35.26°）取 2·Nr+1 = 11 點（±150 µm、間距 30 µm），三維場。
- 模型：兩顆電荷 p̄c1 = û + ê1、p̄c2 = û + ê2，約束 ê1u − ê2u = 0；bᵢ(p) = Sᵢ·G，S 兩欄（兩顆電荷），G 閉式消去，
  min J(ℓ̂, ê, G) = Σbᵢᵀbᵢ − (ΣSᵢᵀbᵢ)ᵀ(ΣSᵢᵀSᵢ)⁻¹(ΣSᵢᵀbᵢ)。
- ℓ̂ 輪廓（固定 ℓ̂、重擬 ê，J/J_min，ℓ̂ = 500 ~ 856 µm）：

| | ℓ̂ | J/J_min @ 500 µm | J/J_min @ 856 µm |
|---|---|---|---|
| 一個激發（P1） | 678.4764 µm | 11.6 | **1.18**（高側幾乎平） |
| 兩個激發（P1、P2） | 623.5410 µm | 6.76 | 4.16 |

### 1.2 一維模型重做（簡報第 4 頁，NTU）

- 模型：H(x) = a·I/(x − ℓ̂)²，資料 |b_x|/I，給定 ℓ̂ 時 a 閉式，對 ℓ̂ 掃描得 profile cost J(ℓ̂)。
- 取樣：比照第 3 頁「中心取樣」，但軸線改水平——以極尖前 500 µm 為中心、±150 µm、17 點（x = −650 ~ −350 µm）。

| 配置 | 激發 | ℓ̂ | cost 曲線 |
|---|---|---|---|
| 單極 | P1 | 416.9 µm | 有明確極小值 |
| 對極（左極尖在 −1000 µm、同軸） | 只激發 P1 | 1560.5 µm | 由 ℓ̂ = 0 單調下降後**平坦**，J ≤ 1.1 J_min 的範圍 1348 ~ 1829 µm |

對照（簡報第 2 頁，取樣貼近極尖 x = −340 ~ −20 µm）：單極 ℓ̂ = 103.5 µm、對極 ℓ̂ = 123.0 µm，兩者都有明確極小值。

## 2. 問題

cost 曲線在 ℓ̂ 大的一側變平，目前只在下面三個條件**同時**成立時出現：

1. 對極結構（另一顆極存在）；
2. 只激發其中一顆極；
3. 取樣點在兩極中間（離被激發極尖遠），而不是貼近極尖。

為什麼？是資料本身（對極的場在中間區域的形狀）造成，還是模型數學（ℓ̂ 與 a、ê、G 的可辨識性）造成，或兩者皆有？

## 3. 可考慮的方向

- 未激發的那顆極是導磁的回流路徑：兩極中間的場是「被激發極 + 另一極回流」的疊加，沿軸的衰減可能比單一點電荷 1/d² 慢，單電荷模型只能用很大的 ℓ̂ 去配平緩的場。
- 可辨識性：取樣點離電荷越遠、取樣範圍相對距離越窄，1/d² 越接近一段緩坡，ℓ̂ 與 a（或 G）的敏感度越共線，谷底越平。
- 第二個激發：兩個激發時曲線兩側都被夾住（第 3 頁紅線），代表第二筆場資料提供了單一激發缺少的約束；待確認一維模型是否有同樣現象。
- 分開驗證：沿整條軸看單一激發的場形狀（最小值位置、局部衰減率），以及只改一個條件（排列、激發數、取樣位置）時曲線如何變。

## 4. 相關檔案（`matlab/Flux/Maxwell/` 下）

| 內容 | 腳本 | 資料 | 圖 |
|---|---|---|---|
| NTU 一維取樣與擬合 | `utils/scripts/NTU_hexapole/{sp_axis_calc,pp_axis_calc,sp_axis_fit}.m` | `utils/data/NTU_hexapole/{sp,pp}_axis[_c500].mat` | `figures/NTU_hexapole/{single_pole,pair_pole}/sampling[_c500].png` |
| NTU profile cost | `utils/scripts/NTU_hexapole/sp_cost_b.m`（`CASE`、`SAMP`） | `utils/data/NTU_hexapole/{sp,pp}_cost_b[_c500].mat` | `figures/NTU_hexapole/{single_pole,pair_pole}/cost_b[_c500].png` |
| NTU 殘差（擬合 vs FEM） | `plot/NTU_hexapole/current/sp_axis_plot.m` | `{sp,pp}_axis.mat` | `figures/NTU_hexapole/{single_pole,pair_pole}/b_axis.png` |
| 志鵬對極取樣、一維擬合 | `utils/scripts/zhi_peng/{zp_axis_fit,zp_slice_y0}.m`、`plot/zhi_peng/{pair_sampling,pair_field_xa}.m` | `utils/data/zhi_peng/zp_*.mat` | `figures/zhi_peng/pair_pole/*.png` |
| 志鵬對極三維輪廓（第 3 頁） | 10-01 的 `fitting_pair.m` 已於 10-03 刪除 | — | `figures/zhi_peng/pair_profile.png` |
