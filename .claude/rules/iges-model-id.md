# IGES = 模擬模型識別規則

當使用者貼一個 `IGES_converted/<...>.iges` 路徑(或 `IGES/<...>.iges`)時,**Claude 必須從路徑識別出這是哪個物理模擬模型**,然後把這個模型當作「當前討論模型」context 進入後續對話。

## 強制行為

✅ **該做**:
- 看到 IGES 路徑 → 直接從第二層 topic 名 + 檔名認出物理模型(例如:Long Fei 半切六極 / Kuo Quadrupole / Zhang Quadrupole)
- 用模型物理特徵回應(「這是 Long Fei 6 極半切 hexapole 的 baseline 幾何」)
- 後續對話預設這個模型是 context

❌ **不該做**:
- 問使用者「這是哪個模型?」/「對應到哪個 topic?」
- 問「要對應到哪個 result dir / analysis dir?」(那是 [[main-workspace]] 規則的事)
- 自動把 IGES 路徑展開成 5 個 sub-folder(使用者沒要求)

## 模型清單(2026-05-30)

依 `magnetic_sim/ANSYS/main/IGES_converted/` 第二層 topic 名:

| IGES 路徑樣本 | 物理模型 |
|---|---|
| `magnetic_sim/ANSYS/main/IGES_converted/long2016_hexapole_halfcut/*.iges` | Long Fei 6 極**下極半切** hexapole(V4 主軸)|
| `magnetic_sim/ANSYS/main/IGES_converted/long2016_hexapole_full/*.iges` | Long Fei 6 極完整 hexapole |
| `magnetic_sim/ANSYS/main/IGES_converted/long2016_p1_only/*.iges` | Long Fei P1 + yoke 部分極 |
| `magnetic_sim/ANSYS/main/IGES_converted/long2016_p1p2_yoke/*.iges` | Long Fei P1+P2 + yoke |
| `magnetic_sim/ANSYS/main/IGES_converted/long2016_p1p2p3_yoke/*.iges` | Long Fei P1+P2+P3 + yoke |
| `magnetic_sim/ANSYS/main/IGES_converted/long2016_dipole_tilted/*.iges` | Long Fei 上下對極 tilted dipole |
| `magnetic_sim/ANSYS/main/IGES_converted/long2016_dipole_lower/*.iges` | Long Fei 下對極 dipole |
| `magnetic_sim/ANSYS/main/IGES_converted/long2016_h1h2/*.iges` | Long Fei H1/H2 sensor 比值(DC)|
| `magnetic_sim/ANSYS/main/IGES_converted/long2016_h1h2_harmonic/*.iges` | Long Fei H1/H2 AC(渦電流)|
| `magnetic_sim/ANSYS/main/IGES_converted/kuo_quadrupole/*.iges` | Kuo 4-pole MEMS Quadrupole(主 Lp046 + V2 ScaleDown 變體)|
| `magnetic_sim/ANSYS/main/IGES_converted/zhang_quadrupole/*.iges` | Zhang 4-pole 0/90/180/270° |
| `magnetic_sim/hung/IGES_converted/*.iges` | Hung Hexapole(獨立 design root)|

新增模型時,在本檔清單加一列即可。

## 無 context 時

使用者說「**這個 / 上次的 / 那個模型**」但對話無前文 → 主動問:

> 請貼 IGES_converted 路徑(或告訴我是哪個 topic),我才知道在講哪個模型。

不要自己猜。

## 觸發片語

- 任何 `.iges` 路徑出現在使用者訊息
- 「看看這個模型」/「這個 IGES」/「這個幾何」
- 「我現在在做 XXX 模型」(XXX 沒對應到清單時問)

## 相關

- [[main-workspace]](`.claude/rules/main-workspace.md`)— 各模型對應的 sub-folder 表(要找 results/analysis 才查這)
- [[kuo-iges-export-workflow]] — IGES_converted 同步規則
