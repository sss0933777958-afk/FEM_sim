# function/single_pole

單極（一片極板、一顆線圈、一次激發）校正用的函式，由 `../../main/single_pole/main.m` 呼叫。
檔名一律 `sp_` 開頭，避免和 `../hexapole/` 的同名函式在 MATLAB path 上互相遮蔽。

| 檔 | 用途 |
|---|---|
| `sp_config.m` | 載 `config/<model>/<geom>/mt_constants.m`（保留 config 的 `N_I`、要求 `strategy='single_pole'`） |
| `sp_load.m` | 讀 `.fld` → 極尖座標系（原點＝極尖、極軸 +x）、mT、all-source；建鐵節點遮罩（尖端倒圓＋尖端楔形）、三線性內插、取樣點的「八個鄰點皆空氣」檢查 |
| `sp_fit.m` | 點電荷擬合：b = G·S(p/ℓ̂ − (1, e_y, e_z))，ℓ̂ 為沿 x 離極尖的距離，e_x 固定 0，G 閉式解 |
| `sp_emit.m` | `.mat` → PDF（ᴮĝ_I、ℓ̂、ê、RMS），落點 `results/<model>/<single|eighteen>/current_R<µm>_<geom>.pdf` |

目前只有 current 版；單極沒有 sensor，voltage 版未做。
