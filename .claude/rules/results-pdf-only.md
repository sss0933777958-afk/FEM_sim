# 校正結果 `results/` 只放 PDF

適用各校正包的 `results/<model>/{single,eighteen}/`。

1. 只放 `.pdf`。`.mat` 放 `data/`；`.tex`、`.aux`、`.log` 在 xelatex 編完後即刪。
2. `main.m` 的最後一步就是呼叫 `emit_results` 產 PDF。
3. Maxwell 分支檔名：`<base>_R<半徑µm>[_N<點數>][_<tag>][_soff<N>mm].pdf`
   - `<base>` = `current` | `voltage`；`_R` 一律標
   - `_N<點數>` 只在降取樣時出現
   - `_<tag>` = variant 去掉分支名與 `convN<n>` 後的殘餘
   - `_soff<N>mm` 只在 voltage 且 sensor 距不是 4.572 mm 時出現
   - 例：`current_R150.pdf`、`voltage_R150_N80_soff3mm.pdf`
