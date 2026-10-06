# reference/Question/ — 待討論 / 待辦的研究問題

每個問題一份 `.md`：背景、問題、可考慮的方向、待決事項。問題解決後在該檔頂端更新狀態。

| 檔案 | 問題 | 狀態 |
|---|---|---|
| `noise_analysis.md` | 一維電荷模型 H = a/(x − b)² 在 FEM 資料加雜訊後，a、b 的變異數 | 已做（2026-10-05，N(0, 1) mT/A × 10000 次）：ℓ̂ 變異數單極 7.03 µm²、對極 5.85 µm²；腳本與結果在 `matlab/Flux/Maxwell/utils/{scripts,data}/NTU_hexapole/`；文件內容仍是原始問題 |
| `force_metrics.md` | 力模型能否定義類似磁通 𝒞 / κ 的 Efficiency index 與 isotropy | 待討論 |
| `cost_flat_pair.md` | 對極結構單一激發、中心取樣時，cost 曲線在 ℓ̂ 大的一側變平（10.6 簡報第 3、4 頁） | 待討論 |
