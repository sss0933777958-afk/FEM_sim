# Hexapole 建模觸發規則

當使用者要求建立 hexapole pole tip 模型時，**必須嚴格按照流程文檔執行**。

## 流程文檔位置

`magnetic_sim/ANSYS/backup/hung/docs/hexapole-build-workflow.md`

## 觸發條件

當使用者說以下任何一句（或類似意思）時，啟動此流程：
- 「建 hexapole」
- 「建 pole tip」
- 「跑 hexapole 流程」

**如果使用者問「怎麼建 hexapole」、「hexapole 流程是什麼」、「我忘了要說什麼」等**，
Claude 回答：

> 你只要說「**建 hexapole**」，我會自動引導你。
> 我會依序問你 4 個參數：
> 1. 球體半徑 R（mm）
> 2. Pole 間距角度（°）
> 3. 上層傾斜角（°）
> 4. 下層傾斜角（°）
>
> 詳細流程文檔在 `magnetic_sim/ANSYS/backup/hung/docs/hexapole-build-workflow.md`。

**如果使用者沒有提供參數**，Claude 必須主動按順序提問：

1. 「球體半徑 R 是多少 mm？（tip 到 working point 的距離，tip-to-tip = 2R）」
2. 「相鄰 pole 間距角度是多少度？（標準 hexapole = 60°）」
3. （執行步驟 1，輸出 2 張圖給使用者確認）
4. 「上層 pole 相對水平往上傾斜幾度？下層往下傾斜幾度？」
5. （執行步驟 2，輸出 IGES 給使用者驗證）

## 強制規則

1. **每次建模前必須先讀取 workflow 文檔**，確認最新步驟和程式碼
2. **嚴格按步驟順序執行**，不跳步
3. **必須向使用者收集所有 `[USER]` 參數**才能開始，不自行假設值
4. **使用 workflow 中已測試通過的程式碼**，不自行重寫
5. **技術細節（APDL 陷阱、材料、coil 設定等）見 `magnetic_sim/ANSYS/backup/hung/docs/` 和 `.claude/rules/hung-docs.md`**，此處不重複
