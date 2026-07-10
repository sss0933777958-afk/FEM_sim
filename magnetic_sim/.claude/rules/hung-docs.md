# Hung Hexapole 必讀文件規則

當工作涉及 `magnetic_sim/ANSYS/backup/hung/` 目錄下的任何檔案時（建模、模擬、後處理、繪圖），**開始動手前必須依序讀取以下文件**：

1. `magnetic_sim/ANSYS/backup/hung/docs/troubleshooting.md` — 已知陷阱，避免重蹈覆轍（最優先）
2. `magnetic_sim/ANSYS/backup/hung/docs/hexapole-build-workflow.md` — 建模 + 模擬完整 SOP（步驟 1~7）
3. `magnetic_sim/ANSYS/backup/hung/docs/hexapole-sphere-geometry.md` — 球面幾何、座標系、tip 座標
4. `magnetic_sim/ANSYS/backup/hung/docs/pole-geometry.md` — pole 尺寸、傾斜角、正交驗證結果
5. `magnetic_sim/ANSYS/backup/hung/docs/hexapole-simulation-reference.md` — 通用模擬規範（材料、元素、求解器）

## 關鍵提醒

- APDL 乘號 `*` 前後**不能有空格**（會被當成註解）
- Block 定位用 **pole end 居中**（`fc_endz ± BLK_T/2`），不用 `BLK_T/sin(TILT_UP)`
- COIL_H（14mm, yoke 間距）和 COIL_DZ（15mm, SOURC36 截面高度）是**不同的參數**
- 9 Volumes：V1-6=tips, V7=sphere air, V8=steel body, V9=cylinder air
- IGES 匯出用 `MM = 1/25.4`（SolidWorks 相容），模擬用 `MM = 1e-3`（MKS）
- 零件尺寸變更時，更新 `magnetic_sim/ANSYS/backup/hung/apdl/geom/export_parts.txt` 裡的對應數值並重新跑 ANSYS 匯出到 `magnetic_sim/ANSYS/backup/hung/IGES/`（直接覆蓋舊檔）
- **`magnetic_sim/ANSYS/backup/hung/IGES/` 和 `magnetic_sim/ANSYS/backup/hung/IGES_converted/` 必須同步更新**：任何一個零件（`.iges`）在 `IGES/` 更新後，必須重新產生對應的 `IGES_converted/` 版本。流程：
  ```bash
  cp magnetic_sim/ANSYS/backup/hung/IGES/Part.iges magnetic_sim/ANSYS/backup/hung/IGES_converted/Part.iges
  sed -i "s/,1.0,6,,/,1.0,1,,/" magnetic_sim/ANSYS/backup/hung/IGES_converted/Part.iges
  ```
  刪除或重命名 `.iges` 時兩邊也要同步處理。不可只改其中一邊。

## 不涉及 magnetic_sim/ANSYS/backup/hung/ 時忽略此規則
