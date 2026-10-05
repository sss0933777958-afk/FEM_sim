# FEM_sim

六極電磁微探針（hexapole magnetic tweezers）的 FEM 模擬與點磁荷模型校正。

- 求解器：ANSYS MAPDL（APDL）、ANSYS Maxwell
- 分析：MATLAB
- 模型：長飛 2016 半切六極（主力）、志鵬平板六極、hung、NTU

## 資料夾

```
magnetic_sim/ANSYS/
├── main/                    活躍工作區
│   ├── CAD_model/           SolidWorks 原檔 + STEP（幾何的依據）
│   ├── model_check/         交付檢查用的 mm STEP
│   ├── apdl/<model>/        APDL 腳本
│   ├── ANSYS_data/<model>/  APDL 結果（.dat、.db，不進 git）
│   ├── matlab/
│   │   ├── Flux/            磁場模型校正（APDL / Maxwell 兩個分支）
│   │   └── Force/           力模型校正
│   ├── figures/             論文圖
│   └── reference/           推導、報告、論文
└── backup/                  歸檔的舊設計
```

Maxwell 專案與匯出的場（`.fld`）放在 repo 外：`D:\Maxwell_sim\<model>\`。

## 校正包的結構

以 `matlab/Flux/Maxwell/` 為例（Force 沒有自己的 `config/`，共用這裡的）：

| 資料夾 | 放什麼 |
|---|---|
| `config/` | 各模型的設定與幾何 |
| `main/main.m` | 校正主程式，改檔頂的參數切換模型與 current / voltage |
| `function/` | 校正用函式 |
| `data/` | 校正結果 `.mat`（只留定案版） |
| `results/` | 校正結果 PDF |
| `utils/scripts/<model>/`、`utils/data/<model>/` | 額外分析的計算腳本與結果（依模型分夾） |
| `plot/` | 只做繪圖的腳本 |
| `figures/<model>/{current,voltage}/` | 圖 |

## 怎麼跑

```matlab
run('magnetic_sim/ANSYS/main/matlab/Flux/Maxwell/main/main.m')   % 磁場校正
run('magnetic_sim/ANSYS/main/matlab/Force/Maxwell/main/main.m')  % 力模型校正
```

跑完會自動存 `.mat` 到 `data/`、出 PDF 到 `results/`。

## 環境

- ANSYS 2025 R2（MAPDL、Maxwell）
- MATLAB R2025b

## 參考

Fei Long, *Design, Fabrication, and Calibration of a Hexapole Magnetic Tweezers*, PhD Dissertation, Ohio State University, 2016.
