# iges-sync-quick

改了單一 IGES 零件後,快速同步到 `IGES_converted/`。
**不是首次 export** — 那走 [cad-export.md](cad-export.md)。

## 何時用

- 改了個別零件 IGES,要同步轉檔(不重新 export 全套)
- 修 IGES 單位 flag(誤抄 hung sed 後復原)

## 輸入

| 參數 | 範例 |
|---|---|
| `{topic}` | `kuo_quadrupole` / `long_fei` |
| `{iges_files}` | 變動的 `.iges` 清單 |

## 前置

- 原始 `magnetic_sim/ANSYS/main/IGES/{topic}/<file>.iges` 已更新

## 步驟

### kuo(MKS-metre — **絕大多數**)

從 SolidWorks **重新 export 該零件**,單位選 mm,直接放到
`magnetic_sim/ANSYS/main/IGES_converted/{topic}/`(取代舊檔)。

或在 SW 不開的情境下,改 IGES header 單位 flag 為 2:
```powershell
# (rough — Linux sed 不直接適用 Windows;用文字編輯器改 header
# 「Unit Flag」整數欄為 2,或重 export)
```

> ⚠ **不可抄 hung 的 sed 公式** `s/,1.0,6,,/,1.0,1,,/` —— 那是 hung 1/25.4 inches 才對。
> kuo MKS-metre 抄了會炸 1000×(memory `feedback_iges_unit_conversion`)。

### hung(inches — 參考用,不在本 workflow scope)

```bash
cp magnetic_sim/hung/IGES/Part.iges magnetic_sim/hung/IGES_converted/Part.iges
sed -i "s/,1.0,6,,/,1.0,1,,/" magnetic_sim/hung/IGES_converted/Part.iges
```

(本 SOP **只管 kuo**;hung 的請參 `.claude/rules/hung-docs.md`)

## 產物

- [ ] `magnetic_sim/ANSYS/main/IGES_converted/{topic}/<file>.iges` 與 `magnetic_sim/ANSYS/main/IGES/{topic}/<file>.iges` 同步
- [ ] 第三方 viewer(FreeCAD / SW)再 import 驗 bbox(走 `model-check.md` §2)

## 常見坑

- 抄錯設計的 sync 公式(hung vs kuo)→ 尺寸 1000× / 0.001×
- 改了 `IGES/` 忘記同步 `IGES_converted/` → APDL import 用舊形狀
- 刪除 / 重命名一邊 → 兩邊要同步處理(`main-workspace.md` 強制規則 #3)

## 適用 topic

所有 kuo topic;**hung 不適用本檔**。
