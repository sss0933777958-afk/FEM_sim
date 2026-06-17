# semulator-process

跑 SEMulator 3D MEMS 製程模擬(PolyMUMPS / 自訂),從 mask layout 產 3D 結構,
給後續 STEP/IGES 匯出 → APDL。

> ⚠ **現況:paused at v2 / GUI-only**(license 缺 `COV_Expeditor`,無 CLI batch;
> 詳見 memory `project_semulator_harrison`)。本 SOP 重點在 GUI 流程 + resume 清單。

## 何時用

- 想驗 idealised APDL 幾何 vs 真實製程後的 3D 形狀差異
- 嘗試新製程組合(光罩 / 沉積 / 蝕刻 step)

## 輸入

| 參數 | 範例 |
|---|---|
| `{topic}` | `kuo_quadrupole`(目前唯一) |
| `{model_file}` | `magnetic_sim/ANSYS/main/semulator/{topic}/model_file/harrison.zam`(或 `tweezer.zam`) |
| `{process_file}` | `magnetic_sim/ANSYS/main/semulator/{topic}/process/harrison/*.vproc` |
| `{layout}` | `*.gds`(可用 Python 生成) |

## 既有資產(`magnetic_sim/ANSYS/main/semulator/kuo_quadrupole/`)

| 元件 | 路徑 | 狀態 |
|---|---|---|
| baseline | `baseline/PolyMUMPS.{zam,vproc,vmpd}` + `output/` + `result.png` | ✅ 跑通(2026-05-12) |
| custom model | `model_file/{harrison,tweezer}.zam`、`__semu_*.py` | ⚠ Harrison v2 paused |
| custom process | `process/harrison/*.vproc` + `*.gds` + `*.py` + notes | 📝 in progress(2026-05-14~29) |
| results | `results/` | 🗑 gitignored |
| scripts | `scripts/` | 🔲 空 |

## 步驟(GUI-only)

1. **開 SEMulator**:Coventor / Lam 預設安裝
2. **載 `.zam`**:File → Open → `{model_file}`
3. **載 `.vproc`**:Process → Import → `{process_file}`(若不在 `.zam` 內)
4. **跑 build**:Tools → Build 3D → 等(分鐘級)
5. ⏸ **檢查點 → [`model-check.md` §1](model-check.md#1-solidworks-視覺檢查)** 變體:
   3D viewer 內目視確認形狀(沉積厚度、蝕刻角度、對稱性)
6. **匯出 STEP / STL**:File → Export → 給後續 `cad-export.md` / `step-to-apdl.md`
7. **更新 notes**:`process/harrison/_NOTES.md`(沒有就建)

## 資料流接續

```
[semulator-process] → STEP/STL → [cad-export] → IGES → [step-to-apdl] → 參數表
                                                                          ↓
                                                            [apdl-geom-build] (對比 idealised)
```

## 產物

- [ ] 3D model 在 SEMulator 內存檔(`.zam` 更新)
- [ ] `magnetic_sim/ANSYS/main/semulator/{topic}/results/<build_tag>/` 輸出
- [ ] (可選)STEP / STL export → `magnetic_sim/ANSYS/main/CAD/{topic}/STEP/`(走 cad-export)

## License / 限制

- 有:`COV_SEMu3D`(GUI build)
- 沒:`COV_Expeditor`(CLI batch / scripting)→ **無法 headless 跑**
- 工作流必須有人在電腦前操 GUI

## 常見坑

- 想 CLI batch → license 不支援,死路
- `.vproc` 跨機器移植要 update material database `.vmpd` 引用
- 量級 µm 設計裡光罩 line width 設錯 → build 出無意義形狀
- 沒先存 `.zam` 就改流程 → SEMulator 不保證 undo 完整

## Resume tasks 清單(本 SOP 落地時)

- [ ] Harrison v2 製程定稿(目前 paused)
- [ ] `zhang_quadrupole` 是否要套同 framework 還沒決定
- [ ] `scripts/` 目錄想放什麼還沒想清楚(GUI-only 限制下可能用不到)

## 適用 topic

目前只 `kuo_quadrupole`;framework 可擴 `zhang_quadrupole` / 未來新設計。
