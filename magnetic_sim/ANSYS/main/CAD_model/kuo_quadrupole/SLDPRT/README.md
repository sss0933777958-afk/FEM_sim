# CAD_model/kuo_quadrupole/SLDPRT/ — Kuo Quadrupole SolidWorks 原檔

**用途**：Kuo 4-pole MEMS Quadrupole 的 SolidWorks 可編輯母檔（元件 + 組合件 + 縮小變體），幾何修改一律改這裡再重出 STEP/IGES。

**內容**：元件 `pole.SLDPRT`、`coil.SLDPRT`、`yoke.SLDPRT`、`guid_post.SLDPRT`；組合件 `quadrupole.SLDASM`；縮小變體 `kuo_scale_down.SLDPRT`、`kuo_scale_down_V2.SLDPRT`、`kuo_scale_down_fabricate.SLDPRT`。（`~$*` 為 SolidWorks 暫存鎖檔，勿動。）

**資料來源 / 流向**：母檔 → `../STEP/` 匯出 → 下游 IGES pipeline。

**命名 / 慣例**：`.SLDPRT` 元件、`.SLDASM` 組合件；變體後綴 `_V2`、`_fabricate`、`ScaleDown` 表示設計修訂/製程版本。

**相關**：`../README.md`、出檔 SOP `doc/workflows/cad-export.md`。
