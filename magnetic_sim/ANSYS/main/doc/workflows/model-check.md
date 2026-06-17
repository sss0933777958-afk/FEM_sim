# model-check

4 節獨立的「使用者檢查模型是否正確」routine,被其他 SOP 在 ⏸ 點 reference。
本文件不重抄,只列**檢項 / 判準 / 怎麼做 / 不通過退回到哪**。

| 節 | 何時跑 | 主場 |
|---|---|---|
| §1 | cad-export.md Step 4 之後 | SolidWorks |
| §2 | cad-export.md Step 6 之後 | FreeCAD / SolidWorks |
| §3 | step-to-apdl.md Step 7 之後 | 文字編輯 + 心算 |
| §4 | apdl-geom-build.md Step 7 之後 | ANSYS APDL |

---

## §1 SolidWorks 視覺檢查

- **檢項**
  - 主尺寸:pole 長 / pole tip 半徑 / yoke 厚度 / WP 距離
  - 拓樸:零件數、mate 完整、無孤立 body
  - 對稱:quadrupole 0/90/180/270°、hexapole 0/60/120/180/240/300°、上下層 60° offset
- **判準**:全部 ± 0.1% 內;對稱旋轉預期角度後幾何重合
- **怎麼做**:SolidWorks → Tools → Measure;View → Display → Rotate 360°
- **不通過**:退回 SLDPRT 修,**不要** re-export STEP / IGES

---

## §2 IGES round-trip 檢查

- **檢項**:`IGES_converted/{topic}/{basename}.iges` re-import 後 bounding box 對得上原 SLDPRT
- **判準**:bbox 誤差 < 0.5%;不會「縮小 1000× 變看不見」或「放大 1000× 爆 viewport」
- **怎麼做**(任一):
  - FreeCAD 開 `.iges` → View → Fit All → 看右下角 bbox
  - SolidWorks open IGES → Tools → Measure 量整體 X/Y/Z
- **不通過**:回 `cad-export.md` Step 5 重做 IGES_converted
  - 常見成因:抄了 hung 的 sed 公式 / IGES unit flag 沒改

---

## §3 參數表 sanity check

- **檢項**
  - 所有 `trust=3D ✓` 行 ± 0.1% 內對得上 SLDPRT 量測值
  - 每行 `trust=2D ⚠` 標註處置:`keep`(其實是 3D 漏判)/ `drop`(2D 邊界,丟)/ `refine`(再 grep)
- **判準**:✓ 行全對 + ⚠ 行全處置
- **怎麼做**:把 `<basename>_params.md` 開在編輯器,左邊 SLDPRT measure,逐行對
- **不通過**:回 `step-to-apdl.md` Step 2-5 補抓 / 修 grep pattern

---

## §4 APDL 模型驗證

**自動檢查**(SOP 內貼 APDL 直接跑):

```apdl
! volume / element 數
VLIST
*GET, n_elem, ELEM, , COUNT
*GET, n_vol,  VOLU, , COUNT

! 鋼體積(對照 CAD)
VSEL, S, MAT, , 1                    ! 假設 MAT=1 為鋼
*GET, vol_steel, VOLU, , COUNT
VSUM
*GET, vol_total, VOLU, , VOLU

! coil 數(SOURC36 = ETYPE 11 or 自定)
ETLIST

! 邊界條件
! 在 .txt 內 grep "D,ALL,MAG,0" 確認存在
```

**判準**:
- `n_vol` 對得上拓樸(例 hexapole halfcut = 9 個 volume)
- `n_elem` 在預期區間(防呆:< 1e3 太粗、> 1e7 跑不完)
- 鋼體積與 CAD 比 < 1% 偏差
- coil 數 = `n_poles`
- `D,ALL,MAG,0` 在腳本內存在(`simulation-constraints.md` §7)

**使用者目視**:

```apdl
/VIEW, 1, 1, 1, 1
/PNUM, MAT, 1
VPLOT                                 ! pole / yoke / air 顏色分群
EPLOT                                 ! 網格密度 vs ESIZE
/PNUM, VOLU, 1
VPLOT                                 ! 編號標示,對照 *SET 區塊註解
```

**判準**:5 項自動全綠 + 使用者口頭批准 VPLOT。

**不通過**:對應退回 `apdl-geom-build.md` 的哪步:
- `n_vol` 錯 → Step 3(幾何建構)
- `n_elem` 太多/少 → Step 5(網格 ESIZE)
- 鋼體積偏 → Step 2(*SET 填錯)或 Step 3
- 缺 BC → Step 4(材料 / BC)補

---

## Pole 配置特有檢項

| pole_config | §1/§4 額外要看 |
|---|---|
| hexapole | alpha = 54.74° 約束;3 對 pole 軸正交;上下層 60° offset |
| quadrupole | 4 極 0/90/180/270° 等角分佈;上下對稱(若雙層) |
| dipole | 2 極軸線共線;極性相反 |
| single | 單極軸對齊預期方向(常見坑:抄樣板沒改 VROTAT 角) |

---

## 立意

reference memory `feedback_verify_before_act` — 動手前先 inspect、dry-run、確認。
避免「腳本噴射」反模式(連跑 fix_v1~v8 然後沒一個對)。
