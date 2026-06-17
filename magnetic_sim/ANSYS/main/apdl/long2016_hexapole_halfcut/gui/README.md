# apdl/long2016_hexapole_halfcut/gui/ — GUI 檢視用建模腳本（非求解）

**用途**：產生「**只建模、不求解**」的 `.db`，供使用者在 MAPDL GUI 開來**檢視整個鐵件幾何 + 每顆線圈的電流繞向**。與 `sim/`（求解）區別：不解（`/SOLU` 全移除），故快且不留求解副產物。

**共同做法**（兩支 deck 都源自 `sim/baseline/MT_Sim_P1.txt` halfcut，只差 `CURR_ARRAY`）：建幾何（9 volumes 半切鐵件）+ 6 SOURC36 線圈後，**不 mesh**（vol 8/yoke 在此幾何 tet-mesh 會 choke，顯示也不需 FE mesh）；用幾何畫 6 條紅色「繞向弧+箭頭」（每顆 ~3/4 圈弧 + 箭頭，**方向依該線圈 `CURR_ARRAY` 正負**：+1 CCW / −1 CW），存成 component `ARROW_LN` 染紅；`/GTYPE VOLU+LINE`+`/SHADE` 讓鐵件以 solid-model 著色實體 + 紅箭頭疊上；`SAVE` 出 `.db`。**改 `CURR_ARRAY` 即可做任何組合**（箭頭方向自動跟著）。

**內容（各組合一支 deck → 一個 db）**：
- `MT_AllSource.txt` — `CURR_ARRAY=[−1,−1,−1,+1,+1,+1]`＝**all-source**（下極線圈 CW、上極 CCW，每極尖端皆射出；場在工作空間抵銷）。→ `db/allsource/allsource.db`，jobname `allsource`。
- `MT_CoilsSameDir.txt` — `CURR_ARRAY=[−1,−1,−1,−1,−1,−1]`＝**6 線圈電流同向**（全 CW；工作空間場最強，即 field_cancellation sweep 的 max 端）。→ `db/coils_same_dir/coils_same_dir.db`，jobname `coils_same_dir`。

**GUI 開法**（在指令列貼上，重現驗證過的視圖；`<job>` = `allsource` 或 `coils_same_dir`，於對應 `db/<job>/` 開）：
```
RESUME,<job>,db
/GRAPHICS,POWER $ /SHADE,,1
VSEL,S,VOLU,,1,6 $ VSEL,A,VOLU,,8      ! 只顯示鐵件(藏空氣 7/9)
CMSEL,S,ARROW_LN                        ! 紅色繞向箭頭線
/GTYPE,1,VOLU,1 $ /GTYPE,1,LINE,1 $ /GTYPE,1,ELEM,0
/COLOR,LINE,RED                          ! 箭頭染紅(若已紅可略)
/VIEW,1,1,-0.7,0.55 $ /ANG,1 $ /AUTO,1
GPLOT                                     ! 鐵件實體 + 6 紅繞向箭頭
```
> 互動 GUI session 會自行產生 `.lock/.log/.err` 等暫存（正常，關閉後可自清）。

**命名 / 慣例**：`MT_*`；改動標 `[MODIFIED]`；不求解故無 `/SOLU`。

**相關**：見 `../README.md`、`../sim/baseline/`（求解版源頭）。
