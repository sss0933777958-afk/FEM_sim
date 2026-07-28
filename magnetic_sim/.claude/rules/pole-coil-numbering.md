# 磁極 / coil 編號慣例（強制讀取）

**使用者拍板（2026-07-17，因為舊規則錯誤害出真 bug）**：**APDL coil 編號 → paper 磁極的對應是 per-model 的事實，不是全域通則。**
舊的 `FEM_sim/CLAUDE.md` 寫「Mapping: APDL {1,2,3,4,5,6} = Paper {P1,P3,P6,P5,P2,P4}」當全 repo 通則 —— **那是錯的**（那只是 long2016 deck 的建構順序），它擴散成 `matlab/` 底下 38 個硬寫點、把 hung 的 K̄_I 欄位靜默錯位。

當工作涉及：
- 寫/改任何 `apdl_to_paper_idx`、`F` 矩陣、`exc_sign`/`coil_sign`/`s_sink`
- 新建 model 的 APDL geom/sim deck（決定 coil 建構順序、激發方向）
- 載入某 model 的 `data/coilN/` 做擬合、矩陣、畫圖
- 看到 K̄_I / D̄ 的對角號誌怪怪的

**動手前先讀完此規則。**

對應 memory：`feedback_pole_coil_numbering.md`
相關規則：`charge-model-source-convention.md`（號誌，**跟編號是兩回事**）、`actuator-frame.md`（座標系）、`simulation-constraints.md`（magic angle 幾何）。

---

## 🔒 第 1 層：paper 極命名 = 全域通則（龍飛 notation，不可改）

| paper 極 | P1 | P2 | P3 | P4 | P5 | P6 |
|---|---|---|---|---|---|---|
| 方位角 | 0° | 180° | 120° | 300° | 60° | 240° |
| 層 | 下 | 上 | 下 | 上 | 上 | 下 |

- 下極 = **P1/P3/P6**（0/120/240°）、上極 = **P2/P4/P5**（60/180/300°）。
- 三個 model（long2016 / hung / NTU）**已經一致**，這層沒有爭議。對外討論、圖、PDF 一律只講 paper 名。

## 🔒 第 2 層：APDL coil 編號 = per-model（**這層才是坑**）

coil 編號純粹是 **deck 迴圈的建構順序**，各 model 不同：

| model | deck 的 coil 1..6 方位 | = paper 極 | `apdl_to_paper_idx` |
|---|---|---|---|
| **long2016** | 0, 120, 240, 60, 180, 300°（下極先、上極後） | P1,P3,P6,P5,P2,P4 | `[1,3,6,5,2,4]` |
| **NTU** | 0, 120, 240, 60, 180, 300°（同上） | P1,P3,P6,P5,P2,P4 | `[1,3,6,5,2,4]` |
| **hung** | 0, 180, 120, 300, 60, 240°（**直接照 paper 序**） | P1,P2,P3,P4,P5,P6 | **`[1,2,3,4,5,6]`（identity）** |

- **既有三個 model 的資料一律不動**（改編號 = 51 個 `coilN` 夾改名 + 所有既有圖/PDF/memory 的標籤靜默改義，風險遠大於收益）。
- **`main.m` 的 `apdl_to_paper_idx` 必須如實描述「該 model 的 deck」**。抄別的 model 的 map = 靜默錯位（K̄_I 欄跑到別欄、`ismember(...,[1 3 6])` 的翻號翻到錯的欄）。
- 查證法：看該 model 的 deck 怎麼擺 coil（long2016/NTU 看 `PANG`/`ROTA_ANG`；hung 看 `mesh/_mesh_R700.txt` 的 `fc_ra(k)` + `REAL,fc_pi`）。旁證：擬合出的 `G`（未經 F）**對角占優** ⇒ map 正確。

## 🔒 第 3 層：**之後所有新 model** 的硬性要求

1. **編號一律 identity（coil k = paper P k，hung 模式）** —— deck 直接照 paper 序建，`apdl_to_paper_idx = [1 2 3 4 5 6]`。不要再製造新的建構順序。
2. **電流激發方向必須讓 raw FEM 全 source** —— 每顆極被激發時，B 從**它自己的尖端射出**（見 `charge-model-source-convention.md`）。
   - 在 deck 端就做對（`R,N,1,±TURNS*...` 的號誌），**不要靠後處理翻號**。
   - ⚠ 連帶：這種資料的 loader **不得再全域 negate**（`-B_FEM`）——那是 hung/long2016 的歷史包袱，對全 source 的 raw 會把它翻成全 sink。
3. **map 寫在該 model 的 `mt_constants.m` 當常數宣告**，不要在每支腳本硬寫 literal（現況 38 個硬寫點就是這樣長出來的）。

## ✅ 驗證判準（使用者拍板：**要嘛就是看 K̄_I 的正負號**）

全 source + 正確 map ⇒ K̄_I 必須同時滿足：
- **對角占優**（`argmax|row|` = 1..6，identity）← 這條驗 **map**
- **對角全正**、**off-diag 全負**、**每列和 ≈ 0**（電荷中性，emergent 非強制）← 這條驗 **號誌**
- `K̄(1,1) = 5/6 = 0.8333`（gauge）

`solve_KI_bar_gain.m` 的 K̂ 是**自由 LS 解**，唯一的 gauge 是純量 `5/(6·g11)`（改不了任何相對號誌/結構）→ 此檢驗**非循環論證**。

## 📐 對外顯示一律 P1~P6 排序（使用者拍板 2026-07-26）

給使用者看 / 印**任何逐極矩陣**（裸 `G`、K̄_I、ᴮĤ、D…）一律照 **P1~P6** 排列，**不用 APDL coil 序**。
- `G`（`solve_current` 回傳，6×6）：**列（電荷=Pc_base）本來就是 paper 序**，只需換**欄（激發=coil 序）** → 用 `apdl_to_paper_idx` 的**反排列** `paper2coil`（`for j, paper2coil(apdl_to_paper_idx(j))=j`；long2016/NTU `[1,5,2,6,4,3]`、hung `1:6`），`Gp = G(:,paper2coil)`。換完 R150 的 G 對角占優（~+7~8.4 / off-diag ~−1.4）。
- K̄_I / ᴮĤ / D 等**已經過 F 重排**成 paper 序 → 直接印即可；**只有裸 `G` 要自己換欄**。

## 踩過的坑（2026-07-17，別重犯）

`hung/current_base/main.m:44` 抄了 long2016 的 `[1,3,6,5,2,4]` 去解讀 hung 的 identity 資料，**又**疊上從 long2016 抄來的 `coil_sign=[1 -1 1 -1 -1 1]`（hung raw 全 sink、全域 negate 就到位，根本不需要）→ `model_results_R{700,300}.pdf` 的 K̄_I 有 3 個自激發變負。修正後對角占優、六個全正。
⚠ `ℓ̂`/`ĝ_I`/`𝒞`/`κ` **不受影響**（`fitting` 不碰 F；`gB=6/5·G(1,1)`；`𝒞`/`κ` 取 `svd(kernel·Ĥ)`，欄置換/欄變號 = 右乘置換或 diag(±1) → 奇異值不變）—— **所以這種錯不會讓數字爆掉，只會靜默錯，必須靠 K̄_I 號誌檢驗抓**。

## 觸發片語
- 「coil 編號 / 磁極編號 / apdl_to_paper_idx / map」
- 「K̄_I 對角怎麼是負的」/「coil2 是哪顆極」
- 新 model 建 deck、決定 coil 順序或激發方向時

## 何時不適用
- 純幾何 / mesh / IGES（不涉及 coil 激發與 paper 對應）。
