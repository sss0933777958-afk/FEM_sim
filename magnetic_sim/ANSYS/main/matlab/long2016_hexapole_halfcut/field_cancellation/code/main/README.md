# …/field_cancellation/code/main/ — driver

**用途**：source/sink 組合掃描主程式。

**內容**：`sweep_field_cancellation.m` —
1. 載 coil1..6 standard `'all'`（同一 mesh，6 個單線圈 1A 解）。
2. 建兩個節點集（WP 框 z_wp=z−SPH_OFST）：**評分球** `x²+y²+z_wp²≤(50µm)²`（348 節點）、**繪圖切面** `|z_wp|<25µm & x²+y²≤R_norm²(500µm)`（13050 節點）。皆真實節點。
3. 掃 64 種 `s∈{−1,+1}⁶`：`B=Σ s_j·[bx,by,bz]_j`（raw、符號=電流方向）；
   **評分＝R=50µm 球內 mean|B|**（min/max 佐證）。
4. 依 **mean|B|** 升冪排序、印 32 代表（s 與 −s 簡併）；取 min cancel / max / all-source(`[−+−++−]`)。
5. 呼叫 `../plot/plot_field_xy` 在 **500µm 切面**畫該三組合（title 顯示 R=50µm 球 mean|B|）→ 3 張圖。

**參數**：`R_avg=50µm`（評分球）、`R_ws=R_norm(500µm)`（切面範圍）、`zslab=25µm`。

**判準說明**：評分區域＝使用者關注的 **R=50µm 操作球**；mean|B|（非 min）避免單點誤判。評分區(50µm 球) 與繪圖區(500µm 切面) 刻意不同（小球量抵銷、切面看脈絡）。

**資料來源 / 流向**：`ansys_path`→`import_ansys_data`→疊加→`plot_field_xy`→`../../figures/`。表格印 console、不存 `.mat`。

**相關**：見上層 `../README.md`、`../plot/README.md`。
