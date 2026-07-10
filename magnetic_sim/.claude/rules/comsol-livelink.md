# COMSOL LiveLink 連線方法（跨設計通用規則）

**這份規則回答一個問題:「怎麼成功連上 COMSOL server 用 MATLAB 跑 `.mph`?」**
任何設計(hung / kuo / hexapole-long2016 / 未來 quadrupole)要跑 COMSOL LiveLink 時,
照這份做,不要重新踩坑。

對應 memory:`feedback_comsol_automation`、`feedback_comsol_multiturn_ac_artifact`
詳細 SOP(kuo 實戰):`magnetic_sim/ANSYS/main/doc/workflows/comsol-livelink.md`
一鍵 launcher + 完整 pitfall:`magnetic_sim/ANSYS/main/comsol/kuo_quadrupole/run_matlab_with_comsol.ps1` + 同層 `README.md`

---

## 🔑 核心連線法:兩個 process,不要用整合式 launcher

成功的關鍵是**把 server 跟 client 拆成兩個 process**:

1. 獨立啟一個 `comsolmphserver.exe`(長駐,聽 port 2036)
2. 另一個 `matlab.exe -batch` 當 client,內部 `mphstart(2036)` 連上去

**絕對不要**用 COMSOL 整合式的 `comsolmphserver matlab` 一鍵啟動器 —— 它在
Windows + MATLAB R2025b 上是壞的(見下方「絕對不要」),拆兩個 process 才乾淨。

### 啟動順序(手動兩步)

```powershell
# 步驟 1:獨立啟 server(-login auto -user 防止卡在 stdin 問帳號)
Start-Process -FilePath "G:\my_workspace\software\COMSOL62\Multiphysics\bin\win64\comsolmphserver.exe" `
    -ArgumentList @('-login','auto','-user','kuo','-port','2036') `
    -RedirectStandardOutput "$env:TEMP\comsol_server.log" `
    -RedirectStandardError  "$env:TEMP\comsol_server.err" `
    -NoNewWindow -PassThru

# 等 port 2036 開始 listen 再進下一步
while (-not (Get-NetTCPConnection -LocalPort 2036 -ErrorAction SilentlyContinue)) { Start-Sleep -Seconds 1 }

# 步驟 2:跑 MATLAB client(單一 token 的 wrapper,見下)
Start-Process -FilePath "C:\Program Files\MATLAB\R2025b\bin\matlab.exe" `
    -ArgumentList @('-batch','run_<task>') `
    -WorkingDirectory "G:\my_workspace\code\FEM_sim\<design>\comsol\<topic>" `
    -NoNewWindow -PassThru
```

### Bootstrap wrapper(每個任務寫一個 `run_<task>.m`)

`-batch` 的引數要保持**單一 token**(避免 PowerShell here-string 換行/引號的坑),
所以真正的指令放進一個 `.m` wrapper:

```matlab
% run_<task>.m
addpath('G:\my_workspace\software\COMSOL62\Multiphysics\mli');   % LiveLink API
addpath('G:\my_workspace\code\FEM_sim\<design>\comsol\<topic>');
cd('G:\my_workspace\code\FEM_sim\<design>\comsol\<topic>');
mphstart(2036);     % <-- 這行就是「連上 server」
<task_script>;      % 真正任務,不加 .m
```

### 一鍵(推薦,自動管 server 生命週期)

已封裝好,直接用:

```powershell
G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\comsol\kuo_quadrupole\run_matlab_with_comsol.ps1 `
    -Script run_<task>            # 跑完自動關 server
# 跨多個任務重用同一個 server 就加 -KeepServer
```

它做的事:port 2036 已開就重用、否則啟 server → 跑 `matlab -batch <Script>` →
**等所有 matlab/MATLAB process 結束** → (沒 `-KeepServer` 就)關 server。
新設計可複製這支 ps1,只改 `-ScriptDir` 預設值。

---

## 路徑表 + 換機器先驗證

本機(2026-06,已驗證存在)的實際路徑:

| 元件 | 路徑 |
|---|---|
| COMSOL 根 | `G:\my_workspace\software\COMSOL62\Multiphysics\` |
| LiveLink API(`mphstart` 等)| `…\Multiphysics\mli\` |
| Server exe | `…\Multiphysics\bin\win64\comsolmphserver.exe` |
| MATLAB | `C:\Program Files\MATLAB\R2025b\bin\matlab.exe` |
| Port | `2036`(預設) |
| Server user | `kuo` |

**換機器 / 換版本前,先驗證這幾條路徑再跑**(否則整段會在第一步就掛):

```powershell
Test-Path "G:\my_workspace\software\COMSOL62\Multiphysics\bin\win64\comsolmphserver.exe"
Test-Path "G:\my_workspace\software\COMSOL62\Multiphysics\mli\mphstart.m"
Test-Path "C:\Program Files\MATLAB\R2025b\bin\matlab.exe"
```

路徑不同就先用 `where.exe matlab` / `Get-ChildItem "…\software\COMSOL*"` 找新位置,
更新上表再跑。COMSOL license 須涵蓋 `LLMATLAB`(LiveLink for MATLAB)。

---

## ❌ 絕對不要

**不要用 COMSOL 整合式 `comsolmphserver matlab` 一鍵啟動器**。它在
Windows + MATLAB R2025b 上有兩個壞點(花過 ~1 小時 debug):

1. `-mlroot "C:\Program Files\..."` 的引號被 PowerShell `-ArgumentList` 吃掉 →
   它在第一個空格處截斷,找不到 MATLAB。
2. 它寫死去找 `bin\win32\matlab.exe`,但 R2025b 沒有這個 legacy 路徑。

→ 解法就是上面的「拆兩個 process + `mphstart(2036)`」,完全繞過它。

---

## 必守的 5 條(不照做會卡住或拿不到結果)

1. **server 啟動加 `-login auto -user <name>`**:預設 `-login info` 會在 stdin 問帳號 →
   batch 下永遠卡住。`-login never` 沒 credentials 會 fail。
2. **任務腳本內用 `diary('out.log'); diary on;` 寫 log**:`matlab -batch` 的 stdout 在
   PowerShell `Start-Process` 下常常是空的(thin launcher 先結束,win64 worker 的輸出不回傳)。
   靠 stdout 會以為「沒輸出 = 失敗」,其實只是沒接到。
3. **wrapper 用單一 `.m` token 當 `-batch` 引數**,不要把多行 MATLAB 塞進 `-ArgumentList`
   (here-string 換行不穩)。
4. **等「所有 matlab/MATLAB process」結束才收尾**:`matlab.exe` 是 thin launcher,真正 worker
   是 `bin\win64\MATLAB.exe`;launcher 會先結束,要等到 worker 也消失(launcher ps1 已處理)。
5. **跨任務重用 server 加 `-KeepServer`**:`comsolmphserver` 預設 client 斷線就退;頻繁啟停慢,
   連續多個任務時保持 server 長駐。

---

## 相關坑(連結,不在此重複)

- COMSOL 從 ANSYS 出來的 STEP/IGES 做 Free Tet 報「邊界處理的內部錯誤」→ 把 `imp1.repairtol`
  從 1e-5 拉到 **1e-3 mm**(memory `feedback_comsol_automation` 末段)。
- Multi-Turn Coil(`CoilType='Numeric'`)自動建 `ccc1`、`mphselectcoords` 回空 → 見 kuo README
  pitfall 5–6。
- 看 kHz 頻率相依別急著說「artifact」,先 `mphinterp` 查 σ(memory
  `feedback_comsol_multiturn_ac_artifact`:導電鋼 σ=7e6 的渦電流是真的)。

完整 pitfall 清單 + Multi-Turn Coil / repairtol 細節在
`magnetic_sim/ANSYS/main/comsol/kuo_quadrupole/README.md` 與 memory,**這份規則只管「怎麼連上 server」的通用骨架**。

---

## 觸發片語(任一即啟動此規則)

- 「跑 COMSOL」/「跑 LiveLink」/「跑 `.mph`」
- 「連 COMSOL server」/「mphserver」/「mphstart」
- 「COMSOL 頻掃」/「COMSOL 抽 I/V」/「COMSOL 場探針」

## 何時不適用

- 純 ANSYS / MATLAB 後處理、不開 COMSOL 的工作
- COMSOL Desktop GUI 手動操作(這份只管 headless LiveLink 自動化)
