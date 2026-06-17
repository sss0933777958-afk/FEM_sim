# IGES_converted/ — 單位轉換後的 IGES（給 ANSYS IGESIN）

`../IGES/` 經單位/旗標轉換後、可直接被 ANSYS `IGESIN` 正確讀入的版本。

## 結構：`<model>/`（與 IGES/ 對應）
```
IGES_converted/
├── kuo_quadrupole/
├── long2016_hexapole_halfcut/
└── zhang_quadrupole/
```

## 檔案類型
- `*.iges` — 轉換後幾何（MKS 公尺、flag 已調好）。

## 規則
- **與 `../IGES/` 一一對應、必須同步**；刪除/改名兩邊一起處理。
- kuo（MKS）轉換用 **flag 2 或重 export**；**不可用 hung 的 `sed ,1.0,6,, → ,1.0,1,,`**（那是 hung 英制專用，會把 mm 讀成 inch）。
- 同步 SOP：`doc/workflows/iges-sync-quick.md`。
