# .claude/ — repo 級 Claude Code 設定

**用途**：存放本 repo 的 Claude Code（claude.ai/code）本地設定與 path-scoped 編輯規則，不是模擬產物，一般工作不需更動。

**內容**：
- `settings.local.json` — 本機 Claude Code 設定（權限 allowlist 等）。
- `rules/` — 規則檔（2026-10-03 由 `magnetic_sim/.claude/rules/` 移回本層，session 開在 repo 根時才會自動載入）。

**相關**：repo 總覽見 `../CLAUDE.md`（Quick Triggers 列出各規則觸發片語）與 `../README.md`。
