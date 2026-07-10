---
globs: ["magnetic_sim/ANSYS/backup/hexapole-long2016/apdl/*.txt"]
---

# APDL Script Editing Rules

- Mark all changes with `[ADDED]` (new code) or `[MODIFIED]` (changed code) comments in English
- Keep 6 Coil scripts synchronized: only `CURR_ARRAY` values should differ between them
- Always verify `D,ALL,MAG,0` boundary conditions exist before the `/SOLU` block
- Preserve original commented-out code (prefixed with `!****`) unless user asks to remove it
- Use tab indentation matching the original file style
- Do not modify geometric parameters (R_norm, POLE_TIP_R, YOKE dimensions, etc.) without explicit approval
- Do not change element types (SOLID96, SOURC36) or material numbers without approval
- When adding POST1 commands, place them after the existing `/POST1` block at end of file
