# Technical Documentation

## Files

| File | Description | Audience |
|------|-------------|----------|
| `troubleshooting.md` | Known APDL pitfalls and solutions. **Read first before any modeling work** | Claude / developer |
| `hexapole-build-workflow.md` | Complete build SOP: Steps 1–7 from tip positioning to solver execution | Claude / developer |
| `hexapole-sphere-geometry.md` | Sphere geometry framework: magic angle derivation, coordinate systems, tip positions | Claude / developer |
| `pole-geometry.md` | Pole dimensions, tilt angles (Hung vs Long2016), orthogonality verification, pole length optimization | Claude / developer |
| `hexapole-simulation-reference.md` | Cross-design simulation specification: constraints, materials, elements, solver, superposition | Claude / developer |

## Reading Order (for new sessions)
1. `troubleshooting.md` — avoid known traps
2. `hexapole-build-workflow.md` — understand the full pipeline
3. `hexapole-sphere-geometry.md` — understand the geometry
4. `pole-geometry.md` — understand the part dimensions
5. `hexapole-simulation-reference.md` — understand the simulation rules

This order is enforced by `.claude/rules/hung-docs.md`.
