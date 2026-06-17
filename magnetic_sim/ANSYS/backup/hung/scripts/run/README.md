# Batch Run Scripts

Shell scripts that wrap MAPDL calls to automate running Coil2-6 simulations sequentially. They **do not contain simulation logic** — that lives in `magnetic_sim/ANSYS/backup/hung/apdl/sim/*.txt` and `magnetic_sim/ANSYS/backup/hung/apdl/postproc/*.txt`. These are just convenience loops that pass the right APDL script to MAPDL for each coil.

## Scripts (current)

| Script | What it does |
|--------|--------------|
| `run_coil2to6_refined_v2.sh` | **Main workflow.** For each Coil2–6: solve (`apdl/sim/MT_Hung_Simulate_Coil${i}_filleted.txt`) then POST1 export (`apdl/postproc/post_export_data_coil${i}.txt`). Clears old `*.rmg` first to prevent stale results. Uses `-m 8000` for memory. |
| `run_post2to6.sh` | **Post-only re-run.** Resumes existing `sim1.db` and just runs POST1 to regenerate `.dat` files. Use this when you change the post-processing script but don't need to re-solve. |

## variants/ — superseded versions

Kept for reference; not part of the current workflow.

| Script | History |
|--------|---------|
| `run_coil2to6.sh` | v1: solve only, no POST1 export. Replaced by `_refined.sh` which added POST1. |
| `run_coil2to6_refined.sh` | v2: added POST1 but used `-m 4000`. Replaced by `_refined_v2.sh` which raised to `-m 8000` + added `rm -f` pre-cleanup. |

## Usage

```bash
# From project root:
bash magnetic_sim/ANSYS/backup/hung/scripts/run/run_coil2to6_refined_v2.sh   # ~40 min total
bash magnetic_sim/ANSYS/backup/hung/scripts/run/run_post2to6.sh              # ~5 min total (POST1 only)
```

Logs go to `magnetic_sim/ANSYS/backup/hung/results/logs/run_*.log`.

## Why these are scripts and not in `apdl/`

APDL (`.txt`) files are the actual simulation code read by MAPDL. These `.sh` files are bash wrappers that call MAPDL repeatedly for each coil — they are build-automation, not simulation.
