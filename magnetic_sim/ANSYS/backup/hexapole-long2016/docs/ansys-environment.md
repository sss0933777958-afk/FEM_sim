# ANSYS Environment

## Software
- **Version:** ANSYS Mechanical APDL 2025 R2
- **Executable:** `C:\Program Files\ANSYS2025R2\v252\ansys\bin\winx64\MAPDL.exe`
- **Desktop shortcut:** `Mechanical APDL 2025 R2`
- **License:** Research/Academic (single seat)

## Batch Mode Syntax
```
MAPDL.exe -b -np <cores> -m <memory_MB> -dir <workdir> -j <jobname> -i <input.txt> -o <output.out>
```
| Flag | Description |
|------|-------------|
| `-b` | Batch mode (no GUI) |
| `-np 4` | 4 CPU cores |
| `-m 24000` | 24 GB memory allocation |
| `-dir` | MAPDL process working directory — temp files (`.err`, `.log`) go here |
| `-j` | Jobname prefix (e.g. `coil1` → `coil1.err` instead of `file.err`) |
| `-i` | Input APDL script (use absolute path when `-dir` is set) |
| `-o` | Output log file |

**Note:** `-dir` controls where MAPDL process temp files land. The APDL script's `/CWD` command separately controls where solver outputs (`.rst`, `.db`) are written. Both should point to `results/coilN/`.

## Hardware
- **CPU:** Intel Core i5-14500 (14 cores / 20 threads)
- **RAM:** 32 GB DDR5
- **Storage:** NVMe SSD
- **OS:** Windows 11 Pro

## Notes
- Typical solve time: 30-60 min per coil (SMRT=5 mesh, ~2.4M elements)
- Output per coil: ~1.5 GB (rst + db + supporting files)
- Total for 6 coils: ~10 GB
