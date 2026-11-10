# GEMMul8_Julia

This repository provides Julia code for validating the numerical accuracy and error analysis of Ozaki-II scheme.
The implementation is not optimized for performance.

## 0. Tested environment

- Julia Version 1.12.6
- Commit 15346901f0 (2026-04-09 19:20 UTC)
- Build Info:
  - Official https://julialang.org release
- Platform Info:
  - OS: Windows (x86_64-w64-mingw32)
  - CPU: 8 × Intel(R) Core(TM) Ultra 9 288V
  - WORD_SIZE: 64
  - LLVM: libLLVM-18.1.7 (ORCJIT, arrowlake-s)
  - GC: Built with stock GC
- Threads: 1 default, 1 interactive, 1 GC (on 8 virtual cores)

## 1. Set the working directory and verify Julia

Navigate to the `GEMMul8_Julia` directory and check your Julia installation using the commands below.
Replace the example path with the actual location of the directory.

If the `julia` command is not found, install Julia if necessary, add the directory containing the Julia executable to your `PATH`, and restart your terminal.

### Windows PowerShell

```powershell
Set-Location "C:\path\to\GEMMul8_Julia"
julia --version
```

### Linux / macOS

```bash
cd "/path/to/GEMMul8_Julia"
julia --version
```

## 2. Run the experiments

Use the following command on Windows (PowerShell), Linux, macOS, or WSL:

```text
julia --project=. test.jl
julia --project=. select_moduli.jl
```

On completion, `test.jl` creates or overwrites `results.csv` in the directory containing the script.
The current configuration uses the following settings:

| Parameter         | Values                           |
| ----------------- | -------------------------------- |
| Matrix dimensions | `m = 8`, `n = 8`, `k = 8192`     |
| Values of `phi`   | `0, 1, 2, 4, 8, 16`              |
| Number of moduli  | `3:40` (3 through 40, inclusive) |

To change these settings, edit the corresponding values in `test.jl`.

`select_moduli.jl` creates or overwrites `select_moduli.csv` in the directory containing the script.
The CSV file contains the data used in Table 1.

The scripts check for the required Julia packages and automatically install any that are not found.

## 3. Generate the plots

Run the following command to generate plots from `results.csv`:

```text
julia --project=. plot.jl
```

`plot.jl` also checks for its required packages and automatically installs any that are not found.
The generated figures are saved in the directory containing `plot.jl`.

| Output file        | Description                      |
| ------------------ | -------------------------------- |
| `plot_S_error.png` | Comparison of SGEMM errors       |
| `plot_D_error.png` | Comparison of DGEMM errors       |
| `plot_S_bound.png` | Comparison of SGEMM error bounds |
| `plot_D_bound.png` | Comparison of DGEMM error bounds |

### Project environment and dependencies

The `--project=.` option selects the current directory as the active Julia project.
Each script uses `Base.find_package` to check whether its required packages are available.
Packages already discoverable on Julia's load path are not added again.
