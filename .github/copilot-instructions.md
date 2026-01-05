# AI assistant instructions for TB-tbG

Purpose
- Short: help an AI coding agent become productive quickly in this Fortran + MPI TB-tbG codebase.
- This file contains concrete, discoverable patterns, build/run commands, and examples from the repo.

Quick facts
- Language: Modern Fortran (free-form .f90 modules).
- Parallelism: MPI (uses `mpi_f08` and `mpiifx` in Makefile).
- Linear algebra: LAPACK (zheevr) + MKL (Makefile uses `-qmkl`).
- Build: `make` → produces `tbsolver`. Clean with `make clean`.

How to run (examples)
- MPI (recommended):
  - `mpirun -n 4 ./tbsolver` (rank 0 performs I/O and result aggregation).
- Serial (single process):
  - `./tbsolver` (works when running with 1 MPI rank).

Key inputs & examples
- Default input file: `input.in` (set in `constants.f90` as `f_input`).
- Minimal band example (from README):
  ```text
  JOB B
  FPOSCAR cont.vasp
  FEIGEN tb_band.dat
  FKPOINTS KPOINTS.band
  ```
- Default POSCAR filename: `cont.vasp` (VASP POSCAR format, `parsePOSCAR` in `ioutils.f90`).
- KPOINTS:
  - Path mode: first integer `L` then pairs of k-points (read by `readKPOINTS`).
  - Mesh mode: `G` or `M` supported to generate Gamma/Monkhorst-Pack grids.
- EELS example (from README): uses `JOB E`, `OMEGA`, `QPOINT`, `QTAG`.
  - `QTAG D` means q in fractional (`direct`) coordinates, `C` means Cartesian.

Important code locations (quick map)
- `src/constants.f90` — global defaults and `init_constants()` (kinds, defaults: `f_input`, `f_poscar`, `f_eig`, `f_kpoint`, `model_type`, `nomega`, `ncache`).
- `src/parser.f90` — input keyword handling (explicit keywords: JOB, FPOSCAR, FEIGEN, FKPOINTS, FWANNIER, IBAND, NCACHE, EELSMODE, OMEGA, QPOINT, QTAG, EFERMI, TEMPERATURE, DELTA, CALC_IQR, TMPI, TDEBUG, MODEL).
- `src/main.f90` — MPI initialization and job dispatch (`B` for band, `E` for EELS).
- `src/mpi_solver.f90` — MPI orchestration, workload partition (`calculate_workload`), broadcasting, and collection patterns.
- `src/solver.f90` — serial solver wrapper (calls `build_H` then LAPACK `zheevr` via `interface.f90`).
- `src/ioutils.f90` — POSCAR/KPOINTS parsing and result writers (`writeBand`, `writeEELS`).
- `src/eels.f90` — physics of EELS (`calculate_IPF_klist`, `calculate_eels`).

Project conventions and patterns to follow
- Module names map to filenames (e.g., `tbmodel` → `tbmodel.f90`). Add new modules similarly.
- Compilation order controlled by `MODULE_ORDER` in `Makefile`. When adding a new module that others depend on, update `MODULE_ORDER` and/or `Makefile` dependency lines.
- Precision: a project-wide kind `prec` is defined in `constants.f90`. Use it for real/complex declarations.
- Indexing and array shapes (be explicit):
  - Eigenvalues: `eig(nbands, nkpts)`.
  - Wavefunctions: `wavef(nions, nbands, nkpts)` when produced.
  - Build H(k): `Hk(nions,nions)`.
- MPI pattern: rank 0 does heavy I/O and aggregates results; workers send groups back to rank 0. Work division uses `ncache` and `calculate_workload`.
- Default behavior: `FEIGEN` triggers writing bands (`job_band`/`write_band` flags set in parser).

Debugging & profiling
- Turn on timers/logs using input flags:
  - `TMPI` (MPI timer output), `TDEBUG` (debug timers). Defaults are set in `init_constants()`.
  - `NCACHE` controls memory vs communication trade-off (group size per MPI rank).
- LAPACK/BLAS: eigenvalue issues show up via `info` return from `zheevr` (check `compute_eig` logs).
- Use `write(*,...)` messages already present; search for `"[Main]"`, `"[Solver]"`, `"[EELS]"` for runtime checkpoints in logs.

Integration & external deps
- Assumes Intel Fortran + MKL by default. To use another compiler, set `FC`/`FFLAGS` in `Makefile`.
- LAPACK routine `zheevr` and BLAS `zgemm` are declared in `interface.f90`; linking must provide these symbols.

Common pitfalls / gotchas
- POSCAR tag: coordinates may be direct (`D`) or cartesian (`C`) — `parsePOSCAR` handles both but ensure the tag line is present.
- KPOINTS path format: `readKPOINTS` expects pairs (start, end) when L mode is used.
- Work distribution intentionally does not give extra remainder load to rank 0 (see `calculate_workload`).

When in doubt (recommended first steps for an agent)
- Build locally with `make` and run a small, single-rank test: `mpirun -n 1 ./tbsolver` with a minimal `input.in` and `cont.vasp`.
- Inspect logs for `"[Main]"`, `"[Solver]"`, `"[EELS]"` to verify the stage reached.
- To add instrumentation, prefer toggles via `constants` flags (`timer_debug`, `timer_mpi`, etc.) and new `write` statements.

If anything above is missing or unclear, ask the maintainers which POSCAR/KPOINTS examples are canonical and whether CI/build variations (non-Intel toolchains) must be supported.