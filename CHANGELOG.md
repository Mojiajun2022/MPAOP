# Changelog

## v0.4.0

### Breaking

- **`SOMPA` was removed.** `MOMPA` is the single entry point:
  `num_objectives = 1` (now the default) runs the classic single-objective MPA
  and returns the same `(best_fitness, best_position, convergence_curve)`
  triple `SOMPA` returned; `num_objectives ≥ 2` runs the multi-objective
  archive-based solver. Migrating a v0.3 script is a rename
  (`sed -i '' 's/\bSOMPA(/MOMPA(/g' *.jl`) — every `SOMPA` keyword
  (`opposition`, `ftol`, `patience`, `round_digits`, …) is accepted by `MOMPA`
  and applies when `num_objectives == 1`.
- Default output filenames unified: `mpa_so_log.csv` / `mompa_log.csv` →
  `mpa_log.csv`, and `mpa_so_history.h5` / `mompa_history.h5` →
  `mpa_history.h5`.

### Added

- A vector-valued objective run with the default `num_objectives = 1` now
  raises an `ArgumentError` that names the keyword to set, instead of failing
  obscurely inside the evaluation loop.
- `LICENSE` file (MIT).

### Unchanged

- The v0.1 positional API (`MPA`, `MPA_MPI`, `confidence_interval`) still
  works and forwards to `MOMPA`.
- On-disk HDF5 layout, CSV log format, and the multi-objective return
  convention (non-root MPI ranks return `(nothing, nothing, nothing)`).

## v0.3.0

Complete rewrite of the compute core. **The keyword interface of `SOMPA` and
`MOMPA` is unchanged**; existing scripts run unmodified and faster.

### Performance

Measured on an Apple M-series machine (12 cores) with Julia 1.11:

| Benchmark | v0.2.1 | v0.3.0 |
|-----------|--------|--------|
| Single objective, 30-D sphere, 100 agents × 200 iter | 0.480 s / 262 MiB | 0.010 s / 0.1 MiB |
| Multi objective ZDT1, 30-D, 100 agents × 200 iter | 0.875 s / 3038 MiB | 0.023 s / 3.5 MiB |
| Multi objective ZDT3, same | 0.807 s / 3019 MiB | 0.021 s / 3.4 MiB |
| Non-dominated sort `n=2000, M=2`, 5 calls | 0.850 s / 3147 MiB | 0.002 s / 0.5 MiB |
| `update_archive`, 20 calls | 0.030 s / 129 MiB | 0.0002 s / 1.0 MiB |
| Thread speed-up, 4 threads | 1.27× | 3.63× |
| Multi objective MPI, 4 ranks | 0.323 s | 0.042 s |
| `using MPAOP` | ≈ 21 s | ≈ 5 s |

Quality on ZDT1/ZDT2/ZDT3 (15 independent runs) is unchanged within one
standard deviation, slightly better on ZDT1 and ZDT3.

### Changed (internals)

- Population stored column-major as `dim × nagents`; one agent is one
  contiguous column, so passing it to the objective is a `memcpy` and every
  kernel is stride-1.
- Movement, FADs, box projection, greedy memory and initialisation rewritten as
  allocation-free in-place kernels. Random numbers are drawn where they are
  used instead of materialising three `N × dim` matrices per half-iteration,
  and phases that do not need Lévy or uniform draws no longer generate them.
- Evaluation moved behind a parametric `EvalCtx{F,B,C}`, making the main loops
  type stable (`comm = nothing` previously forced dynamic dispatch throughout).
- Non-dominated sorting replaced by ENS-SS (lexicographic pre-sort + sequential
  front insertion).
- Archive is now incremental: in-place insertion with eviction of dominated
  members and a single crowding-distance truncation, instead of `vcat` in a
  loop plus a full re-sort of `archive ∪ population` each iteration.
- MPI: 2 `Bcast!` + 2 `Allgatherv!` per iteration, on raw `Float64` buffers,
  replacing `nagents/nprocs` `Scatter`/`Gather` round trips and six serialising
  `MPI.bcast` calls.
- Threading: `Threads.@spawn` over balanced chunks with per-chunk buffers,
  replacing `ThreadPools.@qthreads`.
- Dependencies reduced from 15 to 7: `Distributions`, `FiniteDiff`,
  `PositiveFactorizations`, `ThreadPools`, `StaticArrays`, `Combinatorics`,
  `Glob` and `Distributed` are gone.

### Added

- `parallelism = :mpi_threads` — hybrid MPI across ranks, threads within.
- `fobj_batch` — vectorised/GPU objective evaluated on the whole population.
- `seed` / `rng` — bit-reproducible runs, identical across serial, threads and
  any number of MPI ranks.
- `ftol`, `patience`, `max_time` — early stopping.
- `callback` — per-iteration hook that can stop the run.
- `init = :lhs | :center` — Latin-hypercube and centred initialisation.
- `opposition = true` — opposition-based learning at start-up (`SOMPA`).
- `archive_mode = :pareto | :fronts` and `elite_selection = :crowding | :random`
  (`MOMPA`).
- `hv_ref` — live hypervolume reporting (`MOMPA`).
- `disp_every`, `csv_flush_every`, `hdf_compress`, `nthreads`,
  `chunks_per_thread`, `reuse_buffer`, `beta`, `round_digits`.
- Quality indicators: `hypervolume`, `igd`, `gd`, `spacing_metric`,
  `max_spread`, `pareto_filter`.
- `ReadMPAConvergence`, and a `MPA/ConvergenceCurve` dataset in saved history.
- `MPA`, `MPA_MPI`, `confidence_interval` — the v0.1 positional API, forwarding
  to the new engine. `confidence_interval` is reimplemented on `LinearAlgebra`
  and `SpecialFunctions` only.
- Test suite (`test/runtests.jl`, 250 tests) and an MPI correctness test
  (`test/mpi_test.jl`) that checks a 17-agent population on 4 ranks reproduces
  the serial result bit for bit.
- Objectives returning `NaN` are treated as `Inf` instead of corrupting
  `findmin` and dominance comparisons.
- Input validation with clear error messages (`lb`/`ub` length, positive
  population, known `variant` / `parallelism` / `archive_mode`).

### Fixed

- `SOMPA(parallelism = :mpi)` crashed on MPI.jl 0.20 —
  `MPI.Gather(Ref(fit_val), comm; root=0)` raised
  `ArgumentError: Type must be isbitstype`.
- `SOMPA(saveHDF = true)` with the default `history_save_interval` raised
  `BoundsError` at the end of the run: a 2-D fitness history was indexed with
  three indices, and a scalar/vector was passed where a matrix was expected.
- The MPI population was padded up to a multiple of the rank count, wasting
  objective evaluations; uneven splits are now handled directly.
- `MOMPA(variant = :nmpa)` added the current position twice during phase 1
  (`+=` where `SOMPA` uses `=`).
- `SaveMPAHistory` deleted the target file and slept 50 ms before rewriting it
  on every checkpoint.
- `select_elite_from_archive` printed a warning on every call when the archive
  was empty.

### Behaviour changes

| Item | v0.2 | v0.3 | Restore old behaviour |
|------|------|------|-----------------------|
| `SOMPA` return precision | fitness rounded to 4 digits, position to 8 | full precision | `round_digits = 4` |
| MOMPA archive contents | may contain dominated solutions | strictly non-dominated | `archive_mode = :fronts` |
| MOMPA leader selection | uniform random | crowding-distance tournament | `elite_selection = :random` |
| `initialization(n, dim, ub, lb)` (v0.1 4-argument form) | returned `dim × n` | returns `n × dim` | transpose the result |

## v0.2.1

Keyword-based `SOMPA` / `MOMPA` interface, multi-objective support, HDF5
history, MPI and thread execution.

## v0.1

Initial release: positional `MPA` / `MPA_MPI`, single objective only.
