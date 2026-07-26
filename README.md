# MPAOP.jl

**Marine Predators Algorithm (MPA) for Julia** — single- and multi-objective
global optimisation, running serially, on threads, on MPI, or on both at once.

The algorithm was proposed by Faramarzi *et al.*
(*Expert Systems with Applications* **152** (2020) 113377).
This package started as a port of the reference MATLAB implementation and has
since grown a multi-objective solver, distributed execution, and a
zero-allocation core.

📖 **[Full tutorial → `docs/TUTORIAL.md`](docs/TUTORIAL.md)**

---

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Mojiajun2022/MPAOP")
```

---

## Quick start

### Single objective

```julia
using MPAOP

fobj(x) = abs(abs(x[1] + x[2]) - abs(x[3])) +
          abs(x[1] * x[2] * x[3] + 18) +
          abs(x[1]^2 * x[2] + 3 * x[3])

best_fit, best_pos, curve = SOMPA(
    fobj = fobj,
    lb = fill(-10.0, 3),
    ub = fill(10.0, 3),
    SearchAgents_no = 48,
    Max_iter = 200)

println("best value      = ", best_fit)
println("best parameters = ", best_pos)
```

### Multi objective

```julia
using MPAOP

function zdt1(x)
    f1 = x[1]
    g  = 1.0 + 9.0 / (length(x) - 1) * sum(@view x[2:end])
    return [f1, g * (1 - sqrt(f1 / g))]
end

AP, AO, curve = MOMPA(
    fobj = zdt1,
    lb = zeros(30), ub = ones(30),
    SearchAgents_no = 100,
    Max_iter = 200,
    num_objectives = 2)

println("hypervolume = ", hypervolume(AO, [1.1, 1.1]))
```

### Parallel

```bash
julia -t auto script.jl                 # threads
mpiexec -n 8 julia --project=. script.jl   # MPI
mpiexec -n 4 julia -t 8 --project=. script.jl  # hybrid
```

```julia
SOMPA(...; parallelism = :threads)       # :serial | :threads | :mpi | :mpi_threads
```

---

## What's new in v0.3.0

v0.3 is a full rewrite of the compute core. **The keyword interface is
unchanged** — existing scripts keep working and simply run faster.

### Performance (Apple M-series, 12 cores, Julia 1.11)

| Benchmark | v0.2.1 | v0.3.0 | Gain |
|-----------|--------|--------|------|
| Single objective, 30-D sphere, 100 agents × 200 iter | 0.480 s / 262 MiB | **0.010 s / 0.1 MiB** | **48× faster, 2600× less memory** |
| Multi objective ZDT1, 30-D, 100 agents × 200 iter | 0.875 s / 3038 MiB | **0.023 s / 3.5 MiB** | **38× / 870×** |
| Non-dominated sort, `n=2000, M=2` | 0.850 s / 3147 MiB | **0.002 s / 0.5 MiB** | **425×** |
| `update_archive`, 20 calls | 0.030 s / 129 MiB | **0.0002 s / 1.0 MiB** | **~150×** |
| Thread speed-up (4 threads, heavy objective) | 1.27× | **3.63×** | — |
| Multi objective MPI, 4 ranks | 0.323 s | **0.042 s** | **7.7×** |
| `using MPAOP` load time | ≈ 21 s | **≈ 5 s** | 4× |

Solution quality is unchanged or slightly better (15 independent runs, ZDT1/2/3
IGD and hypervolume all within one standard deviation of v0.2).

### How

* **Column-major population layout** (`dim × nagents`) — one agent per
  contiguous column, so extracting an agent is a `memcpy` and every kernel is
  stride-1.
* **Allocation-free hot loops** — random draws are generated on the fly instead
  of materialising three `N × dim` matrices per half-iteration; `repeat`,
  `clamp.` and `copy` replaced by in-place kernels.
* **Type-stable core** — the evaluation context is parametric
  (`EvalCtx{F,B,C}`), removing the dynamic dispatch caused by `comm = nothing`.
* **ENS-SS non-dominated sorting** instead of the textbook `O(M·n²)` sort with
  quadratic allocation.
* **Incremental Pareto archive** — in-place insertion and a single crowding
  based truncation, instead of `vcat`-in-a-loop plus a full re-sort.
* **Real MPI collectives** — 2 `Bcast!` + 2 `Allgatherv!` per iteration instead
  of `nagents/nprocs` `Scatter`/`Gather` round trips and six serialising
  `MPI.bcast` calls.
* **Seven dependencies instead of fifteen.**

### New features

| Feature | Keyword |
|---------|---------|
| Hybrid MPI + threads | `parallelism = :mpi_threads` |
| Vectorised / GPU objectives | `fobj_batch = f` |
| Bit-reproducible runs (identical across serial / threads / any rank count) | `seed`, `rng` |
| Early stopping | `ftol`, `patience`, `max_time` |
| Per-iteration callback | `callback` |
| Latin-hypercube / centred initialisation | `init = :lhs` / `:center` |
| Opposition-based learning | `opposition = true` |
| Strictly non-dominated archive | `archive_mode = :pareto` (default) |
| Crowding-distance leader selection | `elite_selection = :crowding` (default) |
| Live hypervolume reporting | `hv_ref = [...]` |
| HDF5 compression | `hdf_compress = 0…9` |
| Throttled printing / logging | `disp_every`, `csv_flush_every` |
| Quality indicators | `hypervolume`, `igd`, `gd`, `spacing_metric`, `max_spread`, `pareto_filter` |

### Bug fixes

* `SOMPA(parallelism = :mpi)` **crashed** on MPI.jl 0.20
  (`MPI.Gather(Ref(x), …)` → `ArgumentError: Type must be isbitstype`).
* `SOMPA(saveHDF = true)` with the default `history_save_interval` threw a
  `BoundsError` at the end of every run (a 2-D array indexed with three
  indices), and passed a scalar where a matrix was expected.
* The MPI population was silently padded up to a multiple of the rank count,
  wasting objective evaluations. Uneven splits are now handled by
  `Allgatherv!`.
* `MOMPA(variant = :nmpa)` double-counted the current position in phase 1
  (`+=` where `SOMPA` uses `=`).
* `SaveMPAHistory` deleted the output file and then `sleep(0.05)`-ed before
  rewriting it — a 50 ms stall on every checkpoint, for no benefit.
* The single-objective result was silently rounded to 4 digits, printing a
  genuine residual of `1.76e-9` as `0.0`. Full precision is returned now; pass
  `round_digits = 4` for the old behaviour.

### Behaviour changes

Three defaults changed. Each has a keyword that restores the v0.2 behaviour:

| Item | v0.2 | v0.3 | Restore |
|------|------|------|---------|
| `SOMPA` return precision | rounded (4 / 8 digits) | full precision | `round_digits = 4` |
| MOMPA archive | may contain dominated solutions | strictly non-dominated | `archive_mode = :fronts` |
| MOMPA leader selection | uniform random | crowding tournament | `elite_selection = :random` |

---

## Compatibility

The positional API published in v0.1 still works; it forwards to the new engine:

```julia
Top_predator_fit, Top_predator_pos, CV =
    MPA(SearchAgents_no, Max_iterations, p0, lb, ub, narvs, fobj;
        disp=true, Fixbox=true, Write=false, FADs0=0.2, P0=0.5)

MPA_MPI(SearchAgents_no, Max_iterations, p0, lb, ub, narvs, fobj; disp=true)

ci = confidence_interval(best_pos, fobj, 0.95)
```

`CV` keeps its `1 × Max_iter` shape, so `plot(CV')` still works.
`confidence_interval` was reimplemented on `LinearAlgebra` +
`SpecialFunctions`, so it no longer drags in `FiniteDiff`,
`PositiveFactorizations` and `Distributions`.

---

## Testing

```julia
using Pkg; Pkg.test("MPAOP")           # 250 tests
```

```bash
mpiexec -n 4 julia --project=. test/mpi_test.jl    # MPI correctness
```

The MPI test checks that a 17-agent population on 4 ranks reproduces the serial
result **bit for bit**.

---

## Citation

> Faramarzi, A., Heidarinejad, M., Mirjalili, S., & Gandomi, A. H. (2020).
> Marine Predators Algorithm: A nature-inspired metaheuristic.
> *Expert Systems with Applications*, 152, 113377.
