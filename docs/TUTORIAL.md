# MPAOP.jl Tutorial (v0.4.0)

> Marine Predators Algorithm (MPA) for Julia — **single- and multi-objective**
> optimisation with **serial / multithreaded / MPI / hybrid MPI+threads** execution.
>
> Original algorithm: Faramarzi A., Heidarinejad M., Mirjalili S., Gandomi A.H.,
> *Marine Predators Algorithm: A nature-inspired metaheuristic*,
> Expert Systems with Applications **152** (2020) 113377.

**One entry point: `MOMPA`.** Set `num_objectives = 1` (the default) for
ordinary minimisation, or `≥ 2` for multi-objective optimisation. There is no
separate single-objective function.

---

## Table of contents

1. [Installation](#1-installation)
2. [Five-minute quick start](#2-five-minute-quick-start)
3. [Single-objective optimisation (`num_objectives = 1`)](#3-single-objective-optimisation-num_objectives--1)
4. [Multi-objective optimisation (`num_objectives ≥ 2`)](#4-multi-objective-optimisation-num_objectives--2)
5. [Parallel execution](#5-parallel-execution)
6. [Advanced features](#6-advanced-features)
7. [Quality indicators](#7-quality-indicators)
8. [Saving and reading history](#8-saving-and-reading-history)
9. [Performance and tuning](#9-performance-and-tuning)
10. [Migrating from older versions](#10-migrating-from-older-versions)
11. [Full keyword reference](#11-full-keyword-reference)
12. [FAQ](#12-faq)
13. [Worked examples](#13-worked-examples)

---

## 1. Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Mojiajun2022/MPAOP")
```

or, for local development:

```julia
using Pkg
Pkg.develop(path="/path/to/MPAOP.jl")
```

There are only seven dependencies (`CSV`, `Dates`, `HDF5`, `LinearAlgebra`,
`MPI`, `Random`, `SpecialFunctions`). The unused v0.2 dependencies
(`Distributions`, `FiniteDiff`, `PositiveFactorizations`, `ThreadPools`,
`StaticArrays`, `Combinatorics`, `Glob`, `Distributed`) have all been dropped,
which cut the first `using MPAOP` from about **21 s to about 5 s**.

Verify the installation:

```julia
using Pkg
Pkg.test("MPAOP")     # 257 tests
```

---

## 2. Five-minute quick start

### Single objective — minimise a function

```julia
using MPAOP

# the objective takes a Vector{Float64} and returns a real number
function fobj(x)
    f1 = abs(x[1] + x[2]) - abs(x[3])
    f2 = x[1] * x[2] * x[3] + 18
    f3 = x[1]^2 * x[2] + 3 * x[3]
    return abs(f1) + abs(f2) + abs(f3)
end

lb = fill(-10.0, 3)      # lower bounds
ub = fill(10.0, 3)       # upper bounds

best_fit, best_pos, curve = MOMPA(
    fobj = fobj,
    lb = lb, ub = ub,
    SearchAgents_no = 48,     # population size
    Max_iter = 200,           # iterations
    disp = true               # print progress
)                             # num_objectives defaults to 1

println("best value      = ", best_fit)
println("best parameters = ", best_pos)

using Plots
plot(curve, xlabel="Iterations", ylabel="Best value", leg=false, framestyle=:box)
```

### Multi objective — find a Pareto front

```julia
using MPAOP

function zdt1(x)                       # ZDT1 benchmark
    f1 = x[1]
    g  = 1.0 + 9.0 / (length(x) - 1) * sum(@view x[2:end])
    return [f1, g * (1 - sqrt(f1 / g))]
end

d = 30
AP, AO, curve = MOMPA(
    fobj = zdt1,
    lb = zeros(d), ub = ones(d),
    SearchAgents_no = 100,
    Max_iter = 200,
    num_objectives = 2,        # required
    disp = false
)

# AP: archived parameters (n × dim)
# AO: archived objectives (n × num_objectives)
using Plots
scatter(AO[:, 1], AO[:, 2], xlabel="f1", ylabel="f2", leg=false)
```

---

## 3. Single-objective optimisation (`num_objectives = 1`)

### 3.1 Signature

```julia
best_fit, best_pos, curve = MOMPA(; fobj, lb, ub, SearchAgents_no, Max_iter,
                                    num_objectives = 1, kwargs...)
```

`num_objectives` defaults to `1`, so it can be left out entirely. `MOMPA`
returns a **different triple** depending on it:

| `num_objectives` | returns |
|------------------|---------|
| `1` | `(best_fitness::Float64, best_position::Vector, convergence_curve::Vector)` |
| `≥ 2` | `(archive_positions::Matrix, archive_objectives::Matrix, convergence_curve::Vector)` |

**Returns (single objective)**

| Value | Meaning |
|-------|---------|
| `best_fit` | best objective value found (**full precision since v0.3**, see [§10](#10-migrating-from-older-versions)) |
| `best_pos` | parameters achieving it, length `length(lb)` |
| `curve` | convergence curve of length `Max_iter` (best-so-far per iteration) |

**Required keywords**

| Keyword | Type | Notes |
|---------|------|-------|
| `fobj` | `Function` | `f(x::Vector{Float64}) -> Real` |
| `lb`, `ub` | `Vector{Float64}` | bounds, must have equal length |
| `SearchAgents_no` | `Int` | population size, typically 20–80; for expensive objectives make it a multiple of your core count |
| `Max_iter` | `Int` | iterations. Rule of thumb: `SearchAgents_no × Max_iter ≥ 10000` |

> ⚠️ Each iteration calls the objective **`2 × SearchAgents_no` times** (the
> original MPA evaluates once after the movement phase and once after the FADs
> phase), so the total budget is `2 × SearchAgents_no × Max_iter` evaluations.

### 3.2 Algorithm behaviour

```julia
MOMPA(fobj=f, lb=lb, ub=ub, SearchAgents_no=50, Max_iter=300,
      p0_optional = [1.0, 2.0, 3.0],   # known starting point (becomes agent 1)
      variant = :nmpa,                 # :standard_mpa (default) or :nmpa
      first_stage_ratio  = 1/3,        # phase switch points
      second_stage_ratio = 2/3,
      FADs0 = 0.2,                     # fish aggregating devices probability
      P0 = 0.5,                        # predator behaviour constant
      Fixbox = true)                   # keep parameters inside [lb, ub]
```

* `p0_optional` — supply a good starting guess to speed up convergence
  considerably. Its length must equal `length(lb)`, otherwise it is ignored
  with a warning.
* `variant`
  * `:standard_mpa` — the original paper;
  * `:nmpa` — improved variant with an adaptive inertia weight
    `w = 2·exp(-(6t/T)²)` and a linear `CF`. It usually converges much deeper
    on high-dimensional unimodal problems (30-D sphere: standard reaches
    `3.3e-6`, nmpa reaches `2.8e-159`).
* Larger `FADs0` → more local search; smaller `P0` → more local search.
* Set `Fixbox = false` when the bounds are only a rough guess and you do not
  want the search constrained by them.

### 3.3 Output control

```julia
MOMPA(...; disp = true,          # print each iteration
           disp_param = true,    # also print the current best parameters
           disp_every = 10,      # print every 10 iterations (new — use it for long runs)
           write_csv_log = true, # append the run to a CSV log
           csv_log_filepath = "fit_log.csv",
           csv_flush_every = 20) # flush every 20 iterations (new, cuts I/O)
```

---

## 4. Multi-objective optimisation (`num_objectives ≥ 2`)

### 4.1 Signature

```julia
AP, AO, curve = MOMPA(; fobj, lb, ub, SearchAgents_no, Max_iter,
                        num_objectives = 2, kwargs...)
```

`num_objectives` must be set explicitly to `2` or more (it defaults to `1`), and
`fobj` must return a **vector of that length**, e.g. `[f1, f2]`. Forgetting the
keyword raises an error naming it, rather than failing obscurely.

| Value | Meaning |
|-------|---------|
| `AP` | archived parameters, `n × dim` |
| `AO` | archived objectives, `n × num_objectives` |
| `curve` | archive size after each iteration |

Under MPI, **non-root ranks return `(nothing, nothing, nothing)`** (unchanged
from v0.2).

### 4.2 Archive and leader selection (the main v0.3 improvement)

```julia
MOMPA(...; archive_size_factor = 1.0,      # capacity = factor × SearchAgents_no
           archive_mode = :pareto,         # new, default
           elite_selection = :crowding,    # new, default
           hv_ref = [1.1, 1.1])            # new: report hypervolume live
```

**`archive_mode`**

| Value | Behaviour |
|-------|-----------|
| `:pareto` (default) | The archive *is* the running **non-dominated set**. A candidate enters only if no member dominates it, and every member it dominates is evicted. On overflow the most crowded solutions are dropped. **The returned front is guaranteed to contain no dominated point.** |
| `:fronts` | v0.2 behaviour: sort `archive ∪ population` into fronts and fill front by front. When the first front is smaller than the archive, **dominated** solutions are kept too. |

**`elite_selection`**

| Value | Behaviour |
|-------|-----------|
| `:crowding` (default) | binary tournament on crowding distance — leaders come from sparse regions, giving a more even spread |
| `:random` | v0.2 behaviour: uniform sampling |

Measured over 15 independent runs (30-D, 100 agents × 200 iterations):

| Configuration | ZDT1 HV | ZDT3 HV | time |
|---------------|---------|---------|------|
| `:pareto` + `:crowding` (default) | **0.86796** | **6.71749** | **0.024 s** |
| `:pareto` + `:random` | 0.86713 | 6.71714 | 0.027 s |
| `:fronts` + `:crowding` | 0.86782 | 6.71734 | 0.032 s |
| `:fronts` + `:random` (≈ v0.2) | 0.86727 | 6.71709 | 0.032 s |

How to read this: `:crowding` beats `:random` consistently (better hypervolume
*and* noticeably better spacing — 0.0094 vs 0.0100 averaged over 15 seeds), and
that is why it is the default. The `:pareto` vs `:fronts` difference in
hypervolume is **within run-to-run noise**; `:pareto` is the default because it
is faster and because it *guarantees* the returned front contains no dominated
solution, not because it scores higher.

> To **reproduce v0.2 results exactly**, pass
> `archive_mode = :fronts, elite_selection = :random`.

### 4.3 Three or more objectives

```julia
AP, AO, _ = MOMPA(fobj = dtlz2, lb = zeros(12), ub = ones(12),
                  SearchAgents_no = 100, Max_iter = 300,
                  num_objectives = 3, disp = false)
```

There is no limit on the number of objectives. Non-dominated sorting uses
ENS-SS, which at `n = 2000, M = 2` is over **400× faster** than the v0.2 sort.

---

## 5. Parallel execution

`parallelism` accepts four values:

| Value | Use when |
|-------|----------|
| `:serial` | the objective is cheap (< 10 μs) — the default |
| `:threads` | single machine, objective ≥ 100 μs |
| `:mpi` | multi-node cluster |
| `:mpi_threads` | **new**: MPI across nodes, threads inside each rank (recommended on HPC) |

### 5.1 Threads

```bash
julia -t auto  your_script.jl      # or export JULIA_NUM_THREADS=8
```

```julia
MOMPA(fobj = heavy_f, lb = lb, ub = ub, SearchAgents_no = 64, Max_iter = 100,
      parallelism = :threads)
```

The population is split into `2 × nthreads` chunks dispatched with
`Threads.@spawn`, which is far more effective than the v0.2
`ThreadPools.@qthreads` loop:

| Version | speed-up on 4 threads |
|---------|----------------------|
| v0.2.1 | 1.27× |
| **v0.3.0** | **3.63×** |

> ⚠️ Your objective must be **thread safe**: no writes to globals, no shared
> mutable scratch buffers.

### 5.2 MPI

```julia
# fit_mpi.jl
using MPAOP, MPI

fobj(x) = sum(abs2, x)

best_fit, best_pos, curve = MOMPA(
    fobj = fobj, lb = fill(-10.0, 30), ub = fill(10.0, 30),
    SearchAgents_no = 64, Max_iter = 200,
    parallelism = :mpi)

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    println("best value = ", best_fit)
end
MPI.Finalize()
```

```bash
mpiexec -n 8 julia --project=. fit_mpi.jl
```

**What changed in v0.3 (the largest single win of this rewrite)**

| | v0.2.1 | v0.3.0 |
|---|--------|--------|
| Collectives per sweep | `N/nprocs` `Scatter`+`Gather` round trips, plus up to 6 serialising `bcast` | **2** `Bcast!` + **2** `Allgatherv!` |
| Data movement | `MPI.bcast` round-trips the matrix through Julia `Serialization`, allocating every time | in-place `Bcast!`, zero allocation |
| Population size | padded up to a multiple of `nprocs` (wasted evaluations) | any size, uneven blocks via `Allgatherv!` |
| `MOMPA(parallelism = :mpi)` | **errors out on MPI.jl 0.20** (`Gather(Ref(...))` is not `isbitstype`) | works |

Multi-objective MPI, measured (64 agents × 40 iterations, ~20 μs per evaluation):

| Ranks | v0.2.1 | v0.3.0 | speed-up |
|-------|--------|--------|----------|
| 1 | 0.663 s | 0.150 s | 4.4× |
| 2 | 0.769 s (**slower than 1 rank**) | 0.077 s | 10.0× |
| 4 | 0.323 s | 0.042 s | 7.7× |

Single-objective strong scaling (v0.3, heavy objective):
1 rank 0.333 s → 2 ranks 0.184 s (1.81×) → 4 ranks 0.094 s (3.54×).

For single-objective runs **every rank returns the same triple**, so you no
longer need an `if rank == 0` guard just to read the result.

### 5.3 Hybrid (recommended for large clusters)

```bash
mpiexec -n 4 julia -t 8 --project=. fit_mpi.jl     # 4 ranks × 8 threads = 32 cores
```

```julia
MOMPA(...; parallelism = :mpi_threads)
```

---

## 6. Advanced features

### 6.1 Reproducible runs

```julia
r1 = MOMPA(fobj=f, lb=lb, ub=ub, SearchAgents_no=30, Max_iter=100, seed=2024)
r2 = MOMPA(fobj=f, lb=lb, ub=ub, SearchAgents_no=30, Max_iter=100, seed=2024)
@assert r1 == r2       # bit-identical
```

You can also pass your own generator: `rng = Xoshiro(123)`.

> **Guarantee**: every stochastic decision is taken on rank 0 and broadcast, so
> for a given `seed` the result is **identical across serial, threads and any
> number of MPI ranks**. That makes debugging and paper reproduction easy.

### 6.2 Batch objectives (vectorised / GPU)

If your objective can evaluate a whole batch at once (matrix maths, GPU kernels,
an external solver with batch mode), `fobj_batch` removes the per-agent call
overhead entirely:

```julia
# X is (nagents × dim); return a vector of length nagents
batch_f(X) = vec(sum(abs2, X, dims=2))

MOMPA(fobj = x -> sum(abs2, x),   # scalar version still required as a fallback
      fobj_batch = batch_f,
      lb = lb, ub = ub, SearchAgents_no = 512, Max_iter = 200)
```

The multi-objective form must return an `(nagents × num_objectives)` matrix.

### 6.3 Early stopping

```julia
MOMPA(...; ftol = 1e-10,      # smallest improvement that counts
           patience = 50,     # stop after 50 iterations without improvement
           max_time = 3600.0) # or after one hour of wall clock (seconds)
```

After an early stop `curve` still has length `Max_iter`, padded with the final
best value so plots stay well behaved.

### 6.4 Callbacks

```julia
history = Float64[]
MOMPA(...; callback = (iter, best_fit, best_pos, curve) -> begin
        push!(history, best_fit)
        iter % 50 == 0 && @info "iteration $iter" best_fit
        return best_fit > 1e-8      # return false or :stop to terminate
    end)
```

The multi-objective callback signature is `(iter, AP, AO, curve)`.

### 6.5 Better initial populations

```julia
MOMPA(...; init = :lhs)        # Latin hypercube sampling, much more even coverage
MOMPA(...; init = :center)     # put agent 1 at the centre of the box
MOMPA(...; opposition = true)  # opposition-based learning: evaluate the mirrored
                               # population once and keep the better of each pair
```

`opposition = true` costs one extra generation of evaluations but usually saves
more than that on high-dimensional problems.

### 6.6 Numerical robustness

An objective returning `NaN` is treated as `Inf`, so it can never poison
`findmin` or a dominance comparison. Functions like this work out of the box:

```julia
f(x) = x[1] < 0 ? NaN : sqrt(x[1]) + x[2]^2
```

---

## 7. Quality indicators

v0.3 ships a set of multi-objective indicators (minimisation assumed, one
solution per row):

```julia
using MPAOP

hv = hypervolume(AO, [1.1, 1.1])    # larger is better; the only Pareto-compliant indicator
i  = igd(AO, true_front)            # inverted generational distance, smaller is better
g  = gd(AO, true_front)             # generational distance (convergence only)
sp = spacing_metric(AO)             # Schott spacing, smaller = more uniform
ms = max_spread(AO)                 # maximum spread, larger = wider front

mask  = pareto_filter(AO)           # non-dominated mask (handy for merging runs)
AO_nd = AO[mask, :]
```

`hypervolume` uses an exact dimension-sweep recursion for 2–3 objectives; for
more objectives pass `method = :mc, samples = 500_000` for a Monte-Carlo
estimate.

Merging several independent runs:

```julia
allO  = vcat((MOMPA(...)[2] for _ in 1:10)...)
final = allO[pareto_filter(allO), :]
```

---

## 8. Saving and reading history

```julia
MOMPA(...; saveHDF = true,
           hdf_filepath = "so_history.h5",
           history_save_interval = 50,   # checkpoint every 50 iterations
           hdf_compress = 4)             # new: deflate level 0–9
```

Reading it back:

```julia
positions, performance, is_mo, archive_prey, archive_obj = ReadMPAHistory("so_history.h5")
curve = ReadMPAConvergence("so_history.h5")   # new
```

HDF5 layout (unchanged, so old files stay readable):

```
MPA/PositionsHistory    nagents × dim × nsnapshots
MPA/FitnessHistory      nagents × nsnapshots               (single objective)
MPA/ObjectivesHistory   nagents × nobj × nsnapshots        (multi objective)
MPA/ArchivePrey         narchive × dim                     (multi objective)
MPA/ArchiveObjectives   narchive × nobj                    (multi objective)
MPA/ConvergenceCurve    niter                              (new)
```

`nsnapshots = 2 × iterations actually run` (one snapshot before and one after
the movement phase).

> 🐛 **v0.2 bug**: `MOMPA(saveHDF = true)` with the default
> `history_save_interval` threw a `BoundsError` at the end of the run (a 2-D
> array was indexed with three indices). Fixed in v0.3 with a regression test.

---

## 9. Performance and tuning

### 9.1 v0.2.1 → v0.3.0, measured (Apple M-series, 12 cores, Julia 1.11)

| Benchmark | v0.2.1 | v0.3.0 | Gain |
|-----------|--------|--------|------|
| Single objective, 30-D sphere, 100 agents × 200 iter | 0.480 s / 262 MiB | **0.010 s / 0.1 MiB** | **48× faster, 2600× less memory** |
| Multi objective ZDT1, 30-D, 100 agents × 200 iter | 0.875 s / 3038 MiB | **0.023 s / 3.5 MiB** | **38× / 870×** |
| Multi objective ZDT3, same | 0.807 s / 3019 MiB | **0.021 s / 3.4 MiB** | **38×** |
| Non-dominated sort `n=2000, M=2` (5 calls) | 0.850 s / 3147 MiB | **0.002 s / 0.5 MiB** | **425×** |
| `update_archive` (20 calls) | 0.030 s / 129 MiB | **0.0002 s / 1.0 MiB** | **~150×** |
| Thread speed-up on 4 threads (heavy objective) | 1.27× | **3.63×** | — |
| Multi objective MPI on 4 ranks | 0.323 s | **0.042 s** | **7.7×** |
| `using MPAOP` load time | ≈ 21 s | **≈ 5 s** | 4× |

Solution quality did **not** regress (15 independent runs):

| Problem | v0.2.1 IGD | v0.3.0 IGD | v0.2.1 HV | v0.3.0 HV |
|---------|-----------|-----------|-----------|-----------|
| ZDT1 | 0.00594 | **0.00579** | 0.86756 | **0.86796** |
| ZDT2 | **0.00594** | 0.00609 | **0.53437** | 0.53434 |
| ZDT3 | 0.00669 | **0.00636** | 6.71711 | **6.71749** |

All differences are within one standard deviation, i.e. quality is equal or
slightly better.

### 9.2 Where the speed comes from

1. **Column-major population layout.** Internally the population is stored as
   `dim × N` (one agent per column), so extracting an agent is a contiguous
   `memcpy` and every update loop is stride-1. The old `Prey[i, :]` built a
   fresh `Vector` per evaluation — 40 000 heap allocations for a
   100 agents × 200 iterations run.
2. **Allocation-free hot loops.** The three-phase movement operator draws its
   `randn` / `rand` / Lévy increments on the fly instead of materialising three
   `N × dim` temporaries per half-iteration; `repeat(lb', N, 1)`, `clamp.()`
   and `copy()` were all replaced with in-place equivalents.
3. **Only the random numbers that are used.** Phase 1 needs no Lévy increments
   and phase 3 needs no uniform draws — skipping them saves tens of thousands
   of RNG calls per iteration.
4. **Type stability.** The old code set `comm = nothing` and later reassigned an
   `MPI.Comm`, forcing dynamic dispatch through the whole main loop. Everything
   now lives in a parametric `EvalCtx{F,B,C}` so the inner loops are statically
   dispatched.
5. **A better sorting algorithm.** The textbook `O(M·n²)` sort with `O(n²)`
   allocations was replaced by ENS-SS (lexicographic pre-sort, then sequential
   front insertion) — 400× faster at `n = 2000`.
6. **Incremental archive.** The old version rebuilt `archive ∪ population` with
   `vcat` inside a loop (quadratic copying) and re-sorted it every iteration.
   The new one inserts in place and truncates once.
7. **MPI collectives.** See [§5.2](#52-mpi).
8. **Trimmed dependencies.** Load time 21 s → 5 s.

### 9.3 Tuning advice

* **Cheap objective** (< 10 μs): stay on `:serial`; parallelism costs more than
  it saves.
* **Medium objective** (10 μs – 1 ms): use `:threads`, and make
  `SearchAgents_no` a multiple of the thread count.
* **Expensive objective** (> 1 ms, e.g. a physical simulation): use `:mpi` or
  `:mpi_threads`. A population that is a multiple of the rank count balances
  best, though it is no longer required.
* **Long runs**: set `disp_every = 50` — `println` itself becomes a bottleneck.
* **Calibrating** `SearchAgents_no` / `Max_iter` / `FADs0` / `P0`:

```julia
using MPAOP, Statistics
F = [MOMPA(fobj=fobj, lb=lb, ub=ub, SearchAgents_no=36, Max_iter=200,
           disp=false, seed=i, FADs0=0.3, P0=0.6)[1] for i in 1:500]
println("mean = ", mean(F), "  std = ", std(F))
```

---

## 10. Migrating from older versions

### 10.1 From v0.3 — `SOMPA` was removed

v0.4 has **one solver**, `MOMPA`. The single-objective function `SOMPA` is
gone; single objective is now `num_objectives = 1`, which is the default:

```julia
# v0.3
best, pos, curve = SOMPA(fobj = f, lb = lb, ub = ub,
                         SearchAgents_no = 40, Max_iter = 200)

# v0.4 -- same keywords, same returned triple
best, pos, curve = MOMPA(fobj = f, lb = lb, ub = ub,
                         SearchAgents_no = 40, Max_iter = 200)
```

Every keyword `SOMPA` accepted (`opposition`, `ftol`, `patience`,
`round_digits`, …) is accepted by `MOMPA` and applies when
`num_objectives == 1`. A one-line `sed` migrates a whole project:

```bash
sed -i '' 's/\bSOMPA(/MOMPA(/g' *.jl
```

The default log/history filenames were unified as well:
`mpa_so_log.csv` / `mompa_log.csv` → `mpa_log.csv`, and
`mpa_so_history.h5` / `mompa_history.h5` → `mpa_history.h5`. Pass
`csv_log_filepath` / `hdf_filepath` explicitly if you depend on the old names.

### 10.2 From v0.2.x

Besides the `SOMPA` removal above, three behaviours differ:

| Item | v0.2 | v0.4 | Restore old behaviour |
|------|------|------|-----------------------|
| single-objective return precision | fitness rounded to 4 digits, position to 8 | **full precision** | `round_digits = 4` |
| multi-objective archive | may contain dominated solutions | strictly non-dominated | `archive_mode = :fronts` |
| multi-objective leader selection | uniform random | crowding-distance tournament | `elite_selection = :random` |

> The precision change matters: the old `round(fit, digits=4)` printed a true
> residual of `1.76e-9` as `0.0`, throwing away real information when fitting
> physical data.

### 10.3 From the v0.1 published on GitHub

The v0.1 positional API **still works** — it forwards to the new engine, so old
scripts get the speed-up for free:

```julia
Top_predator_fit, Top_predator_pos, CV =
    MPA(SearchAgents_no, Max_iterations, p0, lb, ub, narvs, fobj;
        disp=true, Fixbox=true, Write=false, FADs0=0.2, P0=0.5)

# MPI variant
MPA_MPI(SearchAgents_no, Max_iterations, p0, lb, ub, narvs, fobj; disp=true, Write=true)

# confidence intervals (no longer needs FiniteDiff / PositiveFactorizations / Distributions)
ci = confidence_interval(best_pos, fobj, 0.95)
```

`CV` is still a `1 × Max_iter` row matrix, so `plot(CV')` keeps working.

New code should call `MOMPA` directly to get seeding, early stopping,
batch evaluation and hybrid MPI.

---

## 11. Full keyword reference

### 11.1 Always available

| Keyword | Default | Meaning |
|---------|---------|---------|
| `fobj` | required | objective function |
| `lb`, `ub` | required | bounds |
| `SearchAgents_no` | required | population size |
| `Max_iter` | required | iterations |
| `p0_optional` | `[]` | starting point for agent 1 |
| `variant` | `:standard_mpa` | or `:nmpa` |
| `first_stage_ratio` | `1/3` | phase 1 fraction |
| `second_stage_ratio` | `2/3` | phase 2 fraction |
| `FADs0` | `0.2` | FADs probability |
| `P0` | `0.5` | predator constant |
| `parallelism` | `:serial` | `:threads` / `:mpi` / `:mpi_threads` |
| `disp` | `true` | print progress |
| `disp_param` | `false` | also print parameters |
| `Fixbox` | `true` | enforce bounds |
| `write_csv_log` | `false` | write a CSV log |
| `csv_log_filepath` | `"mpa_log.csv"` | log path |
| `saveHDF` | `false` | save history |
| `hdf_filepath` | `"mpa_history.h5"` | HDF5 path |
| `history_save_interval` | `typemax(Int)` | checkpoint period |
| **`fobj_batch`** | `nothing` | vectorised objective |
| **`seed`** | `nothing` | RNG seed |
| **`rng`** | `nothing` | custom RNG |
| **`init`** | `:uniform` | `:lhs` / `:center` |
| **`beta`** | `1.5` | Lévy exponent |
| **`max_time`** | `Inf` | wall-clock budget (seconds) |
| **`callback`** | `nothing` | per-iteration callback |
| **`disp_every`** | `1` | print period |
| **`csv_flush_every`** | `1` | CSV flush period |
| **`hdf_compress`** | `0` | deflate level 0–9 |
| **`nthreads`** | `Threads.nthreads()` | thread count |
| **`chunks_per_thread`** | `2` | chunks per thread (load balancing) |
| **`reuse_buffer`** | `true` | reuse the input buffer (see FAQ) |

(**bold** = new in v0.3)

### 11.2 Single objective only (`num_objectives = 1`)

| Keyword | Default | Meaning |
|---------|---------|---------|
| **`opposition`** | `false` | opposition-based initialisation |
| **`ftol`** | `0.0` | improvement threshold for early stopping |
| **`patience`** | `typemax(Int)` | iterations without improvement before stopping |
| **`round_digits`** | `nothing` | round the result (`nothing` = full precision) |

### 11.3 Multi objective only (`num_objectives ≥ 2`)

| Keyword | Default | Meaning |
|---------|---------|---------|
| `num_objectives` | `1` | number of objectives; set to `≥ 2` for multi-objective |
| `archive_size_factor` | `1.0` | archive capacity factor |
| **`archive_mode`** | `:pareto` | or `:fronts` |
| **`elite_selection`** | `:crowding` | or `:random` |
| **`hv_ref`** | `nothing` | hypervolume reference point |

---

## 12. FAQ

**Q1. Threads give a different answer than serial?**
With the same `seed`, `:serial`, `:threads`, `:mpi` and `:mpi_threads` produce
**identical** results. If they differ, your objective is not thread safe
(usually a mutable global).

**Q2. What is `reuse_buffer` for?**
To reach zero allocation, MPAOP copies each agent into a reused buffer before
calling `fobj`. If your objective **stores the vector it was given**
(`push!(saved, x)`, or `return x`), later iterations will overwrite it. Set
`reuse_buffer = false` to pass a fresh array every time (slightly slower).

**Q3. Why is the multi-objective archive always full / always tiny?**
Capacity is `archive_size_factor × SearchAgents_no`. In `:pareto` mode only
non-dominated solutions are stored, so early on it can hold just a few points —
that is expected. Increase `archive_size_factor` for a larger archive.

**Q4. Does MPI print duplicated output?**
No. All printing happens on rank 0 only.

**Q5. Must `SearchAgents_no` be divisible by the number of ranks?**
Not since v0.3. Uneven blocks are handled by `Allgatherv!`. v0.2 padded the
population, which wasted objective evaluations.

**Q6. How do I handle constraints?**
* Box constraints: use `lb` / `ub` with `Fixbox = true`.
* General constraints: add a penalty term, or return `Inf` for infeasible
  points (handled correctly).

**Q7. My objective is slow — how do I watch progress and checkpoint?**

```julia
MOMPA(...; disp_every = 1,
      callback = (it, bf, bp, cv) -> begin
          it % 10 == 0 && open("best.txt","w") do io
              println(io, bf); println(io, bp)
          end
          true
      end)
```

**Q8. How should I cite this?**

> Faramarzi, A., Heidarinejad, M., Mirjalili, S., & Gandomi, A. H. (2020).
> Marine Predators Algorithm: A nature-inspired metaheuristic.
> *Expert Systems with Applications*, 152, 113377.

---

## 13. Worked examples

### Single objective: initial guess, threads, early stopping, full logging

```julia
using MPAOP, Plots

function residual(p)
    t     = 0:0.1:10
    model = @. p[1] * exp(-p[2] * t) + p[3]
    data  = @. 2.5 * exp(-0.8 * t) + 0.3
    return sum(abs2, model .- data)
end

lb = [0.0, 0.0, -1.0]
ub = [10.0, 5.0, 1.0]

best, pos, curve = MOMPA(
    fobj = residual, lb = lb, ub = ub,
    SearchAgents_no = 40, Max_iter = 500,
    p0_optional = [1.0, 1.0, 0.0],
    variant = :nmpa,
    parallelism = :threads,
    init = :lhs, opposition = true,
    seed = 20250101,
    ftol = 1e-14, patience = 80,
    disp = true, disp_every = 25,
    saveHDF = true, hdf_filepath = "fit.h5", hdf_compress = 4,
    history_save_interval = 100)

println("residual   = ", best)
println("parameters = ", pos)
println("95% confidence intervals:")
display(confidence_interval(pos, residual, 0.95))

plot(curve, yscale=:log10, xlabel="Iteration", ylabel="Residual",
     leg=false, framestyle=:box)
```

### Multi objective: run and score in one go

```julia
using MPAOP, Plots

AP, AO, curve = MOMPA(
    fobj = zdt1, lb = zeros(30), ub = ones(30),
    SearchAgents_no = 120, Max_iter = 400, num_objectives = 2,
    seed = 42, hv_ref = [1.1, 1.1], disp = true, disp_every = 50)

println("archive size = ", size(AP, 1))
println("hypervolume  = ", hypervolume(AO, [1.1, 1.1]))
println("spacing      = ", spacing_metric(AO))
println("max spread   = ", max_spread(AO))

scatter(AO[:,1], AO[:,2], xlabel="f₁", ylabel="f₂",
        markersize=3, leg=false, framestyle=:box)
```
