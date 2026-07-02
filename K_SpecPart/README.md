# K-SpecPart

![version](https://img.shields.io/badge/version-1.0.0-blue)
![Julia](https://img.shields.io/badge/Julia-1.6%2B-9558B2)
![license](https://img.shields.io/badge/license-BSD--3--Clause-green)

**K-SpecPart** is a supervised spectral framework for **K-way hypergraph partitioning** (K ≥ 2), written in Julia. It does not partition from scratch: given a hypergraph and an initial *hint* partition (typically from hMETIS), it treats the hint as supervision and iteratively refines it into an equal-or-better solution. K-SpecPart is the multi-way extension of the two-way [SpecPart](https://github.com/TILOS-AI-Institute/HypergraphPartitioning) framework.

On the Titan23 suite, K-SpecPart typically matches or beats the best of many hMETIS runs in cut size.

---

## Table of Contents

- [Features](#features)
- [How It Works](#how-it-works)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Configuration](#configuration)
- [Containers (Docker / Apptainer)](#containers-docker--apptainer)
- [Input / Output Formats](#input--output-formats)
- [Parameters](#parameters)
- [Repository Layout](#repository-layout)
- [Architecture](#architecture)
- [Benchmarks](#benchmarks)
- [Troubleshooting](#troubleshooting)
- [Versioning](#versioning)
- [Citation](#citation)
- [License](#license)

---

## Features

- **Supervised spectral refinement** — improves an existing hint partition rather than starting cold.
- **Two-way and K-way** — bisection (`K = 2`) and true multi-way (`K > 2`) via per-block eigenproblems + Multiclass LDA.
- **No commercial license required** — the ILP partitioner builds against open-source **OR-Tools** by default and auto-detects CPLEX only if present.
- **Parallel** — threaded eigensolves, tree candidate generation, and hMETIS runs, with deterministic reductions.
- **Robust** — B-operator regularization + retry for ill-conditioned eigenproblems, graceful fallbacks when optional external tools are missing or fail, and validity gates that guarantee a balanced, non-degenerate result.
- **Portable** — every external-tool path is set via environment variables; a Docker/Apptainer image bundles the whole toolchain.
- **Simple CLI** — `run_kspecpart.jl` runs the full pipeline from one command.

---

## How It Works

K-SpecPart is *supervised*: a good partition lives in the low-dimensional spectral subspace of the hypergraph Laplacian, and biasing that subspace toward a known-good hint yields embeddings from which an even better partition can be recovered.

One refinement iteration:

```
hint partition
      │
      ▼
generalized eigenvectors  ──►  embedding X  (supervised by the hint)
      │
      ▼
reweight graph by |Xi - Xj|  ──►  build spanning tree (LSST / MST / path)
      │
      ▼
distill exact hyperedge cuts on the tree (LCA + RMQ)
      │
      ▼
sweep / METIS-partition the tree  ──►  many candidate partitions
      │
      ▼
FM-refine each candidate (OpenROAD/TritonPart)
      │
      ▼
overlay-cluster the best candidates + hint  ──►  contracted hypergraph
      │
      ▼
exact ILP (small) or 10× hMETIS (large)  ──►  partition of clusters
      │
      ▼
project back to original vertices + FM-refine  ──►  new partition
```

The new partition becomes the hint for the next iteration. After `refine_iters` iterations a final global overlay combines all per-iteration solutions, and the **global lowest-cut valid partition** seen anywhere is returned.

- **K = 2** (`two_way_spectral_refine`): a single generalized eigenproblem separates the two blocks via a bi-clique supervision term; the tree is cut with a linear sweep + METIS.
- **K > 2** (`k_way_spectral_refine`): one eigenproblem per block (block vs. rest), embeddings concatenated and reduced with **Multiclass LDA**, then the tree is partitioned by recursive bisection.

---

## Requirements

### Julia

- [Julia](https://julialang.org/) **1.6+** (developed and validated on 1.7.2).
- Packages are pinned in `Project.toml` / `Manifest.toml` (`IterativeSolvers`, `Laplacians`, `LinearMaps`, `MultivariateStats`, `Metis`, `Graphs`/`LightGraphs`, `SimpleWeightedGraphs`, `LDLFactorizations`, `Combinatorics`, `Clustering`, `DataStructures`, `Shuffle`, …).

### External binaries

| Tool | Role | Required? |
|------|------|-----------|
| [hMETIS](http://glaros.dtc.umn.edu/gkhome/metis/hmetis/overview) (`hmetis`, optionally `khmetis`) | Overlay solve on the contracted hypergraph; also the usual source of the hint. | **Yes** |
| [METIS](https://github.com/KarypisLab/METIS) (`gpmetis`) | Partitions the reweighted tree during candidate generation. | **Yes** |
| ILP partitioner (`ilp_partitioner/`, links [OR-Tools](https://developers.google.com/optimization) or [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio)) | Exact min-cut on the small contracted hypergraph. Falls back to hMETIS if unavailable. | Recommended |
| [OpenROAD / TritonPart](https://github.com/ABKGroup/TritonPart_OpenROAD) (`openroad`) | FM refinement of each candidate. Failures are caught and the unrefined candidate kept. | Recommended |

Only `gpmetis` and `hmetis` are hard requirements; the rest degrade gracefully.

---

## Installation

### 1. Julia packages

```bash
cd K_SpecPart
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### 2. Build the ILP partitioner

It links OR-Tools by default and **auto-detects CPLEX** (if CMake finds a CPLEX install it compiles the CPLEX backend via `-DUSE_CPLEX`; otherwise it uses OR-Tools). No source edits needed either way.

```bash
# CMAKE_PREFIX_PATH points at your OR-Tools install (dir containing lib/cmake/ortools)
cmake -S ilp_partitioner -B ilp_partitioner/build -DCMAKE_PREFIX_PATH=/path/to/ortools
cmake --build ilp_partitioner/build -j        # produces ilp_partitioner/build/ilp_part
```

> The ILP and OpenROAD binaries are dynamically linked against `libortools.so.9`. K-SpecPart prepends the OR-Tools `lib`/`lib64` directory to `LD_LIBRARY_PATH` at module load (override with `KSPECPART_ORTOOLS_LIB`) so child processes resolve it.

### 3. External tools & threads

Install `gpmetis`, `hmetis`, and (optionally) the OpenROAD refiner, then point K-SpecPart at them via environment variables (see [Configuration](#configuration)). Run Julia multithreaded with `-t auto` (or `JULIA_NUM_THREADS`).

If building the whole toolchain is inconvenient, use the [container](#containers-docker--apptainer) instead.

---

## Quick Start

```bash
# 1. produce a hint with hMETIS (K=4, 2% imbalance)
hmetis design.hgr 4 2 10 1 1 1 0 0          # -> design.hgr.part.4

# 2. refine it with K-SpecPart
julia --project=. -t auto run_kspecpart.jl design.hgr \
    --hint design.hgr.part.4 --k 4 --imb 2 \
    --out design.hgr.kspecpart.part.4
```

Sanity-check the Julia environment without any external binaries:

```bash
julia --project=. smoke_test.jl ilp_partitioner/test/ibm10.hgr
```

---

## Usage

### Command line

```bash
julia --project=. -t auto run_kspecpart.jl <hypergraph.hgr> [options]
```

| Option | Description | Default |
|--------|-------------|---------|
| `--hint FILE` | Initial hint partition (hMETIS `.part.K` format). **Strongly recommended.** | — |
| `--fixed FILE` | Fixed-vertex file. | — |
| `--k K` | Number of blocks. | `2` |
| `--imb PCT` | Per-block imbalance (%). | `2` |
| `--eigvecs N` | Eigenvectors per eigenproblem. | `2` |
| `--refine-iters N` | Refinement iterations. | `2` |
| `--solver-iters N` | LOBPCG iterations. | `40` |
| `--best-solns N` | Candidates kept for overlay clustering. | `5` |
| `--ncycles N` | Random cycles for the spectral sparsifier. | `1` |
| `--seed N` | RNG seed. | `0` |
| `--out FILE` | Write the partition (one block id per line). | — |
| `-V, --version` | Print version and exit. | — |
| `-h, --help` | Show help and exit. | — |

### Julia API

```julia
include("specpart.jl")

refined_partition, cutsize = Main.SpecPart.specpart_run(
    "design.hgr";
    hint_file    = "design.hgr.part.4",
    imb          = 2,
    num_parts    = 4,
    eigvecs      = 2,
    refine_iters = 2,
    solver_iters = 40,
    best_solns   = 5,
    ncycles      = 1,
    seed         = 0,
)
```

For batch experiments see `titan_sweep.jl` (resumable Titan sweep across K/seeds, writing best-cut CSVs via `collect_best_cuts.jl`), `regression.jl`, and `ablation.jl`.

---

## Configuration

External-tool and scratch paths default to locations under the repo directory (`@__DIR__`) and are overridable by environment variables — no source edits required.

| Variable | Purpose | Default |
|----------|---------|---------|
| `KSPECPART_SOURCE_DIR` | Scratch dir for temporary files | repo directory |
| `KSPECPART_METIS_DIR` | Directory holding `metis_script.sh` | repo directory |
| `KSPECPART_HMETIS` | `hmetis` binary (recursive bisection) | built-in fallback |
| `KSPECPART_KHMETIS` | Optional `khmetis` (direct k-way); used for the overlay solve at k>2 when set | `""` (fall back to `hmetis`) |
| `KSPECPART_ILP` | ILP partitioner binary | `ilp_partitioner/build/ilp_part` |
| `KSPECPART_REFINER` | OpenROAD/TritonPart FM refiner (`openroad`) | built-in fallback |
| `GPMETIS` | `gpmetis` binary (used by `metis_script.sh`) | built-in fallback |
| `KSPECPART_ORTOOLS_LIB` | OR-Tools lib dir prepended to `LD_LIBRARY_PATH` (empty disables) | site default |

Behavior knobs:

| Variable | Effect | Default |
|----------|--------|---------|
| `JULIA_NUM_THREADS` / `-t` | Thread count for parallel stages | 1 |
| `KSPECPART_LOG_LEVEL` | `info` or `debug` (verbose per-candidate logs) | `info` |
| `KSPECPART_SERIAL_KERNELS` | `1` forces serial matrix-free kernels for bit-for-bit reproducibility | `0` |
| `KSPECPART_MAX_TREE_CANDIDATES` | Upper bound on tree candidates per iteration (caps cost at large k) | `24` |
| `KSPECPART_LSST_TREES` | Randomized low-stretch tree replicates (small k) | `2` |

---

## Containers (Docker / Apptainer)

The `docker/` directory bundles Julia and the full external toolchain, so K-SpecPart runs anywhere Docker, Podman, or Apptainer/Singularity does. Because hMETIS is a closed binary and the OpenROAD refiner is a specific (unpublished) build, those are supplied out of band; see [`docker/README.md`](docker/README.md) for the full guide.

```bash
# Docker / Podman
cp /path/to/hmetis docker/third_party/hmetis          # closed binary, not committed
docker build -f docker/Dockerfile -t kspecpart .
docker run --rm -v "$PWD/data:/data" kspecpart \
  /data/design.hgr --hint /data/hint.part.4 --k 4 --out /data/out.part.4

# Apptainer / Singularity (a .sif is read-only, so bind a writable scratch dir)
apptainer build --fakeroot kspecpart.sif docker/kspecpart.def
apptainer run --bind "$PWD/data:/data" --env KSPECPART_SOURCE_DIR=/data/scratch \
  kspecpart.sif /data/design.hgr --hint /data/hint.part.4 --k 4 --out /data/out.part.4
```

The image builds the ILP against **OR-Tools** (no commercial license baked in).

---

## Input / Output Formats

**Hypergraph** (hMETIS format). First line: `num_hyperedges num_vertices [fmt]`. Each subsequent line lists the vertices (1-indexed) of one hyperedge; an optional leading weight and trailing vertex-weight lines are controlled by `fmt` (`1`/`11`/`10`). Unit weights by default.

```
3 4
1 2
2 3 4
1 4
```

**Hint / partition** (hMETIS `.part.K`): one **0-indexed** block id per line, in vertex order.

**Fixed vertices**: one line per vertex — a block id to pin it to, or `-1` if free.

**Output**: with `--out`, the refined partition is written as one 0-indexed block id per line.

---

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `hypergraph_file` | Input hypergraph (hMETIS format). | — |
| `hypergraph_fixed_file` | Fixed-vertex file. | `""` |
| `hint_file` | Initial hint partition file. | `""` |
| `imb` | Per-block imbalance (%), real-valued. | `2` |
| `num_parts` | Number of blocks K. | `2` |
| `eigvecs` | Eigenvectors per eigenproblem. | `2` |
| `refine_iters` | Refinement iterations. | `2` |
| `solver_iters` | LOBPCG iterations. | `40` |
| `best_solns` | Candidate solutions kept for overlay clustering. | `5` |
| `ncycles` | Random cycles for the spectral sparsifier. | `1` |
| `seed` | Random seed. | `0` |

---

## Repository Layout

The project is a single Julia module (`SpecPart` in `specpart.jl`) that `include`s the other files; external solvers run as out-of-process binaries.

| Path | Role |
|------|------|
| `specpart.jl` | Module entry point: `specpart_run`, `two_way_spectral_refine`, `k_way_spectral_refine`, shared helpers, incumbent tracking. |
| `run_kspecpart.jl` | Command-line front-end. |
| `definitions.jl` | Core structs, matrix-free spectral operators (`hypl`/`clique`/`bi_clique`), configurable paths. |
| `hypergraph.jl`, `io.jl` | Hypergraph construction and hMETIS-format I/O. |
| `isolate_islands.jl` | Union-find; keep the largest connected component. |
| `graphification.jl` | Hypergraph → weighted graph via random-cycle clique expansion. |
| `embedding.jl` | LOBPCG generalized eigensolver with CMG preconditioner (`solve_eigs`). |
| `dimensionality_reduction.jl` / `projection.jl` | LDA / PCA / random-projection reducers. |
| `tree_partition.jl` | Graph reweighting, tree construction, sweeps, recursive K-way bisection. |
| `cut_distillation.jl`, `rmq/` | Exact tree-cut computation via LCA/RMQ. |
| `degree_aware_prims.jl` | Degree-bounded Prim's MST. |
| `metis.jl` | METIS wrappers. |
| `overlay.jl`, `extract_hypergraph.jl` | Overlay clustering + hypergraph contraction. |
| `optimal_attempt_partitioner.jl` | Exact ILP / parallel hMETIS on the contracted hypergraph. |
| `run_triton_part_refiner.jl` | Generates scripts to call the OpenROAD FM refiner. |
| `golden_evaluator.jl` | Ground-truth cut-size and balance evaluator. |
| `cmg/` | Combinatorial Multigrid preconditioner. |
| `ilp_partitioner/` | C++ ILP partitioner (OR-Tools / CPLEX), built with CMake. |
| `logging_setup.jl` | Structured logger. |
| `smoke_test.jl` | Dependency-free pipeline self-check. |
| `titan_sweep.jl`, `collect_best_cuts.jl`, `regression*.jl`, `ablation.jl` | Batch experiments. |
| `docker/` | Container build (Dockerfile + Apptainer def + guide). |

---

## Architecture

### Core data structures (`definitions.jl`)

- **`__hypergraph__`** — the hypergraph in compressed-sparse form both ways: `eptr`/`eind` (hyperedge → pins) and `vptr`/`vind` (vertex → hyperedges), plus `fixed`, `vwts`, `hwts`.
- **`__pindex__`** — the two supervision sides `(p1, p2)` biasing the eigenproblem.
- **`__cut_profile__`** — per-vertex bookkeeping (prefix sums, free/forced edge cuts, predecessors) that lets the tree sweep evaluate the exact hypergraph cut of every tree edge in O(1) after an O(n) precompute.
- **`__least_common_ancestor__`** — Euler tour + RMQ sparse table for O(1) LCA queries.

### Spectral operators

The generalized eigenproblem is `A x = λ B x`, solved matrix-free with LOBPCG + a Combinatorial Multigrid preconditioner:

- `A` = hyperedge Laplacian (`hypl`, clique-model with `(⌊k/2⌋·⌈k/2⌉)/(k−1)` normalization).
- `B` = `clique + 500·bi_clique + reg·I` (mass + hint supervision + regularization for numerical stability).

### Pipeline stages

1. **Read & preprocess** — `read_hypergraph_file`, `isolate_islands` (largest component), `process_hint`.
2. **Supervised eigenvectors** — `solve_eigs`; K-way solves one problem per block and reduces with `lda`.
3. **Candidate generation** — `tree_partition` reweights the graph by embedding distance, builds low-stretch / MST / path trees, distills exact tree cuts, and produces candidates by tree sweep and by METIS-on-tree. The candidate count is auto-scaled with k.
4. **FM refinement** — each candidate is refined by the OpenROAD/TritonPart refiner (failures are caught).
5. **Overlay & exact solve** — `overlay` contracts vertices never separated by any good solution; `optimal_partitioner` solves the small contracted hypergraph exactly (ILP) or with 10× parallel hMETIS.
6. **Project & repeat** — project back, refine, record; the global best valid partition is returned.

---

## Benchmarks

Latest results on the [Titan23](https://www.eecg.utoronto.ca/~kmurray/titan.html) suite: **22 designs × K ∈ {2,3,4} × 3 seeds = 198 runs**. Each run starts from an **hMETIS hint** (the `ub_factor_2` solutions) and refines it with `imb=2`, `eigvecs=2`, `best_solns=5`, `solver_iters=40`, `refine_iters=2`. Columns: **Hint cut** = hMETIS starting point; **K-SpecPart best** = lowest cut over the 3 seeds; **Improvement** = reduction vs. the hint; **Mean runtime** = mean over 3 seeds on 8 threads.

K-SpecPart never worsens the hint (incumbent tracking guarantees a valid, equal-or-better result). Average improvement grows with K — **+2.4% (K=2), +4.6% (K=3), +6.2% (K=4)** — with per-design gains up to **~29%** (`bitcoin_miner`, K=4) and **~41%** (`gsm_switch`, K=3). Reproduce with [`titan_sweep.jl`](titan_sweep.jl) — see [Reproducing the benchmarks](#reproducing-the-benchmarks).

### K = 2  (improved 10/22 designs, mean +2.4%)

| Benchmark | Hint cut (hMETIS) | K-SpecPart best | Improvement % | Mean runtime (s) |
|-----------|-------------------|-----------------|---------------|------------------|
| LU230 | 3363 | 3363 | 0.00 | 386.98 |
| LU_Network | 524 | 524 | 0.00 | 849.96 |
| SLAM_spheric | 1061 | 1061 | 0.00 | 151.97 |
| bitcoin_miner | 1562 | 1336 | 14.47 | 1246.82 |
| bitonic_mesh | 584 | 583 | 0.17 | 439.91 |
| cholesky_bdti | 1156 | 1156 | 0.00 | 725.74 |
| cholesky_mc | 282 | 282 | 0.00 | 366.98 |
| dart | 805 | 785 | 2.48 | 190.68 |
| denoise | 418 | 418 | 0.00 | 91.92 |
| des90 | 374 | 374 | 0.00 | 411.91 |
| directrf | 515 | 494 | 4.08 | 383.67 |
| gsm_switch | 1833 | 1811 | 1.20 | 585.10 |
| mes_noc | 633 | 632 | 0.16 | 118.46 |
| minres | 207 | 207 | 0.00 | 323.41 |
| neuron | 244 | 244 | 0.00 | 154.20 |
| openCV | 434 | 434 | 0.00 | 327.23 |
| segmentation | 120 | 107 | 10.83 | 218.50 |
| sparcT1_chip2 | 876 | 875 | 0.11 | 1455.71 |
| sparcT1_core | 977 | 976 | 0.10 | 107.65 |
| sparcT2_core | 1188 | 1188 | 0.00 | 210.85 |
| stap_qrd | 464 | 371 | 20.04 | 177.82 |
| stereo_vision | 169 | 169 | 0.00 | 120.02 |

### K = 3  (improved 13/22 designs, mean +4.6%)

| Benchmark | Hint cut (hMETIS) | K-SpecPart best | Improvement % | Mean runtime (s) |
|-----------|-------------------|-----------------|---------------|------------------|
| LU230 | 4548 | 4548 | 0.00 | 555.41 |
| LU_Network | 882 | 785 | 11.00 | 789.15 |
| SLAM_spheric | 2720 | 2689 | 1.14 | 137.15 |
| bitcoin_miner | 1917 | 1695 | 11.58 | 2416.52 |
| bitonic_mesh | 895 | 895 | 0.00 | 469.65 |
| cholesky_bdti | 1701 | 1701 | 0.00 | 568.90 |
| cholesky_mc | 864 | 833 | 3.59 | 101.91 |
| dart | 1243 | 1243 | 0.00 | 137.46 |
| denoise | 915 | 837 | 8.52 | 191.60 |
| des90 | 535 | 534 | 0.19 | 88.19 |
| directrf | 762 | 605 | 20.60 | 689.21 |
| gsm_switch | 3694 | 2194 | 40.61 | 300.18 |
| mes_noc | 1125 | 1125 | 0.00 | 474.21 |
| minres | 309 | 309 | 0.00 | 456.06 |
| neuron | 396 | 395 | 0.25 | 169.46 |
| openCV | 525 | 525 | 0.00 | 376.42 |
| segmentation | 453 | 444 | 1.99 | 238.07 |
| sparcT1_chip2 | 1404 | 1397 | 0.50 | 1193.94 |
| sparcT1_core | 1889 | 1881 | 0.42 | 79.35 |
| sparcT2_core | 2249 | 2249 | 0.00 | 262.27 |
| stap_qrd | 497 | 495 | 0.40 | 316.15 |
| stereo_vision | 336 | 336 | 0.00 | 118.07 |

### K = 4  (improved 19/22 designs, mean +6.2%)

| Benchmark | Hint cut (hMETIS) | K-SpecPart best | Improvement % | Mean runtime (s) |
|-----------|-------------------|-----------------|---------------|------------------|
| LU230 | 6310 | 6267 | 0.68 | 1003.94 |
| LU_Network | 1417 | 1360 | 4.02 | 1153.85 |
| SLAM_spheric | 3241 | 2995 | 7.59 | 251.11 |
| bitcoin_miner | 2737 | 1955 | 28.57 | 5018.48 |
| bitonic_mesh | 1311 | 1306 | 0.38 | 679.47 |
| cholesky_bdti | 1865 | 1865 | 0.00 | 972.04 |
| cholesky_mc | 984 | 978 | 0.61 | 193.71 |
| dart | 1401 | 1391 | 0.71 | 279.18 |
| denoise | 1001 | 813 | 18.78 | 398.37 |
| des90 | 747 | 679 | 9.10 | 182.83 |
| directrf | 1092 | 1092 | 0.00 | 1871.79 |
| gsm_switch | 4404 | 3647 | 17.19 | 580.14 |
| mes_noc | 1346 | 1339 | 0.52 | 1286.94 |
| minres | 407 | 407 | 0.00 | 758.99 |
| neuron | 431 | 418 | 3.02 | 245.42 |
| openCV | 522 | 521 | 0.19 | 614.68 |
| segmentation | 490 | 480 | 2.04 | 422.44 |
| sparcT1_chip2 | 1601 | 1590 | 0.69 | 1459.78 |
| sparcT1_core | 2492 | 2309 | 7.34 | 122.61 |
| sparcT2_core | 3558 | 3136 | 11.86 | 549.10 |
| stap_qrd | 674 | 535 | 20.62 | 436.81 |
| stereo_vision | 475 | 465 | 2.11 | 196.60 |

### Reproducing the benchmarks

The full sweep is driven by [`titan_sweep.jl`](titan_sweep.jl); per-run best cuts are aggregated by [`collect_best_cuts.jl`](collect_best_cuts.jl) (invoked automatically at the end):

```bash
julia -t 8 --project=. titan_sweep.jl > sweep.stdout.log 2>&1
```

Paths and scope are environment-overridable:

| Variable | Meaning | Default |
|----------|---------|---------|
| `SWEEP_BENCH_DIR` | directory of `*.hgr` benchmarks | site path (edit or override) |
| `SWEEP_HINT_DIR` | directory of hint partitions | `hints/Titan23/ub_factor_2` |
| `SWEEP_OUT` | output dir (CSVs + partitions) | `titan_sweep_results/` |
| `SWEEP_DESIGNS` | comma-separated designs to restrict to | all `*.hgr` |
| `SWEEP_KS` | comma-separated K values | `2,3,4` |
| `SWEEP_SEEDS` | comma-separated seeds | `0,1,2` |

Hints follow the layout in [`hints/`](hints/) (`<k>_way/<design>.hgr.specpart.ubfactor.2.part.<k>`). Results stream to `titan_sweep_results/sweep_results.csv` (one row per run) and are summarized into `best_cuts.csv` / `best_cuts_pivot.csv`. The sweep is **resumable** (existing partition files are skipped) and **fault-tolerant** (one failure never aborts the run).

---

## Troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| `matrix is not positive definite` (large k, some KaHyPar instances) | The B operator is only positive semi-definite. | Handled automatically (diagonal regularization + retry + fallback to unsupervised eigenvectors). If persistent, lower `--k` or increase `--imb`. |
| `libortools.so.9: cannot open shared object file` | OR-Tools not on the loader path for the refiner/ILP. | Set `KSPECPART_ORTOOLS_LIB` to your OR-Tools `lib`/`lib64` dir. |
| Refiner aborts / `openroad` error | Missing/incompatible OpenROAD build. | Failures are caught (candidate kept unrefined). Fix the path with `KSPECPART_REFINER` or use the container. |
| Degenerate / empty blocks | Imbalance too loose for the chosen k. | K-SpecPart warns when `imb > 100/k`; tighten `--imb`. |
| Result differs run to run | Parallel float reductions. | Set `KSPECPART_SERIAL_KERNELS=1` and a fixed `--seed` for bit-for-bit output. |

---

## Versioning

The current version is in [`VERSION`](VERSION) (`1.0.0`); changes are tracked in [`CHANGELOG.md`](CHANGELOG.md). K-SpecPart follows [Semantic Versioning](https://semver.org/). Check the running version with `julia run_kspecpart.jl --version`.

---

## Citation

If you use K-SpecPart in academic work, please cite the SpecPart framework:

```bibtex
@inproceedings{specpart2022,
  title     = {SpecPart: A Supervised Spectral Framework for Hypergraph Partitioning Solution Improvement},
  author    = {Bustany, Ismail and Kahng, Andrew B. and Koutis, Ioannis and Pramanik, Bodhisatta and Wang, Zhiang},
  booktitle = {IEEE/ACM International Conference on Computer-Aided Design (ICCAD)},
  year      = {2022}
}
```

(Please also cite the K-SpecPart multi-way extension where applicable — update this entry with its reference.)

---

## License

BSD 3-Clause License — Copyright © 2022 TILOS-AI-Institute. See the repository [`LICENSE`](../LICENSE). Bundled/third-party tools (hMETIS, METIS, OR-Tools, CPLEX, OpenROAD) are licensed separately by their respective owners.
