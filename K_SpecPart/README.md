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

Results on the [Titan23](https://www.eecg.utoronto.ca/~kmurray/titan.html) suite. hMETIS times/cuts are averaged over 50 samplings of 20 runs. K-SpecPart used an hMETIS hint with `eigvecs=2`, `best_solns=5`, `solver_iters=40`, `refine_iters=2`.

### K = 2

| Benchmark | Vertices | Hyperedges | Avg hMETIS time (s) | K-SpecPart time (s) | Avg hMETIS cut | K-SpecPart cut |
|-----------|----------|------------|---------------------|---------------------|----------------|----------------|
| sparcT1_core  | 91976   | 92827   | 9.5  | 32.52  | 982.2  | 977  |
| neuron        | 92290   | 125305  | 6.3  | 33.14  | 245.0  | 244  |
| stereovision  | 94050   | 127085  | 8.3  | 22.19  | 171.0  | 169  |
| des90         | 111221  | 139557  | 10.2 | 18.22  | 376.8  | 374  |
| SLAM_spheric  | 113115  | 142408  | 12.7 | 45.18  | 1061.0 | 1061 |
| cholesky_mc   | 113250  | 144948  | 8.6  | 30.36  | 282.0  | 282  |
| segmentation  | 138295  | 179051  | 13.2 | 47.12  | 120.1  | 120  |
| bitonic_mesh  | 192064  | 235328  | 17.6 | 64.24  | 585.2  | 584  |
| dart          | 202354  | 223301  | 14.4 | 51.22  | 837.0  | 805  |
| openCV        | 217453  | 284108  | 11.7 | 42.17  | 435.4  | 434  |
| stap_qrd      | 240240  | 290123  | 15.9 | 56.31  | 377.4  | 464  |
| minres        | 261359  | 320540  | 18.2 | 82.92  | 207.0  | 207  |
| cholesky_bdti | 266422  | 342688  | 19.8 | 91.33  | 1156.0 | 1136 |
| denoise       | 275638  | 356848  | 26.3 | 95.75  | 496.9  | 418  |
| sparcT2_core  | 300109  | 302663  | 29.7 | 112.31 | 1220.7 | 1188 |
| gsm_switch    | 493260  | 507821  | 27.2 | 110.25 | 4235.3 | 1833 |
| mes_noc       | 547544  | 577664  | 42.3 | 121.92 | 634.6  | 633  |
| LU230         | 574372  | 669477  | 46.6 | 162.32 | 3333.3 | 3363 |
| LU_Network    | 635456  | 726999  | 64.1 | 182.23 | 524.0  | 524  |
| sparcT1_chip2 | 820886  | 821274  | 76.8 | 216.82 | 914.2  | 876  |
| directrf      | 931275  | 1374742 | 85.6 | 223.47 | 602.6  | 515  |
| bitcoin_miner | 1089284 | 1448151 | 89.7 | 369.19 | 1514.1 | 1562 |

### K = 3

| Benchmark | Vertices | Hyperedges | Avg hMETIS time (s) | K-SpecPart time (s) | Avg hMETIS cut | K-SpecPart cut |
|-----------|----------|------------|---------------------|---------------------|----------------|----------------|
| sparcT1_core  | 91976   | 92827   | 12.2 | 60.39  | 2187.9 | 1889 |
| neuron        | 92290   | 125305  | 16.7 | 70.55  | 371.6  | 396  |
| stereovision  | 94050   | 127085  | 14.1 | 64.06  | 332.7  | 336  |
| des90         | 111221  | 139557  | 14.8 | 79.95  | 536.5  | 535  |
| SLAM_spheric  | 113115  | 142408  | 19.2 | 72.23  | 2797.1 | 2720 |
| cholesky_mc   | 113250  | 144948  | 21.3 | 82.79  | 886.5  | 864  |
| segmentation  | 138295  | 179051  | 22.4 | 86.86  | 476.1  | 453  |
| bitonic_mesh  | 192064  | 235328  | 23.7 | 107.24 | 895.0  | 895  |
| dart          | 202354  | 223301  | 29.1 | 112.64 | 1189.9 | 1243 |
| openCV        | 217453  | 284108  | 28.8 | 137.72 | 501.8  | 525  |
| stap_qrd      | 240240  | 290123  | 39.3 | 153.39 | 501.2  | 497  |
| minres        | 261359  | 320540  | 42.4 | 167.89 | 309.0  | 309  |
| cholesky_bdti | 266422  | 342688  | 41.8 | 171.23 | 1769.2 | 1755 |
| denoise       | 275638  | 356848  | 37.9 | 189.67 | 952.8  | 915  |
| sparcT2_core  | 300109  | 302663  | 46.7 | 176.26 | 2827.2 | 2249 |
| gsm_switch    | 493260  | 507821  | 52.1 | 172.92 | 4148.6 | 3694 |
| mes_noc       | 547544  | 577664  | 49.8 | 181.42 | 1164.3 | 1125 |
| LU230         | 574372  | 669477  | 48.6 | 192.18 | 4549.5 | 4548 |
| LU_Network    | 635456  | 726999  | 52.3 | 210.04 | 787.1  | 882  |
| sparcT1_chip2 | 820886  | 821274  | 55.8 | 212.23 | 1453.4 | 1404 |
| directrf      | 931275  | 1374742 | 59.2 | 222.18 | 728.2  | 762  |
| bitcoin_miner | 1089284 | 1448151 | 60.3 | 265.86 | 1944.8 | 1917 |

### K = 4

| Benchmark | Vertices | Hyperedges | Avg hMETIS time (s) | K-SpecPart time (s) | Avg hMETIS cut | K-SpecPart cut |
|-----------|----------|------------|---------------------|---------------------|----------------|----------------|
| sparcT1_core  | 91976   | 92827   | 29.1 | 74.41  | 2532.3 | 2492 |
| neuron        | 92290   | 125305  | 23.7 | 84.64  | 431.5  | 431  |
| stereovision  | 94050   | 127085  | 21.5 | 79.42  | 440.2  | 475  |
| des90         | 111221  | 139557  | 20.8 | 89.25  | 695.5  | 747  |
| SLAM_spheric  | 113115  | 142408  | 22.9 | 97.31  | 3371.1 | 3241 |
| cholesky_mc   | 113250  | 144948  | 21.2 | 101.24 | 982.2  | 984  |
| segmentation  | 138295  | 179051  | 33.2 | 116.52 | 496.3  | 490  |
| bitonic_mesh  | 192064  | 235328  | 37.8 | 130.87 | 1304.4 | 1311 |
| dart          | 202354  | 223301  | 47.3 | 145.51 | 1429.9 | 1401 |
| openCV        | 217453  | 284108  | 52.9 | 192.14 | 526.2  | 522  |
| stap_qrd      | 240240  | 290123  | 55.2 | 207.81 | 714.5  | 674  |
| minres        | 261359  | 320540  | 56.9 | 184.54 | 407.0  | 407  |
| cholesky_bdti | 266422  | 342688  | 59.1 | 174.65 | 1874.4 | 1865 |
| denoise       | 275638  | 356848  | 66.3 | 182.19 | 1172.1 | 1001 |
| sparcT2_core  | 300109  | 302663  | 69.1 | 198.42 | 3523.5 | 3558 |
| gsm_switch    | 493260  | 507821  | 67.2 | 202.65 | 5169.2 | 4404 |
| mes_noc       | 547544  | 577664  | 62.7 | 201.18 | 1314.7 | 1346 |
| LU230         | 574372  | 669477  | 66.9 | 204.73 | 6325.3 | 6310 |
| LU_Network    | 635456  | 726999  | 71.2 | 265.67 | 1495.6 | 1417 |
| sparcT1_chip2 | 820886  | 821274  | 72.8 | 395.13 | 1609.8 | 1601 |
| directrf      | 931275  | 1374742 | 85.6 | 301.45 | 1103.6 | 1092 |
| bitcoin_miner | 1089284 | 1448151 | 89.7 | 312.23 | 2605.4 | 2737 |

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
