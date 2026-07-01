# Changelog

All notable changes to K-SpecPart are documented here. This project adheres to
[Semantic Versioning](https://semver.org/). The current version is in the
[`VERSION`](VERSION) file.

## [1.0.0] - 2026-06-30

First tagged, distributable release.

### Added
- **Command-line entry point** `run_kspecpart.jl` (`--k`, `--imb`, `--hint`,
  `--fixed`, `--out`, `--seed`, …, plus `--version` and `--help`).
- **Containerized distribution** under `docker/`: a `Dockerfile` and an
  Apptainer/Singularity `kspecpart.def`, bundling Julia, METIS, OR-Tools, the
  OpenROAD FM refiner, and the ILP partitioner. Runs under Docker, Podman, and
  Apptainer/Singularity. See `docker/README.md`.
- **Versioning**: `VERSION` file as the single source of truth, surfaced in the
  startup banner and via `run_kspecpart.jl --version`.
- **Structured logging** with levels (`KSPECPART_LOG_LEVEL=debug`) and a
  timestamped, component-tagged format.
- **Best-of-all incumbent tracking** with validity gates, so the returned
  partition is always balanced, non-degenerate, and the global lowest-cut seen.
- **Candidate auto-scaling for large k** (`KSPECPART_MAX_TREE_CANDIDATES`,
  default 24), bounding the ~90+ candidates at k>=8 to ~16-24 while leaving
  k<=4 behavior unchanged.
- Parallelism across eigensolves, tree candidate generation, and hMETIS runs,
  with deterministic reductions and BLAS thread management.

### Changed
- **ILP partitioner builds without CPLEX.** It now links open-source OR-Tools by
  default and auto-detects CPLEX only when present (`USE_CPLEX` build gate). No
  commercial license is required to build or run.
- External-tool paths are configurable via `KSPECPART_*` / `GPMETIS` environment
  variables (no source edits needed to relocate binaries).
- `imb` / `ub_factor` accept real (fractional) values.

### Fixed
- **"matrix is not positive definite" at large k and some KaHyPar instances.**
  The generalized eigenproblem's B operator is only positive semi-definite;
  fixed with diagonal regularization, a retry loop with increasing
  regularization, a fallback to unsupervised Laplacian eigenvectors,
  empty-supervision-block skipping, and an LDA try/catch.
- **hMETIS imbalance not scaled by k.** The desired k-way imbalance is now
  translated to a valid integer `UBfactor` accounting for hMETIS's
  recursive-bisection compounding (and `k*epsilon` for direct k-way `khmetis`);
  `check_balance` uses the k-way `(100/k + ub)%` capacity.
- Dropped refiner flags (`-balance_constraint`, `-seed`), no-op fixed-vertex
  propagation in overlay contraction, an undefined variable when writing
  fixed-vertex files, a threading race in the K-way refine loop, and masked
  external-refiner failures (`set -o pipefail`).
- `LD_LIBRARY_PATH` handling so the OR-Tools-linked refiner/ILP resolve
  `libortools.so.9` regardless of the launching shell.

### Pinned versions
- Julia 1.7.2 · OR-Tools 9.4.1874 · OpenROAD refiner `a3f3153c` ·
  OpenSTA `2a2b5cd` · abc `95b3543`.

[1.0.0]: https://github.com/TILOS-AI-Institute/HypergraphPartitioning/releases/tag/v1.0.0
