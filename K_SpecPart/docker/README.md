# Running K-SpecPart in a container

This image bundles Julia and every external solver K-SpecPart shells out to, so
it runs anywhere Docker / Podman / Apptainer does. It is the recommended way to
distribute K-SpecPart given its heavy, partly non-redistributable toolchain.

- **Docker/Podman** users build from `docker/Dockerfile`.
- **Apptainer/Singularity** (HPC) users build from `docker/kspecpart.def`, or
  convert a prebuilt OCI image (see below).
- If you only want to *run* a prebuilt image (no build rights), skip to
  [For recipients](#for-recipients-quickstart).

## What is / isn't in the image

| Component | Source | In the image? |
|-----------|--------|---------------|
| Julia + packages | official tarball + `Pkg.instantiate()` | Yes |
| `gpmetis` (METIS) | built from GKlib + METIS | Yes |
| OR-Tools | official prebuilt C++ release | Yes |
| ILP partitioner | built here against OR-Tools (no CPLEX) | Yes |
| OpenROAD FM refiner | local source snapshot `docker/third_party/openroad-src.tar.gz` (commit `a3f3153c`, submodules included) | Yes |
| `hmetis` / `khmetis` | closed binary | **You supply** under `docker/third_party/` |
| CPLEX | commercial | No (OR-Tools used instead) |

`hmetis` cannot be legally committed, so drop it into `docker/third_party/`
before building (it is git-ignored). If you don't, mount it at runtime and set
`KSPECPART_HMETIS` instead.

**OpenROAD source:** the exact refiner commit (`a3f3153c`, originally the
`master` branch of `TILOS-AI-Institute/HypergraphPartitioning`) was never pushed
to a remote, so it is bundled as a self-contained tarball
(`docker/third_party/openroad-src.tar.gz`, provenance in `openroad-commit.txt` /
`openroad-submodules.txt`). Both are git-ignored. Regenerate the tarball with:

```bash
tar czf docker/third_party/openroad-src.tar.gz \
  --exclude=./build --exclude=./.git --exclude='*/.git' \
  -C /path/to/TritonPart_OpenROAD .
```

If you later publish the fork, build from git instead by swapping the OpenROAD
`RUN` in the Dockerfile for
`bash /tmp/build_openroad.sh --git <repo-url> <commit> /opt/openroad`.

## Build

From the repository root:

```bash
# 1. provide the closed hMETIS binary (required)
cp /path/to/hmetis docker/third_party/hmetis
# cp /path/to/khmetis docker/third_party/khmetis   # optional

# 2. ensure the OpenROAD source snapshot is present (already generated):
ls docker/third_party/openroad-src.tar.gz

# 3. build
docker build -f docker/Dockerfile -t kspecpart .
```

Build args (all optional):

| Arg | Default | Meaning |
|-----|---------|---------|
| `JULIA_VERSION` | `1.7.2` | Julia version |
| `ORTOOLS_VERSION` | `9.4.1874` | OR-Tools release |
| `DEBIAN_VERSION` | `bullseye` | base image (match the OR-Tools/hMETIS ABI) |

The OpenROAD refiner is built from `docker/third_party/openroad-src.tar.gz` (see
"OpenROAD source" above), so no repo/commit build args are needed.

### Apptainer / Singularity (HPC)

`docker/kspecpart.def` mirrors the Dockerfile. Building it needs **fakeroot**
(a configured `subuid`/`subgid` range for your account):

```bash
cp /path/to/hmetis docker/third_party/hmetis
apptainer build --fakeroot kspecpart.sif docker/kspecpart.def
```

If you *don't* have fakeroot/subuid on the cluster (common), build the OCI image
on a Docker-capable host and **convert** it — conversion needs no root:

```bash
# on a Docker host:
docker build -f docker/Dockerfile -t kspecpart . && docker save kspecpart -o kspecpart.tar
# on the cluster:
apptainer build kspecpart.sif docker-archive://kspecpart.tar
```

## Run

Mount a directory with your hypergraph + hint and pass CLI flags (they go
straight to `run_kspecpart.jl`):

```bash
docker run --rm -v "$PWD/data:/data" kspecpart \
  /data/design.hgr --hint /data/design.hgr.part.4 --k 4 --imb 2 \
  --out /data/design.hgr.kspecpart.part.4
```

Use all cores (the entrypoint already passes `-t auto`); cap them with
`docker run --cpus=N`. Increase shared memory for large designs with
`--shm-size`.

For **Apptainer/Singularity**, a `.sif` is read-only, so point the scratch dir
at a writable bind mount:

```bash
apptainer run --bind "$PWD/data:/data" --env KSPECPART_SOURCE_DIR=/data/scratch \
  kspecpart.sif /data/design.hgr --hint /data/hint.part.4 --k 4 --out /data/out.part.4
```

## For recipients (quickstart)

You received the K-SpecPart repository (or a release). To run it:

1. **Get hMETIS** (closed binary, not shipped). Download it from the
   [hMETIS site](http://glaros.dtc.umn.edu/gkhome/metis/hmetis/overview) and
   place it at `docker/third_party/hmetis`.
2. **Get the OpenROAD source snapshot.** If it isn't already at
   `docker/third_party/openroad-src.tar.gz`, download it from the project's
   GitHub **Release assets** and put it there.
3. **Build** with whatever you have:
   - Docker: `docker build -f docker/Dockerfile -t kspecpart .`
   - Podman: `podman build -f docker/Dockerfile -t kspecpart .`
   - Apptainer: `apptainer build --fakeroot kspecpart.sif docker/kspecpart.def`
4. **Run** on your hypergraph + hint (see [Run](#run) above).

Nothing else is host-specific: every tool path is overridable via the
`KSPECPART_*` / `GPMETIS` environment variables, so you can also mount your own
hMETIS/OpenROAD instead of baking them in.

## Maintainer: cutting a release

Because the OpenROAD snapshot (372M) and hMETIS are git-ignored, ship them out
of band:

1. Regenerate the OpenROAD snapshot + provenance if the refiner changed:
   ```bash
   git -C /path/to/TritonPart_OpenROAD rev-parse HEAD > docker/third_party/openroad-commit.txt
   git -C /path/to/TritonPart_OpenROAD submodule status > docker/third_party/openroad-submodules.txt
   tar czf docker/third_party/openroad-src.tar.gz \
     --exclude=./build --exclude=./.git --exclude='*/.git' -C /path/to/TritonPart_OpenROAD .
   ```
2. Tag the K-SpecPart commit and create a GitHub Release.
3. **Attach** `openroad-src.tar.gz` (+ `openroad-commit.txt`,
   `openroad-submodules.txt`) as release assets. Do **not** attach hMETIS or
   CPLEX (not redistributable).
4. In the release notes, record the pinned versions (OpenROAD commit
   `a3f3153c`, OR-Tools `9.4.1874`, Julia `1.7.2`) and link this file.

## Notes

- **32-bit hMETIS:** the bundled hMETIS is a 32-bit ELF, so the image enables
  i386 multiarch and installs `libc6:i386`. If you supply a 64-bit hMETIS
  instead, that dependency is simply unused.
- The OR-Tools prebuilt archive must match the base OS ABI. The default targets
  Debian 11 (`bullseye`); if you change `DEBIAN_VERSION`, pick the matching
  OR-Tools asset.
- The OpenROAD build is the slow, fragile step. If it fails, build OpenROAD once
  on the host and either (a) `COPY` the binary in, or (b) mount it and set
  `KSPECPART_REFINER` at runtime.
- Everything is path-configurable via the same `KSPECPART_*`/`GPMETIS`
  environment variables documented in the top-level README, so you can override
  any bundled tool with a mounted one.
