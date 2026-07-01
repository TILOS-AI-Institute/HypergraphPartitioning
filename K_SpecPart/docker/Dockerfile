# syntax=docker/dockerfile:1
#
# K-SpecPart container image: bundles Julia + all external solvers so the tool
# runs on any machine with Docker.
#
# Redistribution note: hMETIS is a closed binary and CPLEX is commercial, so
# neither is baked into a public image. The ILP partitioner is built here
# against the open-source OR-Tools backend (no CPLEX needed). hMETIS must be
# supplied at build time (see docker/README.md); the OpenROAD FM refiner is
# built from the local source snapshot at docker/third_party/openroad-src.tar.gz
# (its exact commit is not published to any remote).
#
# Build (from the repo root, after placing hmetis under docker/third_party/):
#   docker build -f docker/Dockerfile -t kspecpart .
#
# Run:
#   docker run --rm -v "$PWD/data:/data" kspecpart \
#     /data/design.hgr --hint /data/design.hgr.part.4 --k 4 --out /data/out.part.4

ARG DEBIAN_VERSION=bullseye
ARG JULIA_VERSION=1.7.2
ARG ORTOOLS_VERSION=9.4.1874

FROM debian:${DEBIAN_VERSION}
ENV DEBIAN_FRONTEND=noninteractive
# NOTE: the supplied hMETIS is a 32-bit ELF, so enable i386 multiarch and pull
# in 32-bit libc (provides libc.so.6 + libm.so.6) or hmetis won't run.
RUN dpkg --add-architecture i386 \
    && apt-get update && apt-get install -y --no-install-recommends \
      build-essential cmake git wget ca-certificates pkg-config \
      libgomp1 flex bison tcl-dev libreadline-dev python3 \
      libc6:i386 \
    && rm -rf /var/lib/apt/lists/*

# --- Julia ---------------------------------------------------------------
ARG JULIA_VERSION
RUN wget -qO /tmp/julia.tgz \
      "https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_VERSION%.*}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz" \
    && tar -xzf /tmp/julia.tgz -C /opt && rm /tmp/julia.tgz \
    && ln -s "/opt/julia-${JULIA_VERSION}/bin/julia" /usr/local/bin/julia

# --- OR-Tools (prebuilt C++, provides libortools + cmake config) ---------
ARG ORTOOLS_VERSION
RUN wget -qO /tmp/ortools.tgz \
      "https://github.com/google/or-tools/releases/download/v${ORTOOLS_VERSION%.*}/or-tools_amd64_debian-11_cpp_v${ORTOOLS_VERSION}.tar.gz" \
    && mkdir -p /opt/ortools \
    && tar -xzf /tmp/ortools.tgz -C /opt/ortools --strip-components=1 \
    && rm /tmp/ortools.tgz
ENV LD_LIBRARY_PATH=/opt/ortools/lib

# --- METIS (gpmetis) -----------------------------------------------------
RUN git clone --depth 1 https://github.com/KarypisLab/GKlib /tmp/gklib \
    && cd /tmp/gklib && make config prefix=/opt/metis && make -j"$(nproc)" install \
    && git clone --depth 1 https://github.com/KarypisLab/METIS /tmp/metis \
    && cd /tmp/metis && make config prefix=/opt/metis gklib_path=/opt/metis \
    && make -j"$(nproc)" install \
    && rm -rf /tmp/gklib /tmp/metis

# --- OpenROAD / TritonPart FM refiner ------------------------------------
# Built from the self-contained source snapshot in docker/third_party/
# (commit a3f3153c, recorded in openroad-commit.txt) because that exact refiner
# commit is not published to any remote. Heavy build; delegated to a script so
# it caches independently. To build from a published fork instead, swap the RUN
# for: bash /tmp/build_openroad.sh --git <repo-url> <commit> /opt/openroad
COPY docker/build_openroad.sh /tmp/build_openroad.sh
COPY docker/third_party/openroad-src.tar.gz /tmp/openroad-src.tar.gz
RUN bash /tmp/build_openroad.sh --tarball /tmp/openroad-src.tar.gz /opt/openroad \
    && rm -f /tmp/openroad-src.tar.gz

# --- hMETIS (non-redistributable; supply under docker/third_party/) ------
# Copy just the hMETIS binaries into place (not the whole third_party dir,
# which also holds the large OpenROAD tarball). If you supplied neither, the
# `|| true` keeps the build going and you must mount hMETIS at runtime and set
# KSPECPART_HMETIS.
COPY docker/third_party/ /tmp/third_party/
RUN mkdir -p /opt/hmetis \
    && (cp /tmp/third_party/hmetis  /opt/hmetis/hmetis  2>/dev/null || true) \
    && (cp /tmp/third_party/khmetis /opt/hmetis/khmetis 2>/dev/null || true) \
    && rm -rf /tmp/third_party

# --- K-SpecPart itself ---------------------------------------------------
WORKDIR /opt/kspecpart
COPY . /opt/kspecpart
# Drop the bulky build-context artifacts so they aren't baked into the image
# (the OpenROAD tarball was already consumed above; hMETIS lives in /opt/hmetis).
RUN rm -f /opt/kspecpart/docker/third_party/openroad-src.tar.gz \
          /opt/kspecpart/docker/third_party/hmetis \
          /opt/kspecpart/docker/third_party/khmetis \
    && rm -rf /opt/kspecpart/ilp_partitioner/build
RUN cmake -S ilp_partitioner -B ilp_partitioner/build -DCMAKE_PREFIX_PATH=/opt/ortools \
    && cmake --build ilp_partitioner/build -j"$(nproc)"
RUN julia --project=/opt/kspecpart -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

ENV GPMETIS=/opt/metis/bin/gpmetis \
    KSPECPART_HMETIS=/opt/hmetis/hmetis \
    KSPECPART_ILP=/opt/kspecpart/ilp_partitioner/build/ilp_part \
    KSPECPART_REFINER=/opt/openroad/openroad \
    KSPECPART_ORTOOLS_LIB=/opt/ortools/lib

ENTRYPOINT ["julia", "--project=/opt/kspecpart", "-t", "auto", "/opt/kspecpart/run_kspecpart.jl"]
CMD ["--help"]
