#!/usr/bin/env bash
#
# Build the OpenROAD/TritonPart FM refiner used by K-SpecPart and install the
# `openroad` binary to <install-dir>/openroad.
#
# Two source modes:
#
#   Tarball mode (default for this repo):
#       build_openroad.sh --tarball <path/to/openroad-src.tar.gz> <install-dir>
#     Builds from a self-contained source snapshot (submodules included). Used
#     because the exact refiner commit is not published to any remote.
#
#   Git mode:
#       build_openroad.sh --git <repo-url> <commit> <install-dir>
#     Clones a published fork and checks out a pinned commit.
#
set -euo pipefail

MODE="${1:?usage: build_openroad.sh --tarball <tgz> <dest> | --git <repo> <commit> <dest>}"
SRC=/tmp/openroad-src

case "${MODE}" in
  --tarball)
    TARBALL="${2:?path to openroad-src.tar.gz}"
    DEST="${3:-/opt/openroad}"
    mkdir -p "${SRC}"
    tar xzf "${TARBALL}" -C "${SRC}"
    cd "${SRC}"
    ;;
  --git)
    REPO="${2:?repo url}"
    COMMIT="${3:?commit or branch}"
    DEST="${4:-/opt/openroad}"
    git clone "${REPO}" "${SRC}"
    cd "${SRC}"
    git checkout "${COMMIT}"
    git submodule update --init --recursive
    ;;
  *)
    echo "unknown mode: ${MODE}" >&2
    exit 2
    ;;
esac

# OpenROAD ships an installer for its (many) C++ dependencies. Try it, but don't
# fail the whole build if a subset is already satisfied by the base image.
if [ -x ./etc/DependencyInstaller.sh ]; then
  ./etc/DependencyInstaller.sh -all || ./etc/DependencyInstaller.sh || true
fi

cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)" --target openroad

mkdir -p "${DEST}"
cp build/src/openroad "${DEST}/openroad"

cd /
rm -rf "${SRC}"
echo "openroad installed to ${DEST}/openroad"
