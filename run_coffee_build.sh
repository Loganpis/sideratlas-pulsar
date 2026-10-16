#!/bin/sh
# Project-local native CPU COFFEE build. All upstream files retain GPL-2.0.
set -eu
COFFEE_PROJECT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
COFFEE_REV=9986522d08b4d887e5ef9a53e12cecae4e376d22
COFFEE_VENDOR="$COFFEE_PROJECT/third_party"
COFFEE_DEP="$COFFEE_VENDOR/coffee-tools/env"
COFFEE_SOURCE="$COFFEE_VENDOR/coffee-portable"
COFFEE_BUILD="$COFFEE_VENDOR/coffee-build-clang"
if [ ! -x "$COFFEE_DEP/bin/cmake" ]; then
  echo 'Missing project-local dependencies. See research/COFFEE.md for the exact setup command.' >&2
  exit 1
fi
if [ ! -d "$COFFEE_VENDOR/coffee-src/.git" ]; then
  git init "$COFFEE_VENDOR/coffee-src"
  git -C "$COFFEE_VENDOR/coffee-src" remote add origin https://github.com/fizban007/CoffeeGPU.git
  git -C "$COFFEE_VENDOR/coffee-src" fetch --depth 1 origin "$COFFEE_REV"
  git -C "$COFFEE_VENDOR/coffee-src" checkout --detach FETCH_HEAD
fi
if [ "$(git -C "$COFFEE_VENDOR/coffee-src" rev-parse HEAD)" != "$COFFEE_REV" ]; then
  echo 'Upstream checkout does not match pinned revision.' >&2
  exit 1
fi
if [ ! -f "$COFFEE_SOURCE/CMakeLists.txt" ]; then
  mkdir -p "$COFFEE_SOURCE"
  git -C "$COFFEE_VENDOR/coffee-src" archive "$COFFEE_REV" | tar -x -C "$COFFEE_SOURCE"
  (cd "$COFFEE_SOURCE" && patch -p1 < "$COFFEE_VENDOR/coffee-portable.patch")
fi
"$COFFEE_DEP/bin/cmake" -S "$COFFEE_SOURCE" -B "$COFFEE_BUILD" \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_OSX_DEPLOYMENT_TARGET=14.0 \
  -DCMAKE_OSX_SYSROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk \
  -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ \
  -Duse_cuda=0 -Duse_double=1 -Dportable_scalar=ON -Dbuild_tests=OFF \
  "-DCOFFEE_DEPENDENCY_PREFIX=$COFFEE_DEP"
"$COFFEE_DEP/bin/cmake" --build "$COFFEE_BUILD" --target coffee-3d coffee-verification -j 4
