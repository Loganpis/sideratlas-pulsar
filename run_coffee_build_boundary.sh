#!/bin/sh
# Separate experimental source/binary: never replaces the upstream-boundary build.
set -eu
COFFEE_PROJECT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
COFFEE_DEP="$COFFEE_PROJECT/third_party/coffee-tools/env"
COFFEE_BASE="$COFFEE_PROJECT/third_party/coffee-portable"
COFFEE_VARIANT="$COFFEE_PROJECT/third_party/coffee-surface-anchored"
if [ ! -f "$COFFEE_BASE/CMakeLists.txt" ]; then
  "$COFFEE_PROJECT/run_coffee_build.sh"
fi
if [ ! -f "$COFFEE_VARIANT/CMakeLists.txt" ]; then
  python3 - "$COFFEE_BASE" "$COFFEE_VARIANT" <<'PY'
import shutil,sys
shutil.copytree(sys.argv[1],sys.argv[2],ignore=shutil.ignore_patterns('bin','lib'))
PY
  (cd "$COFFEE_VARIANT" && patch -p1 < "$COFFEE_PROJECT/third_party/coffee-surface-anchored.patch")
fi
"$COFFEE_DEP/bin/cmake" -S "$COFFEE_VARIANT" \
  -B "$COFFEE_PROJECT/third_party/coffee-build-surface-anchored" \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_OSX_DEPLOYMENT_TARGET=14.0 \
  -DCMAKE_OSX_SYSROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk \
  -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ \
  -Duse_cuda=0 -Duse_double=1 -Dportable_scalar=ON -Dbuild_tests=OFF \
  "-DCOFFEE_DEPENDENCY_PREFIX=$COFFEE_DEP"
"$COFFEE_DEP/bin/cmake" --build \
  "$COFFEE_PROJECT/third_party/coffee-build-surface-anchored" \
  --target coffee-3d -j 4
