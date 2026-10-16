#!/bin/sh
# Separate experimental source/binary: never replaces the upstream-boundary build.
set -eu
COFFEE_PROJECT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
COFFEE_DEP="$COFFEE_PROJECT/third_party/coffee-tools/env"
COFFEE_BASE="$COFFEE_PROJECT/third_party/coffee-surface-anchored"
COFFEE_VARIANT="$COFFEE_PROJECT/third_party/coffee-curl-conductor"
if [ ! -f "$COFFEE_BASE/CMakeLists.txt" ]; then
  "$COFFEE_PROJECT/run_coffee_build.sh"
fi
if [ ! -f "$COFFEE_VARIANT/CMakeLists.txt" ]; then
  python3 - "$COFFEE_BASE" "$COFFEE_VARIANT" <<'PYBUILD'
import os,shutil,sys
from pathlib import Path
base,dst=map(Path,sys.argv[1:])
dst.mkdir(parents=True)
for filename in ('CMakeLists.txt','LICENSE','README.md'):
    shutil.copyfile(base/filename,dst/filename)
# Dependencies and unchanged source remain read-only shared inputs.
(dst/'deps').symlink_to(base.resolve()/'deps',target_is_directory=True)
shutil.copytree(base/'src',dst/'src',copy_function=lambda a,b:os.symlink(Path(a).resolve(),b))
for filename in ('src/algorithms/field_solver_EZ.cpp','src/main-pulsar3d.cpp',
                 'src/main-verification.cpp','src/user_pulsar3d_EZ.hpp'):
    f=dst/filename; content=f.read_bytes(); f.unlink(); f.write_bytes(content)
PYBUILD
  (cd "$COFFEE_VARIANT" && patch -p1 < "$COFFEE_PROJECT/third_party/coffee-curl-conductor.patch")
fi
"$COFFEE_DEP/bin/cmake" -S "$COFFEE_VARIANT" \
  -B "$COFFEE_PROJECT/third_party/coffee-build-curl-conductor" \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_OSX_DEPLOYMENT_TARGET=14.0 \
  -DCMAKE_OSX_SYSROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk \
  -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ \
  -Duse_cuda=0 -Duse_double=1 -Dportable_scalar=ON -Dbuild_tests=OFF \
  "-DCOFFEE_DEPENDENCY_PREFIX=$COFFEE_DEP"
"$COFFEE_DEP/bin/cmake" --build \
  "$COFFEE_PROJECT/third_party/coffee-build-curl-conductor" \
  --target coffee-3d coffee-verification -j 4
