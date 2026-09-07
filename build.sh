#!/usr/bin/env bash
# Configures and builds numsolve with CMake + MinGW g++, using the vcpkg
# instance at C:/vcpkg for nlohmann-json and Catch2 (same toolchain set up
# for the sibling firewall-rule-engine-cpp project -- see its README for
# how vcpkg was installed).
set -euo pipefail

cd "$(dirname "$0")"

VCPKG_ROOT="${VCPKG_ROOT:-C:/vcpkg}"

cmake -S . -B build \
    -G "MinGW Makefiles" \
    -DCMAKE_TOOLCHAIN_FILE="$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" \
    -DVCPKG_TARGET_TRIPLET=x64-mingw-static \
    -DVCPKG_HOST_TRIPLET=x64-mingw-static \
    -DCMAKE_BUILD_TYPE=Release

cmake --build build --parallel

echo ""
echo "Built: build/numsolve.exe"
echo "Test binary: build/numsolve_tests.exe"
