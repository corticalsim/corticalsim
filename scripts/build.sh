#!/bin/bash

# Variables
# =========

# This will extract an invariant path regardless of
# where the script has been called from.
SCRIPT_DIR="$(cd -- "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd)"

# The repository root.
ROOT_DIR=$(dirname "$SCRIPT_DIR")

# Build type.
# See https://mesonbuild.com/Running-Meson.html#configuring-the-build-directory
# Possible values: "plain", "debug", "debugoptimized", "release"
BUILD_TYPE="release"

# Aguments passed to `meson setup`.
SETUP_ARGS=""

# Aguments passed to `meson compile`.
COMPILE_ARGS=""

# Default build directory.
BUILD_DIR="$ROOT_DIR/build"

# CLI arguments
# =============
ARG=( "${@}" )
for i in ${ARG[@]}; do
    case "$i" in

        "--reconfigure")
            SETUP_ARGS="$SETUP_ARGS --reconfigure"
            ;;

        "--builddir")
            shift 1
            BUILD_DIR="$1"
            ;;

        "--buildtype")
            shift 1
            BUILD_TYPE="$1"
            SETUP_ARGS="$SETUP_ARGS --buildtype=$BUILD_TYPE"
            ;;

        "--clean")
            COMPILE_ARGS="$COMPILE_ARGS --clean"
            ;;
    esac
    shift 1
done

# Create the build directory
BUILD_DIR="$(readlink -m $BUILD_DIR)/$BUILD_TYPE"
mkdir -p $BUILD_DIR
echo ">>> Build directory: '$BUILD_DIR'"

# Setup
# =====
echo ">>> Running setup..."
SETUP_CMD=$(echo "meson setup $BUILD_DIR $ROOT_DIR $SETUP_ARGS" | xargs)
echo ">>> Setup command: '$SETUP_CMD'"
# $SETUP_CMD
echo ">>> Running setup: done!"

# Compilation
# ===========
echo ">>> Compiling CorticalSim..."
COMPILATION_CMD=$(echo "meson compile -C $BUILD_DIR $COMPILE_ARGS" | xargs)
echo ">>> Compilation command: '$COMPILATION_CMD'"
# $COMPILATION_CMD
echo ">>> Compiling CorticalSim: done!"
