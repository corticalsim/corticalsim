#!/bin/bash

# Functions
# =========

function run_or_bail() {
    local MESSAGE="$1"
    local COMMAND="$2"
    
    echo ">>> ${MESSAGE}..."
    $COMMAND
    if [ $? != 0 ]; then
        echo "!!! ${MESSAGE}: failed."
        exit $?
    fi
    echo ">>> ${MESSAGE}: done!"
} 

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
ARGS=( "${@}" )
for i in ${ARGS[@]}; do
    case "$i" in

        "-r" | "--reconfigure")
            SETUP_ARGS="$SETUP_ARGS --reconfigure --wipe --clearcache"
            ;;

        "-d" | "--builddir")
            shift
            BUILD_DIR="$1"
            ;;

        "-t" | "--buildtype")
            shift
            BUILD_TYPE="$1"
            ;;

        "-c" | "--clean")
            CLEAN_BUILD=true
            ;;
            
        "-j" | "--jobs")
            shift
            BUILD_JOBS="$1"
            ;;
    esac
    shift
done

# Set the build type
echo ">>> Build type: '${BUILD_TYPE}'"
SETUP_ARGS="${SETUP_ARGS} --buildtype=${BUILD_TYPE}"

# Create the build directory
BUILD_DIR="$(readlink -m ${BUILD_DIR})/${BUILD_TYPE}"
mkdir -p $BUILD_DIR
echo ">>> Build directory: '${BUILD_DIR}'"

if [ -n "${BUILD_JOBS}" ]; then
    COMPILE_ARGS="${COMPILE_ARGS} -j ${BUILD_JOBS}"
fi

# Setup
# =====
SETUP_CMD=$(echo "meson setup ${BUILD_DIR} ${ROOT_DIR} ${SETUP_ARGS}" | xargs)
run_or_bail "Running setup" "${SETUP_CMD}" 

# Compilation
# ===========

# Clean the build directory if requested
if [ -n "${CLEAN_BUILD}" ]; then
    CLEAN_CMD="meson compile -C ${BUILD_DIR} --clean"
    run_or_bail "Cleaning the build directory" "${CLEAN_CMD}"
fi

COMPILATION_CMD=$(echo "meson compile -C ${BUILD_DIR} ${COMPILE_ARGS}" | xargs)
run_or_bail "Compiling CorticalSim" "${COMPILATION_CMD}"
