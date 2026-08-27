#!/bin/bash

# Variables
SCRIPT_DIR=`pwd`
ROOT_DIR=`dirname "$SCRIPT_DIR"`
BUILD_TYPE="release"
RECONFIGURE=""
CLEAN=false

# Move to the root directory (one up from this one)
cd $ROOT_DIR

# Get all arguments
ARG=( "${@}" )
for i in ${ARG[@]}; do
    if [ "$i" = "--reconfigure" ]; then
        RECONFIGURE='--reconfigure'
    fi
    if [ "$i" = "--clean" ]; then
        CLEAN=true
    fi
    if [ "$i" = "--buildtype" ]; then
        shift 1
        BUILD_TYPE="$1"
    fi
    shift 1
done

BUILD_DIR="build/$BUILD_TYPE"
echo ">>> Build directory: '$BUILD_DIR'"

if [ $CLEAN == true ]; then
    echo ">>> Cleaning..."
    rm -rf $BUILD_DIR
    echo ">>> Cleaning: done!"
fi

mkdir -p $BUILD_DIR

# Install dependency wraps
echo ">>> Installing dependencies..."
meson wrap install eigen 5.0.1-1
echo ">>> Installing dependencies: done!"

# Configure
echo ">>> Configuring the build..."
CONF_CMD="meson setup $BUILD_DIR . --buildtype $BUILD_TYPE $RECONFIGURE"
echo ">>> Configuration command: '$CONF_CMD'"
$CONF_CMD
echo ">>> Configuring the build: done!"

# Compile
echo ">>> Compiling CorticalSim..."
COMP_CMD="meson compile -C $BUILD_DIR"
echo ">>> Compilation command: '$COMP_CMD'"
$COMP_CMD
echo ">>> Compiling CorticalSim: done!"
