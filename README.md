# corticalsim3D (cortical microtubule dynamics)

![alt text](doc/assets/logo/logomodelCMT-1.0.png)

CorticalSim is a simulator for cortical microtubule dynamics (CMT) on experimentally extracted microscopic images of cells.

Application:

https://www.sciencedirect.com/science/article/pii/S0960982218309242

https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005959

## Note

Testing and contributing is very welcome, especially if you can contribute with new algorithms and features.

# Installation

NOTE: You need a compiler that supports the C++ 14 standard as a minimum.

## Install the dependencies

Currently, `corticalsim3D` depends on Eigen and Boost (the latter is going to be removed as a dependency at some point). These dependencies are managed automatically by [Meson](https://mesonbuild.com/index.html), which is also used to build the executables. You do need to install a fairly recent C++ compiler (GCC or Clang) supporting the C++ 14 standard. The following instructions outline the process of installing the core dependencies on different platforms, as well as compiling and installing `corticalsim3D` itself.

### Linux

The dependencies are handled by Meson, so the only packages that you need to install manually are GCC, Meson, Ninja and CMake.

#### Ubuntu

```bash
sudo apt update && apt upgrade && apt install gcc meson ninja-build cmake
```

#### Arch

```bash
sudo pacman -Sy gcc meson ninja cmake
```

### Windows

1. Meson

    Use the official [Meson MSI installer for Windows](https://github.com/mesonbuild/meson/releases). The installer also bundles the Ninja dependency, so you don't have to install that separately.

1. GCC

    1. Download [MinGW](https://www.mingw-w64.org/)
    1. Run `mingw-xx.exe` (this will create a folder called `MinGw`)
    1. Copy the `MinGw` folder to `C:\`
    1. Add `C:\MinGW` and `C:\MinGW\bin` to the system path

## Compilation

First, clone the `corticalsim3D` repository:

```bash
git clone https://github.com/corticalsim/corticalsim3D
cd corticalsim3D
```

### Build script

You can use the build script located under the `scripts` directory:

```bash
cd scripts
./build.sh
```

### Manual build

Alternatively, you can run all the setup and compilation steps manually by following the steps below:

1. From the root directory, first install the dependency wraps:

```bash
    meson wrap install eigen 5.0.1-1
```

1. Run the setup to configure the build:

```bash
mkdir -p build/release
meson setup build/release . --buildtype release"
```

This command uses the Meson build system to bootstrap everything that we need for the compilation into a directory called `build/release` from the sources located in the current directory.

NOTE: The possible options for the `--buildtype` switch are `plain`, `debug`, `debugoptimized` and `release`. Please check the [Meson documentation page](https://mesonbuild.com/Running-Meson.html#configuring-the-build-directory) on configuring the build directory for information about the meaning of these options. We default to `release` here.

Different build types should go in their own directory, so if you would like to build corticalSim with debug symbols, first you should create a subdirectory for that build:

```bash
mkdir -p build/debug
meson setup build/debug . --buildtype debug
```

1. Compile the source by passing the relevant directory via the `-C` switch.

```bash
meson compile -C build/release
```

This will generate an executable file called `corticalsim3d` in the build directory (e.g., `build/release`).

NOTE: If you are using a very old Meson version (e.g., `<0.54.0`), you need to use `ninja` instead of `meson`, with a slightly different syntax:

```bash
ninja -C build/release
```

## Run

Run the `corticalsim3d` executable from the build directory with the `./config/parameters_ARRAY.txt` parameter file as the first argument:

```bash
./build/release/corticalsim3d ./config/parameters_ARRAY.txt
```

You should see the output of the run, with some statistics at the end. For instance, with a `debug` build:

```bash
...

Division plane area: 72.097664
corticalSim vBC.2018
Total running time: 8 seconds.
stochastic events: 1784, deterministic events: 56882.
```

Compare this with the `release` build:

```bash
...

Division plane area: 72.097664
corticalSim vBC.2018
Total running time: 3 seconds.
stochastic events: 1784, deterministic events: 56882.
```
