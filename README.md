# corticalsim3D (cortical microtubule dynamics)

![alt text](assets/logo/logomodelCMT-1.0.png)

This is a program for simulating cortical microtubule dynamics (CMT) on experimentally extracted microscopic images of cells

Application:

https://www.sciencedirect.com/science/article/pii/S0960982218309242

https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005959

## Note

Testing and contributing is very welcome, especially if you can contribute with new algorithms and features.

# Installation

NOTE: You need a compiler that supports the C++ 14 standard as a minimum.

## Install the dependencies

`corticalsim3D` depends on Eigen and Boost (however, the latter will be dropped soon). In addition, you will need [Meson](https://mesonbuild.com/index.html) to build the executables, as well as a fairly recent C++ compiler (GCC or Clang) supporting the C++ 14 standard. The following instructions outline the process of installing the dependencies on different platforms, as well as compiling and installing the `corticalsim3D` package itself.

### Windows

Most dependencies need to be installed manually on Windows.

1. GCC

    1. Download [MinGW](https://www.mingw-w64.org/)
    1. Run `mingw-xx.exe` (this will create a folder called `MinGw`)
    1. Copy the `MinGw` folder to `C:\`
    1. Add `C:\MinGW` and `C:\MinGW\bin` to the system path

1. Eigen

    1. Download [Eigen](https://libeigen.gitlab.io/)
    1. Extract `eigen-3.4.0` folder and rename it to `Eigen`
    1. Copy the `Eigen` folder to `C:\Program Files`

1. Boost

    Download and install [Boost](https://www.boost.org/releases/latest/).

    <span style="color: red;">NOTE</span>: The Boost dependency will be deprecated and removed soon.

1. Meson

    Use the official [Meson installer for Windows](https://github.com/mesonbuild/meson/releases). The installer also bundles the Ninja dependency, so you don't have to install that separately.

### Linux

Please use your package manager to install all the dependencies. Most Linux distributions provide packages for GCC, Meson, Ninja, Boost and Eigen.

Instructions for Ubuntu (*to be confirmed, may be incomplete*):

```bash
sudo apt update && apt upgrade
sudo apt install gcc libboost-filesystem-dev libeigen3-dev meson ninja-build cmake
```

## Compile

1. First, clone the `corticalsim3D` repository:

    ```bash
    git clone https://github.com/NLeSC/corticalsim3D
    ```

    NOTE: It is recommended that you do **_not_** enter the `corticalsim3D` directory. The next steps assume an out-of-source build setup, meaning that the directory (where the package will be compiled) will be outside the source tree. This has several advantages, such as precluding the possibility of accidentally committing build artefacts to the repository.

1. Configure the build with `meson`:

    ```bash
    meson setup build corticalsim3D
    ```

    This command uses the Meson build system to bootstrap everything that we need for the compilation into a directory called `build` from the sources located in the `corticalsim3D` directory.

1. Enter the build directory and compile the source:

    ```bash
    meson compile -C build
    ```

    This will generate an executable file called `corticalsim3d` in the current directory (`build`).

## Run

Run the `corticalsim3d` executable with the `corticalsim3D/config/parameters_ARRAY.txt` parameter file as the first argument:

```bash
./build/corticalsim3d ./corticalsim3D/config/parameters_ARRAY.txt
```
