# corticalsim development documentation


> [!NOTE]
> Testing and contributing is very welcome, especially if you can contribute new algorithms and features.

> [!WARNING]
> For the development setup we strongly recommend a Linux system. While we try to use cross-platform tools, some things
> might not work as expected.


# Development

## Building the documentation locally

Documentation is generated with [Sphinx](https://www.sphinx-doc.org/) and [Doxygen](https://www.doxygen.nl/) (for automatically generating the technical documentation), with [Breathe](https://breathe.readthedocs.io/) providing a bridge between the Sphinx and Doxygen documentation systems. To build the documentation locally, follow these steps:

You will need to install both Doxygen and Sphinx to build the documentation locally.

### Doxygen

For most Linux distributions, doxygen is available as part of the official package repository

Ubuntu:
```bash
sudo apt install doxygen
```

Fedora:
```bash
sudo dnf install doxygen
```

macOS:
```bash
brew install doxygen
```
Or download the disk image from [the official website](https://www.doxygen.nl/download.html).

Windows:
To obtain the installation binaries for windows follow [the official documentation](https://www.doxygen.nl/download.html).

Building from source:
To build from source follow [the official documentation](https://www.doxygen.nl/manual/install.html).

### Sphinx

Then, install Sphinx and the other Python dependencies from the root of the repository:

```
python -m venv .venv
source .venv/bin/activate
pip install -r doc/requirements.txt
```

Then generate the Doxygen technical documentation XML files and use Sphinx to generate the HTML files with the full documentation:

```
doxygen doc/Doxyfile
sphinx-build doc doc/build
```

And we can view the generated documentation by opening `doc/build/index.html` in a browser.

## Code Style

Our code is automatically formatted using `clang-format`. We plan to enforce the formatting in the CI. This means that
merging will be blocked if the code is not fully formatted. Either set up your editor to automatically format code
while editing or when saving or run clang format manually before submitting a PR. The formatting guidelines are definded in the `.clang-format` file.

> [!WARNING]
> This method only works if you can use ninja as the meson backend. This works by default on Linux. And was not tested
> on other platforms.

For the manual formatting run:

```shell
ninja -C build/release clang-format
```

## Static code analysis

To improve code quality, we are introducing `clan-tidy` for static code analysis.
No `clang-tidy` checks are enforced at the moment, but the project root already contains a `.clang-tidy` file, which
will probably generate new warnings in your IDE or editor. It is good practice to not generate new clang tidy
warnings during development, but to fix them before committing.

### Installing `clang-tidy`

`clang-tidy` is part of [clang-tools-extra](https://clang.llvm.org/extra/). So depending on the operating system it
will be installed either as part of the llvm or clang-tools-extra package.

On fedora use the following command:
```shell
sudo dnf install clang-tools-extra
sudo dnf install llvm
```

### Running `clang-tidy`

Run `clang-tidy` on the project you first need to build the project, using the meson build system.

So the standard workflow is to first build the system

```bash
meson setup build/release . --buildtype release --reconfigure
```
and then to run the `tidy` target
```bash
ninja -C build/release tidy
```
This will run clang-tidy on all source files located in the `src/` folder and output all the warnings to the
terminal as well as to the `clang-tidy.log` file, which will be created in the project root.
In the CI workflow, we add a configuration step to promote analyzer warnings to errors, by setting

```bash
meson configure build/release  -Dwarnings-as-errors='clang-analyzer-*'
```

This step is happening before the `tidy` target is run.
