# corticalsim development documentation

> [!NOTE]
> Testing and contributing is very welcome, especially if you can contribute new algorithms and features.

> [!WARNING]
> For the development setup we strongly recommend a Linux system. While we try to use cross-platform tools, some things
> might not work as expected.


# Development

## Core dependencies

In order to be able to build the main dependencies of the project you will require:
 - c++ compiler (preferably GCC or clang)
 - python (this is a dependency of meson)
 - meson (build system)
 - ninja (this is the build system backend, other backends work with meson as well, but we do not test them)
 - CMake (`corticalsim` still depends on third party libraries like `boost` and `Eigen`. Meson may need `CMake` to search if they are already installed on the system. Dependencies that can not be found will then be installed by meson).

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

In the active virtual environment (i.e. once `source .venv/bin/activate` has been run) you will have access to
additional meson targets. Which build the documentation.
In order to simply build the documentation run

```
meson setup build/debug . --buildtype debug
ninja -C build/debug sphinx-build
```

This will generate a `doc_build` directorey in the build directory which will contain the generated static website of
the documentation. Launching the `index.html` will display the documentation.
To obtain the exact path to this `index.html` file, run
```
ninja -C build/debug get-sphinx-build-path
```
The documentation contains a dashboard, which tracks  the quality metrics of the project.
The dashboard is automatically generated and displayed in the online version of the documentation.
If you want to generate the dashboard locally as well then you have to run
```
ninja -C build/debug sphinx-build-with-dashboard
```
instead of the
```
ninja -C build/debug sphinx-build
```
However, be aware that this will rerun all the metrics generators, like `clang-tidy` and `coverage`, which can take significantly longer.
```

## Code Style

Our code is automatically formatted using `clang-format`. We plan to enforce the formatting in the CI. This means that
merging will be blocked if the code is not fully formatted. Either set up your editor to automatically format code
while editing or when saving or run clang format manually before submitting a PR. The formatting guidelines are defined in the `.clang-format` file.

First install clang-format which is easiest to do with pip (on any platform):

```
pip install clang-format
```


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

On any Linux or MacOS system you can also use [homebrew](https://brew.sh/) to install llvm:
```shell
brew install llvm
```
This will also instll `clang-tidy` automatically.

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


## Testing

## Unit-Tests

We use [gtest](https://google.github.io/googletest/) for unit testing. The dependency is automatically installed by
the meson build system.

To run all test targets you can run:
```bash
meson test -C <path/to/build> --print-errorlogs
```
after setting up the build folder.

To run only a specific test target, use the following command:
```bash
meson test -C <path/to/build> --print-errorlogs <target_name>
```

### How to write unit tests

We use the Given/When/Then pattern for unit tests. This means that each unit test case has three sections (clauses).
Each clause starts with a comment containing the clause name and description, followed by the code that accomplishes
the described behavior.

Example:

```cpp
TEST(DLListTest, DLListSizeReportsCorrectSize)
{
    // Given a newly created DLList
    DLList<DummyTestType> dl_list{};

    // When we create several new elements
    constexpr std::size_t num_new_elements{ 17 };

    for (std::size_t i = 0; i < num_new_elements; ++i)
    {
        std::ignore = dl_list.create();
    }

    // Then the list size corresponds to the number of new elements
    EXPECT_EQ(dl_list.size(), num_new_elements);
}
```

**Given clause**

Establish the initial known state of the system. To keep the clause simple only address the initial state relevant to the test.

**When clause**

Describe the action, which you want to test.

**Then clause**

Describe the expected outcome of the tested action.


**Naming the test**

In gTest each test case starts with a `TEST` macro which takes two arguments, `test_suite_name` and `test_name`, both of which have to be valid c++ identifiers. `test_suite_name` should be shared between tests that test the same interface, and each `test_name` must be unique.
Since these names will later appear in test logs it is useful to make them descriptive.
The suggested convention for the `test_suite_name` is to use the related class or function name followed by the word "test".
The `test_name` should contain a full sentence which states the tested behavior and expected outcome.

Example:

Say we have Parameters class which contains all simulation parameters, and this class has a `load_from_file` method
with the following signature:
```cpp
static std::optional<Parameters> Parameters::load_from_file(std::string_view path);
```
We will want to test what happens if we provide a file path that does not exist.
If we can not load the parameters from the file then we want to return an empty optional.
A possible test case would be named as follows:

```cpp
TEST(ParametersTest, load_from_file_returns_empty_optional_when_file_does_not_exist)
```

## Test Coverage

Requires [`gcovr`](https://gcovr.com/en/stable/installation.html) to be installed. If `gcovr` is not installed then the
coverage target will not be created.

To get generate the local coverage metrics run the test suite first:
```bash
ninja test -C build/debug/
```
Then run the `coverage` target:
```bash
ninja -C build/debug/ coverage
```
This will generate the coverage and print the location of the logs to the consol. Open the index.html at the printed
location. Which will be probably at `corticalsim/build/debug/meson-logs/coveragereport/index.html`.

