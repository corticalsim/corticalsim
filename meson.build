project(
    'corticalsim',
    'cpp',
    meson_version: '>=1.3.0',
    default_options: ['cpp_std=c++20'],
)

# check if a feature flag is defined
feature_flag_conditioned_dependencies = []
cpp_args = get_option('cpp_args')
if not cpp_args.contains('-DFEATURE_FLAG_BOOSTLESS')
    feature_flag_conditioned_dependencies = [
        # Boost pool and its sub-dependencies.
        dependency('boost-throw_exception'),
        dependency('boost-assert'),
        dependency('boost-integer'),
        dependency('boost-type_traits'),
        dependency('boost-config'),
        dependency('boost-pool'),
    ]
else
endif

# Sources
# =======
# Directories containing header files
incdir = include_directories('include')

# Directories containing source files
subdir('src')
subdir('tests')

# Dependencies
# ============
# Eigen
eigen = dependency('eigen3')

corticalsim_deps = [
    eigen,
] + feature_flag_conditioned_dependencies

# Executable(s)
# =============
corticalsim_lib = library(
    'libcorticalsim',
    lib_sources,
    include_directories: incdir,
    dependencies: corticalsim_deps,
)

corticalsim = executable(
    'corticalsim',
    exec_sources,
    include_directories: incdir,
    link_with: corticalsim_lib,
    dependencies: [
        corticalsim_deps,
    ],
)

sources = lib_sources + exec_sources

############################################################################
############################## Tooling #####################################
############################################################################

# Google Test

gtest = dependency('gtest_main')

if gtest.found()
    test_executable = executable(
        'all_tests',
        test_sources,
        include_directories: incdir,
        link_with: corticalsim_lib,
        dependencies: [
            gtest,
            corticalsim_deps,
        ],
    )

    test('all_tests', test_executable)
endif

# link compile_commands.json to the projec root
meson.add_postconf_script(
    find_program('python3'),
    meson.current_source_dir() / 'scripts/link_files.py',
    meson.current_build_dir() / 'compile_commands.json',
    meson.current_source_dir() / 'compile_commands.json',
)

# Clang-format

# copy .clang-format file to the build directory
configure_file(input: '.clang-format', output: '.clang-format', copy: true)

run_target(
    'format',
    command: [
        'clang-format',
        '-i',
        '-style=file',
        sources,
    ],
)

absolute_paths_to_sources = []

clang_tidy_config_file = get_option('config-file')
warnings_as_errors = get_option('warnings-as-errors')

run_target(
    'tidy',
    command: [
        find_program('python3'),
        meson.current_source_dir() / 'scripts/run_tidy.py',
        meson.current_build_dir(),
        meson.current_source_dir(),
        clang_tidy_config_file,
        warnings_as_errors,
        sources,
    ],
    depends: [corticalsim],
)
