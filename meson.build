project(
    'corticalsim',
    'cpp',
    meson_version: '>=1.3.0',
    default_options: [
        'cpp_std=c++20',
        'b_coverage=true', # enable coverage targets if gcovr is available
    ],
)

# check if a feature flag is defined
feature_flag_conditioned_dependencies = []
if not get_option('boostless').enabled()
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
    message('Compiling without boost')
    add_global_arguments('-DFEATURE_FLAG_BOOSTLESS', language: 'cpp')
endif

# Sources
# =======
# Directories containing header files
incdir = include_directories('include')

# Directories containing source files
subdir('src')

# Dependencies
# ============
# Eigen
eigen = dependency('eigen3')

corticalsim_deps = [
    eigen,
] + feature_flag_conditioned_dependencies

corticalsim_lib = library(
    'libcorticalsim',
    lib_sources,
    include_directories: incdir,
    dependencies: corticalsim_deps,
)

# Executable(s)
# =============
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

# to be platform independent, we implement all our dev tooling scripts in python. This however means that python is a
# dev dependency
python3 = find_program('python3', required: false)
if not python3.found()
    warning('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    warning('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    warning(
        'no python3 found. This means that some dev dependencies will not work. Main target can still be built. Please install python3. See README.dev.md for further information.',
    )
    warning('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    warning('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    subdir_done()
endif
# link compile_commands.json to the projec root
meson.add_postconf_script(
    python3,
    meson.current_source_dir() / 'scripts/link_file.py',
    meson.current_build_dir() / 'compile_commands.json',
    meson.current_source_dir() / 'compile_commands.json',
)

# Google Test

subdir('tests')


# Clang-format
# clang format target is automatically added by meson if clang-format is installed

clang_tidy_config_file = get_option('config-file')
warnings_as_errors = get_option('warnings-as-errors')

clang_tidy_bin = find_program('clang-tidy', required: false)
run_clang_tidy = find_program('run-clang-tidy', required: false)

if run_clang_tidy.found() and clang_tidy_bin.found()

    tidy_script = meson.current_source_dir() / 'scripts/run_tidy.py'

    clang_tidy_output_log_name = 'clang-tidy.log'
    clang_tidy = custom_target(
        'tidy',
        command: [
            python3,
            '@INPUT0@',
            '@INPUT1@',
            '@INPUT2@',
            meson.current_build_dir(),
            meson.current_source_dir(),
            '@OUTPUT@',
            '@INPUT3@',
            warnings_as_errors,
            sources,
        ],
        input: [
            tidy_script,
            clang_tidy_bin.full_path(),
            run_clang_tidy.full_path(),
            meson.project_source_root() / clang_tidy_config_file,
        ],
        output: clang_tidy_output_log_name,
    )

endif

subdir('doc')
