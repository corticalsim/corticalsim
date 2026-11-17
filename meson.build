project(
    'corticalsim',
    'cpp',
    meson_version: '>=1.3.0',
    default_options: ['cpp_std=c++20'],
)

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

# Boost pool and its sub-dependencies.
boost_throw_exception = dependency('boost-throw_exception')
boost_assert = dependency('boost-assert')
boost_integer = dependency('boost-integer')
boost_type_traits = dependency('boost-type_traits')
boost_config = dependency('boost-config')
boost_pool = dependency('boost-pool')

# Executable(s)
# =============
corticalsim = executable(
    'corticalsim',
    sources,
    include_directories: incdir,
    dependencies: [
        eigen,
        boost_assert,
        boost_integer,
        boost_throw_exception,
        boost_type_traits,
        boost_config,
        boost_pool,
    ],
)

# Tooling

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
