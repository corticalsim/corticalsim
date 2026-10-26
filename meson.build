project('corticalsim', 'cpp')
add_project_arguments(
  '-std=c++20',
  '-Wall',
  '-Wpedantic',
  '-Wextra',
  '-Wconversion',
  language: 'cpp',
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
# link compile_commands.json to the projec root
meson.add_postconf_script(
  find_program('python3'),
  meson.current_source_dir() / 'scripts/link_files.py',
  meson.current_build_dir() / 'compile_commands.json',
  meson.current_source_dir() / 'compile_commands.json',
)

