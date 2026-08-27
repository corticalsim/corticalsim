project('corticalsim', 'cpp')

# Sources
# =======
# Directories containing header files
incdir = include_directories('include')

# Directories containing source files
subdir('src')

# Dependencies
# ============
# Eigen is actually version 5, despite the name of the exposed dependency (eigen3).
eigen = dependency('eigen3')

# Custom Boost dependency wraps.
boost_throw_exception = dependency('boost-throw_exception')
boost_assert = dependency('boost-assert')
boost_integer = dependency('boost-integer')
boost_type_traits = dependency('boost-type_traits')
boost_config = dependency('boost-config')
boost_pool = dependency('boost-pool')

# Executable(s)
# =============
executable(
  'corticalsim3d',
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
