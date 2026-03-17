project('corticalsim', 'cpp')

# Compiler flags
# add_global_arguments('-O2', language : 'cpp')

# Directories containing header files
incdir = include_directories('include')

# Directories containing source files
subdir('src')

# Dependencies
eigen = dependency('eigen3')
boost = dependency('boost', modules : ['filesystem'])

# Executable(s)
executable('corticalsim3d', sources, include_directories: incdir, dependencies: [eigen, boost])
