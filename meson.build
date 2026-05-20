project('corticalsim', 'cpp')

cc = meson.get_compiler('cpp')
if cc.get_id() == 'gcc'
    assert(cc.version() <= '16', 'Please use GCC<16. Support for GCC>=16 is work in progress.')
endif

# Directories containing header files
incdir = include_directories('include')

# Directories containing source files
subdir('src')

# Dependencies
eigen = dependency('Eigen3', modules: ['Eigen3::Eigen'])
boost = dependency('boost', modules : ['filesystem'])

# Executable(s)
executable('corticalsim3d', sources, include_directories: incdir, dependencies: [eigen, boost])
