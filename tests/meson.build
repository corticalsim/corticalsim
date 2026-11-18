test_sources = files(
    'DLList_test.cpp',
)

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
