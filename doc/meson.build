###############################################################################################
################################################## constants ##################################
documentation_build_directory_name = 'doc_build'
copied_doc_source_directory_name = 'doc_src'
extra_data_path = copied_doc_source_directory_name/'_extra'
doxygen_output_folder_name = 'doxygen_xml'
doxyfile_path = meson.current_source_dir() / 'Doxyfile'
################################################ end constants ################################
###############################################################################################

###############################################################################################
############################## Update Doxygen configuration file ##############################
keyval = import('keyval')
doxygen_config = keyval.load(doxyfile_path)

doxyfile_lines = []
foreach k, v : doxygen_config
    if k == 'OUTPUT_DIRECTORY'
        v = '"@0@"'.format(documentation_build_directory_name)
    endif
    if k == 'XML_OUTPUT'
        v = '"@0@"'.format(doxygen_output_folder_name)
    endif
    doxyfile_lines += '@0@ = @1@\n'.format(k, v)
endforeach

new_config_content = '\n'.join(doxyfile_lines)

doxygen_conf = configuration_data()
doxygen_conf.set('CONTENT', new_config_content)

tpm_doxygen_conf = configure_file(
    output: 'Doxyfile.ini',
    configuration: doxygen_conf,
)
############################ end Update Doxygen configuration file ############################
###############################################################################################

doxygen = find_program('doxygen', required: false)
if not doxygen.found()
    warning('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    warning('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    warning('doxygen could not be found, thus the documentation can not be built.')
    warning('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    warning('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    subdir_done()
endif

doxygen_target = custom_target(
    'doxygen',
    input: tpm_doxygen_conf,
    output: doxygen_output_folder_name,
    command: [doxygen, '@INPUT@'],
)

sphinx_build = find_program('sphinx-build', required: false)
if not sphinx_build.found()
    warning('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    warning('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    warning(
        'doxygen could not be found. You may need to activate the virtual environment.',
    )
    warning('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    warning('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    subdir_done()
endif

documentation_build_directory = meson.project_build_root() / documentation_build_directory_name

env = environment()
env.append('DOXYGEN_BUILD_DIR', documentation_build_directory)
env.append('DOXYGEN_BUILD_SUB_DIR', doxygen_output_folder_name)

copy_doc_to_build = custom_target(
    command: [
        python3,
        meson.project_source_root() / 'scripts/copy_dir.py',
        meson.project_source_root() / 'doc',
        meson.project_build_root() / copied_doc_source_directory_name,
    ],
    output: [copied_doc_source_directory_name],
)

copy_clang_tidy_log_to_doc = custom_target(
    command: [
        python3,
        meson.project_source_root() / 'scripts/copy_file.py',
        meson.project_build_root() / clang_tidy_output_log_name,
        meson.project_build_root() / extra_data_path / clang_tidy_output_log_name,
    ],
    output: [clang_tidy_output_log_name],
    depends: [clang_tidy, copy_doc_to_build],
)

copy_coverage_report_to_doc = []

ninja = find_program('ninja')
if ninja.found()
    coverage_run = custom_target(
        'coverage-run',
        output: 'coverage-run.log',
        command: [ninja, '-C', meson.project_build_root(), 'coverage'],
        capture: true,
        build_always_stale: true,
        build_by_default: false,
    )
    coverage_report_path = 'meson-logs/coveragereport'
    coverage_destination = 'coveragereport'

    copy_coverage_report_to_doc += custom_target(
        'copy_coverage_report_to_doc',
        command: [
            python3,
            meson.project_source_root() / 'scripts/copy_dir.py',
            meson.project_build_root() / coverage_report_path,
            meson.project_build_root() / extra_data_path / coverage_destination,
        ],
        capture: true,
        output: ['copy_coverage_report_to_doc.log'],
        depends: [coverage_run, copy_doc_to_build],
    )

endif

shpinx_target_deps = [doxygen_target, copy_doc_to_build]

run_target(
    'sphinx-build',
    command: [
        sphinx_build,
        meson.project_build_root() / copied_doc_source_directory_name,
        documentation_build_directory,
    ],
    depends: shpinx_target_deps,
    env: env,
)


run_target(
    'sphinx-build-with-dashboard',
    command: [
        sphinx_build,
        meson.project_build_root() / copied_doc_source_directory_name,
        documentation_build_directory,
    ],
    depends: shpinx_target_deps + [copy_clang_tidy_log_to_doc] + copy_coverage_report_to_doc,
    env: env,
)

run_target(
    'get-sphinx-build-path',
    command: [
        python3,
        '-c',
        'import sys; print(80*"=");print("");print(sys.argv[1]);print("");',
        meson.project_build_root() / documentation_build_directory / 'index.html',
    ],
)
