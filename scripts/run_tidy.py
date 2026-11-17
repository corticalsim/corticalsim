#!/usr/bin/env python3
"""
This is an internal script and not part of the public API.
It runs clang-tidy for the build sytem.
It is more usefull to run the clang-tidy through this script instead of
directly through meson since we can first validate that the configuration file
is in avalid state, and then post process the output.
"""

import subprocess
import sys
import os
from typing import List


def get_core_count(files: List[str]):
    max_available_cores = max(1, os.cpu_count()-1)
    use_cores = min(max_available_cores, len(files))
    return use_cores


def check_config_integrity(clang_tidy_path: str, clang_tidy_config_file: str):
    if not os.path.exists(clang_tidy_config_file):
        print(f"Config file '{clang_tidy_config_file}' does not exist.")
        sys.exit(1)
    proc = subprocess.run(
        [clang_tidy_path,  '--verify-config', clang_tidy_config_file],
        capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"Config file '{clang_tidy_config_file}' is not valid.")
        print(proc.stderr)
        print(proc.stdout)
        sys.exit(1)


def fix_paths(text: str, build_dir: str, source_dir: str):
    lines = text.splitlines()
    out = []
    for line in lines:
        # clang-tidy format: path:line:col: message
        parts = line.split(':', 3)
        if len(parts) >= 3 and parts[0].startswith('../'):
            rel_path = parts[0]
            abs_path = os.path.abspath(os.path.join(build_dir, rel_path))
            source_rel_path = os.path.relpath(abs_path, source_dir)

            parts[0] = source_rel_path
            line = ':'.join(parts)

        out.append(line)
    return '\n'.join(out)


def main():
    sys.argv.pop(0)

    if len(sys.argv) < 6:
        print("Usage: python run_tidy.py",
              "<clang-tidy path>",
              "<run-clang-tidy path>",
              "<build dir>",
              "<source dir>",
              "<output fiel name>",
              "<path/to/clang/tidy/config>",
              "<warning-as-errors-setting>",
              "<files...>")
        sys.exit(1)

    clang_tidy_path = sys.argv[0]
    run_clang_tidy_path = sys.argv[1]
    build_dir = sys.argv[2]
    source_dir = sys.argv[3]
    output_fiel_name = sys.argv[4]
    clang_tidy_config_file = sys.argv[5]
    warnings_as_errors = sys.argv[6]

    files = sys.argv[7:]

    cmd = [run_clang_tidy_path, '-p', build_dir, '-j',
           str(get_core_count(files)),
           '-config-file', clang_tidy_config_file,
           ]

    if warnings_as_errors:
        cmd.extend(['-warnings-as-errors', warnings_as_errors])

    abs_path_files = []
    for file in files:
        abs_path_files.append(os.path.abspath(os.path.join(build_dir, file)))

    cmd.extend(abs_path_files)

    check_config_integrity(clang_tidy_path, clang_tidy_config_file)
    print(f"Provided clang-tidy config at {clang_tidy_config_file} is valid!")

    proc = subprocess.run(cmd, capture_output=True, text=True)
    stdout = fix_paths(proc.stdout, build_dir, source_dir)
    stderr = fix_paths(proc.stderr, build_dir, source_dir)

    sys.stdout.write(stdout)
    sys.stderr.write(stderr)

    with open(output_fiel_name, 'w') as f:
        f.write(stdout)
        f.write(stderr)

    sys.exit(proc.returncode)


if __name__ == "__main__":
    main()
