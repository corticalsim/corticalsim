#!/usr/bin/env python3

import subprocess
import sys
import os


def get_core_count():
    max_available_cores = max(1, os.cpu_count()-1)
    use_cores = min(max_available_cores, len(files))
    return use_cores


def check_config_integrity(clang_tidy_config_file):
    if not os.path.exists(clang_tidy_config_file):
        print(f"Config file '{clang_tidy_config_file}' does not exist.")
        sys.exit(1)

    proc = subprocess.run(
        ['clang-tidy',  '--verify-config', clang_tidy_config_file],
        capture_output=True, text=True)

    if proc.returncode != 0:
        print(f"Config file '{clang_tidy_config_file}' is not valid.")
        print(proc.stderr)
        print(proc.stdout)
        sys.exit(1)


if __name__ == "__main__":
    sys.argv.pop(0)
    if len(sys.argv) < 3:
        print("Usage: python run_tidy.py <build dir> <source dir> <files...>")
        sys.exit(1)

    build_dir = sys.argv[0]
    source_dir = sys.argv[1]
    clang_tidy_config_file = sys.argv[2]
    warnings_as_errors = sys.argv[3]

    files = sys.argv[4:]

    cmd = ['run-clang-tidy', '-p', build_dir, '-j', str(get_core_count()),
           '-config-file', clang_tidy_config_file,
           ]

    if warnings_as_errors:
        cmd.extend(['-warnings-as-errors', warnings_as_errors])

    abs_path_files = []
    for file in files:
        abs_path_files.append(os.path.abspath(os.path.join(build_dir, file)))

    cmd.extend(abs_path_files)

    os.chdir(source_dir)
    proc = subprocess.run(cmd, capture_output=True, text=True)

    def fix_paths(text):
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

    stdout = fix_paths(proc.stdout)
    stderr = fix_paths(proc.stderr)

    sys.stdout.write(stdout)
    sys.stderr.write(stderr)

    with open('clang-tidy.log', 'w') as f:
        f.write(stdout)
        f.write(stderr)

    sys.exit(proc.returncode)
