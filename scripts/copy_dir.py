#!/usr/bin/env python3

import sys
import os
import shutil


def main():
    sys.argv.pop(0)
    if len(sys.argv) != 2:
        print("Usage: copy_file.py <source> <destination>")
        sys.exit(1)

    source = sys.argv[0]
    destination = sys.argv[1]

    if os.path.isdir(destination):
        print("Destination directory exists. Overwriting!")
        shutil.rmtree(destination)
    shutil.copytree(source, destination)


if __name__ == "__main__":
    main()
