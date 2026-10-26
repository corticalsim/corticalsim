#!/usr/bin/env python3

import sys
import os


def link(source, destination):
    if not os.path.exists(source):
        print(f"Source file '{source}' does not exist.")
        sys.exit(1)

    if os.path.exists(destination):
        print(f"Destination file '{destination}' already exists. overwriting!")
        os.remove(destination)

    os.symlink(source, destination)


if __name__ == "__main__":
    sys.argv.pop(0)
    if len(sys.argv) != 2:
        print("Usage: link_files.py <source> <destination>")
        sys.exit(1)
    source = sys.argv[0]
    destination = sys.argv[1]

    link(source, destination)
