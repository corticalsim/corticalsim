#!/usr/bin/env python3

import sys
import shutil


def main():
    sys.argv.pop(0)
    if len(sys.argv) != 2:
        print("Usage: copy_file.py <source> <destination>")
        sys.exit(1)

    source = sys.argv[0]
    destination = sys.argv[1]

    shutil.copy(source, destination)


if __name__ == "__main__":
    main()
