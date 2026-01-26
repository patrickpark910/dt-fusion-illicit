#!/usr/bin/env python3
import os
import sys
import subprocess

def main():
    target = "./inputs/"
    if not os.path.isdir(target):
        print(f"Error: '{target}' is not a directory", file=sys.stderr)
        sys.exit(1)
    for entry in (os.listdir(target)): # reversed
        # Only process files ending with .inp
        if not entry.endswith('.inp'):
            continue
        path = os.path.join(target, entry)
        if not os.path.isfile(path):
            continue
        if os.path.exists(f"./outputs/o_{entry.rsplit('.', 1)[0]}.out"):
            continue
        cmd_to_run = f"mcnp6 i={path} o=./outputs/o_{entry.rsplit('.', 1)[0]}.out tasks 32"
        try:
            subprocess.run(cmd_to_run, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Command failed on '{path}' with exit code {e.returncode}", file=sys.stderr)

if __name__ == "__main__":
    main()