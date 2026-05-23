#!/usr/bin/env python3
"""
Walks all subdirectories of the script's location.
- If a subdirectory contains statepoint.25.h5, deletes all statepoint.XX.h5 where XX < 25.
- If it does NOT contain statepoint.25.h5, prints the subdirectory name.
"""
 
import os
import re
import glob
 
script_dir = os.path.dirname(os.path.abspath(__file__))
 
for dirpath, dirnames, filenames in os.walk(script_dir):
    # Skip the root directory itself
    if dirpath == script_dir:
        continue
 
    has_25 = "statepoint.25.h5" in filenames
 
    if has_25:
        for fname in filenames:
            match = re.match(r"^statepoint\.(\d+)\.h5$", fname)
            if match and int(match.group(1)) < 25:
                filepath = os.path.join(dirpath, fname)
                # print(f"Deleting: {filepath}")
                os.remove(filepath)
    else:
        nums = []
        for fname in filenames:
            match = re.match(r"^statepoint\.(\d+)\.h5$", fname)
            if match:
                nums.append(int(match.group(1)))
        highest = max(nums) if nums else None
        if highest is not None:
            print(f"Missing statepoint.25.h5: {dirpath}  (highest: statepoint.{highest}.h5)")
        # else:
        #    print(f"Missing statepoint.25.h5: {dirpath}  (no statepoint files found)")
