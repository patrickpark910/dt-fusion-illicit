#!/usr/bin/env python3
"""Convert X column from MeV to eV (multiply by 1e6). Y column unchanged.
Usage: python convert_MeV_to_eV.py input.txt output.txt
"""
import sys

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} input.txt output.txt")
    sys.exit(1)

with open(sys.argv[1]) as f:
    lines = f.readlines()

out = []
for line in lines:
    stripped = line.strip()
    if stripped.startswith("#"):
        out.append(line.replace("MeV", "eV"))
        continue
    if not stripped:
        out.append(line)
        continue
    parts = stripped.split()
    if len(parts) >= 2:
        try:
            x_ev = float(parts[0]) * 1e6
            out.append(f"  {x_ev:>14.6g}      {parts[1]} \n")
        except ValueError:
            out.append(line)

with open(sys.argv[2], "w") as f:
    f.writelines(out)

print(f"Converted {len(out)} lines. X column multiplied by 1e6 (MeV -> eV).")