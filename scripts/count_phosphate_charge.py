#!/usr/bin/env python3
"""
Count phosphate groups in PDB files by counting P atoms and estimate negative charge.

Usage:
    python3 scripts/count_phosphate_charge.py file1.pdb file2.pdb

If no files are given it will default to:
  files_afternoon_1/part_1.pdb files_afternoon_1/part_0.pdb

Behavior:
 - Iterates over lines starting with ATOM or HETATM.
 - Uses fixed-column PDB parsing to extract atom name (cols 13-16 -> [12:16]).
 - Counts atoms with name 'P' and reports per-file counts and estimated charge (-1 per P).
"""
import sys
import os
from typing import List

DEFAULT_FILES = [
    os.path.join(os.path.dirname(__file__), '..', 'files_afternoon_1', 'part_1.pdb'),
    os.path.join(os.path.dirname(__file__), '..', 'files_afternoon_1', 'part_0.pdb'),
]


def count_p_atoms(pdb_path: str) -> int:
    count = 0
    with open(pdb_path, 'r') as fh:
        for line in fh:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            atomname = line[12:16].strip()
            # Standard phosphate atom name in nucleic acids is 'P'
            if atomname.upper() == 'P':
                count += 1
    return count


def main(files: List[str]):
    total_p = 0
    total_charge = 0.0
    for f in files:
        fpath = os.path.normpath(f)
        if not os.path.isfile(fpath):
            print(f"Warning: file not found: {fpath}")
            continue
        try:
            pcount = count_p_atoms(fpath)
        except Exception as e:
            print(f"Error reading {fpath}: {e}")
            continue
        charge = -1.0 * pcount  # each phosphate ~ -1 at pH 7
        print(f"{os.path.basename(fpath)}: P atoms = {pcount}, estimated negative charge = {charge:.1f}")
        total_p += pcount
        total_charge += charge

    print(f"\nTotal P atoms across files considered: {total_p}")
    print(f"Total estimated negative charge from phosphates: {total_charge:.1f}")


if __name__ == '__main__':
    if len(sys.argv) > 1:
        files = sys.argv[1:]
    else:
        files = DEFAULT_FILES
    main(files)
