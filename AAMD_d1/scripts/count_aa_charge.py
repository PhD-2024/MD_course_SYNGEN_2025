#!/usr/bin/env python3
"""
Count amino-acid residues in a PDB and estimate net charge at pH 7.

Usage:
    python3 scripts/count_aa_charge.py /path/to/part_2.pdb

If no path is given the script will try to use
"files_afternoon_1/part_2.pdb" relative to the repository root.

Behavior:
 - Iterates over all lines and picks up lines starting with ATOM or HETATM.
 - Extracts chain ID, residue sequence number and residue name using PDB fixed-column parsing.
 - Builds a dict with counts of unique residues per residue name.
 - Uses a second dict mapping standard 3-letter amino-acid codes to their sidechain charge at pH 7
   (HIS treated as partially protonated: +0.1). Non-standard residues default to 0 charge.
 - Adds N-terminus (+1) and C-terminus (-1) charges for each chain.
 - Prints the counts, the charge-map used for found residues, and the estimated total charge.
"""
from collections import defaultdict
import sys
import os

# Standard amino acids 3-letter -> approximate sidechain charge at pH 7
# Note: HIS is partially protonated at physiological pH; here we use +0.1 to reflect that.
SIDECHAIN_CHARGE_P7 = {
    'ALA': 0.0,
    'ARG': +1.0,
    'ASN': 0.0,
    'ASP': -1.0,
    'CYS': 0.0,   # pKa ~8.3 -> mostly uncharged at pH 7
    'GLN': 0.0,
    'GLU': -1.0,
    'GLY': 0.0,
    # 'HIS': +0.1,  # partial protonation approximation
    'HIS': 0.0, 
    'ILE': 0.0,
    'LEU': 0.0,
    'LYS': +1.0,
    'MET': 0.0,
    'PHE': 0.0,
    'PRO': 0.0,
    'SER': 0.0,
    'THR': 0.0,
    'TRP': 0.0,
    'TYR': 0.0,
    'VAL': 0.0,
}


def parse_pdb_residues(pdb_path):
    """Parse PDB and return mapping of (chain, resid) -> resname (3-letter)."""
    residues = {}
    try:
        with open(pdb_path, 'r') as fh:
            for line in fh:
                if not (line.startswith('ATOM') or line.startswith('HETATM')):
                    continue
                # PDB fixed column format
                # resname is columns 18-20 (1-indexed) -> python indices [17:20]
                resname = line[17:20].strip()
                chain = line[21].strip() or '_'
                try:
                    resid = int(line[22:26].strip())
                except ValueError:
                    # fallback: use whole field as string if not integer
                    resid = line[22:26].strip()
                key = (chain, resid)
                # only set once (first occurrence of that residue)
                if key not in residues:
                    residues[key] = resname
    except FileNotFoundError:
        raise
    return residues


def count_by_resname(residues):
    """Given residues dict (chain,resid)->resname return counts per resname and chains' residue orders."""
    counts = defaultdict(int)
    chains = defaultdict(list)
    for (chain, resid), resname in residues.items():
        counts[resname] += 1
        chains[chain].append(resid)
    # sort residue numbers for each chain to detect termini
    for ch in chains:
        try:
            chains[ch] = sorted(chains[ch], key=lambda x: int(x))
        except Exception:
            chains[ch] = sorted(chains[ch])
    return counts, chains


def estimate_net_charge(counts, chains):
    """Estimate net charge from counts and add terminal contributions per chain.
    Returns (net_charge, used_charge_map, unknown_residues)
    """
    net = 0.0
    used_charge_map = {}
    unknown = set()
    for resname, count in counts.items():
        charge = SIDECHAIN_CHARGE_P7.get(resname)
        if charge is None:
            # unknown or non-standard residue: assume 0 but keep note
            unknown.add(resname)
            charge = 0.0
        used_charge_map[resname] = charge
        net += charge * count
    # Add terminal charges per chain: N-term +1, C-term -1
    # This assumes standard free termini (not capped).
    for ch, reslist in chains.items():
        if not reslist:
            continue
        net += 1.0  # N-terminus
        net -= 1.0  # C-terminus
        # net contribution is zero for a single continuous chain, but for completeness we keep it explicit
    return net, used_charge_map, unknown


if __name__ == '__main__':
    if len(sys.argv) > 1:
        pdb_file = sys.argv[1]
    else:
        # default relative path in repo
        pdb_file = os.path.join(os.path.dirname(__file__), '..', 'files_afternoon_1', 'part_2.pdb')
        pdb_file = os.path.normpath(pdb_file)

    try:
        residues = parse_pdb_residues(pdb_file)
    except FileNotFoundError:
        print(f"Error: pdb file not found: {pdb_file}")
        sys.exit(2)

    counts, chains = count_by_resname(residues)
    net, used_map, unknown = estimate_net_charge(counts, chains)

    # Output results
    print("Residue counts (3-letter codes):")
    for r, c in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
        print(f"  {r}: {c}")

    print("\nCharge map used (per residue, pH7):")
    for r, ch in sorted(used_map.items()):
        print(f"  {r}: {ch}")

    if unknown:
        print("\nUnknown / non-standard residues detected (assumed 0 charge):")
        for r in sorted(unknown):
            print(f"  {r}")

    print("\nChains and termini info:")
    for ch, reslist in sorted(chains.items()):
        start = reslist[0] if reslist else None
        end = reslist[-1] if reslist else None
        print(f"  Chain '{ch}': N-term resid={start}, C-term resid={end}, total residues={len(reslist)}")

    print(f"\nEstimated total net charge (approx.) at pH 7: {net:.2f}")

    # Also export small dict results if caller wants to import this script
    # Example usage from python:
    #   from scripts.count_aa_charge import parse_pdb_residues, count_by_resname, estimate_net_charge
    
