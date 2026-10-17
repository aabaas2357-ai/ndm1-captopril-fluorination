"""
crystal_check.py - baseline contacts in the CRYSTAL pose (original PDB numbering, no MD).
Usage: python crystal_check.py [protein_pdb] [ligand_pdb]
default: 5ZJ2_oxt.pdb  mco_crystal.pdb   (run it inside C:\\NDM1_Project\\md)
"""
import sys
import numpy as np


def parse(path):
    out = []
    for line in open(path, encoding="utf-8-sig", errors="replace"):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[16] not in (" ", "A"):
            continue
        name = line[12:16].strip()
        el = line[76:78].strip().capitalize() or "".join(c for c in name if c.isalpha())[:1].upper()
        out.append((name, line[17:20].strip(), line[22:26].strip(),
                    np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]), el))
    return out


prot = parse(sys.argv[1] if len(sys.argv) > 1 else "5ZJ2_oxt.pdb")
lig = [a for a in parse(sys.argv[2] if len(sys.argv) > 2 else "mco_crystal.pdb") if a[4] != "H"]
zn = [a for a in prot if a[1] == "ZN" or a[4] == "Zn"]
pheavy = [a for a in prot if a[4] not in ("H", "Zn") and a[1] not in ("HOH", "WAT", "ZN")]
print("ligand heavy atoms:", ", ".join("%s(%s)" % (a[0], a[4]) for a in lig))
print("zinc:", ", ".join("ZN%s" % a[2] for a in zn))

print("\nLigand O/S to Zn in the crystal (A):")
for z in zn:
    for a in lig:
        if a[4] in ("O", "S"):
            print("  ZN%-4s %-5s %.2f" % (z[2], a[0], np.linalg.norm(z[3] - a[3])))

print("\nResidues within 4.0 A of ligand (crystal, original numbering):")
res = {}
for p in pheavy:
    d = min(np.linalg.norm(p[3] - a[3]) for a in lig)
    k = "%s %s" % (p[1], p[2])
    res[k] = min(res.get(k, 99.0), d)
for k, d in sorted(res.items(), key=lambda kv: kv[1]):
    if d <= 4.0:
        print("  %-10s %.2f" % (k, d))

print("\nPolar contacts N/O ... ligand N/O <= 3.5 A (distance only):")
for a in lig:
    if a[4] in ("N", "O"):
        for p in pheavy:
            if p[4] in ("N", "O"):
                d = np.linalg.norm(p[3] - a[3])
                if d <= 3.5:
                    print("  %-4s %-14s %.2f" % (a[0], "%s %s:%s" % (p[1], p[2], p[0]), d))