"""
zn_o_timeline.py - where do the ligand carboxylate O atoms sit relative to Zn, Lys170 NZ and Asn179 N?
Prints the values in the starting structures (solvated / minimized / after equilibration / final)
and every 20th frame (200 ps) of each production DCD. Run inside C:\\NDM1_Project\\md.
Needs analyze_contacts.py in the same folder. MD names: O1, O2 = carboxylate, O3 = amide C=O.
Residue numbers are those of the MD files (Lys170 = Lys211, Asn179 = Asn220 in PDB numbering).
"""
import os
import numpy as np
import analyze_contacts as ac

OUT = "md_5ZJ2"
PDBS = ["solvated.pdb", "minimized.pdb", "prod_start_K40000.pdb", "final_K40000.pdb",
        "prod_start_K40000_r2.pdb", "final_K40000_r2.pdb"]
DCDS = ["prod_K40000.dcd", "prod_K40000_r2.dcd"]


def load_pdb(path):
    atoms, box = [], None
    for line in open(path, encoding="utf-8-sig", errors="replace"):
        if line.startswith("CRYST1"):
            box = np.array([float(line[6:15]), float(line[15:24]), float(line[24:33])])
        elif line.startswith(("ATOM", "HETATM")):
            atoms.append((line[12:16].strip(), line[17:20].strip(), line[22:26].strip(),
                          np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                          line[76:78].strip().capitalize()))
    return atoms, box


def dist(X, box, i, j):
    d = X[i] - X[j]
    if box is not None:
        d = d - box * np.round(d / box)
    return float(np.linalg.norm(d))


ref, _ = load_pdb(os.path.join(OUT, "minimized.pdb"))
lig = {a[0]: i for i, a in enumerate(ref) if a[1] == "CAP"}
zn = [i for i, a in enumerate(ref) if a[4] == "Zn" or a[1] == "ZN"]
lys = [i for i, a in enumerate(ref) if a[1] == "LYS" and a[2] == "170" and a[0] == "NZ"][0]
asn = [i for i, a in enumerate(ref) if a[1] == "ASN" and a[2] == "179" and a[0] == "N"][0]
print("columns (A):  O1-Zn  O2-Zn  O3-Zn | O1-Lys NZ  O2-Lys NZ | O1-Asn N  O2-Asn N  O3-Asn N")
print("(O-Zn = nearest of the two Zn)")


def row(X, box):
    z = [min(dist(X, box, lig[o], k) for k in zn) for o in ("O1", "O2", "O3")]
    k = [dist(X, box, lig[o], lys) for o in ("O1", "O2")]
    n = [dist(X, box, lig[o], asn) for o in ("O1", "O2", "O3")]
    return "%6.2f %6.2f %6.2f | %9.2f %9.2f | %8.2f %8.2f %8.2f" % tuple(z + k + n)


print("\nstructures:")
for f in PDBS:
    p = os.path.join(OUT, f)
    if os.path.exists(p):
        at, box = load_pdb(p)
        X = np.array([a[3] for a in at])
        print("  %-26s %s" % (f, row(X, box)))

for f in DCDS:
    p = os.path.join(OUT, f)
    if not os.path.exists(p):
        continue
    xyz, boxes = ac.read_dcd(p)
    print("\n%s  (%d frames, every 20th = 200 ps):" % (f, len(xyz)))
    for i in range(0, len(xyz), 20):
        print("  frame %3d  %s" % (i, row(xyz[i].astype(float), boxes[i])))