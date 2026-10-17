import glob
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

EL = {"C": "C", "A": "C", "N": "N", "NA": "N", "OA": "O", "SA": "S", "S": "S", "F": "F"}
xyz = lambda l: np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
zn = [xyz(l) for l in open("5ZJ2.pdbqt") if l.startswith(("ATOM", "HETATM")) and l[77:79].strip().upper() == "ZN"]
ref = {}
for l in open("mco_crystal.pdb"):
    if l.startswith(("ATOM", "HETATM")):
        e = (l[76:78].strip() or l[12:16].strip()[:1]).upper()
        if e != "H": ref.setdefault(e, []).append(xyz(l))
ref = {k: np.array(v) for k, v in ref.items()}

def load(f):
    models, cur = [], None
    for l in open(f):
        if l.startswith("MODEL"): cur = {"score": None, "atoms": []}
        elif l.startswith("REMARK VINA RESULT"): cur["score"] = float(l.split()[3])
        elif l.startswith(("ATOM", "HETATM")): cur["atoms"].append((l[77:79].strip(), xyz(l)))
        elif l.startswith("ENDMDL"): models.append(cur)
    return models

def rmsd(atoms):
    d = {}
    for t, p in atoms:
        if t in EL: d.setdefault(EL[t], []).append(p)
    tot = n = 0
    for e in ref:
        c = cdist(np.array(d[e]), ref[e]); r, k = linear_sum_assignment(c)
        tot += (c[r, k] ** 2).sum(); n += len(r)
    return (tot / n) ** 0.5

for name in ("captopril", "F2captopril"):
    best_all, best_br, found = [], [], []
    for f in sorted(glob.glob(f"seeds/{name}_*.pdbqt")):
        ms = load(f); top = ms[0]["score"]
        rows = []
        for i, m in enumerate(ms, 1):
            S = [p for t, p in m["atoms"] if t in ("S", "SA")][0]
            br = all(np.linalg.norm(S - z) < 3.0 for z in zn)
            rows.append((rmsd(m["atoms"]), br, i, m["score"] - top))
        best_all.append(min(r[0] for r in rows))
        b = [r for r in rows if r[1]]
        if b: best_br.append(min(r[0] for r in b))
        c = [r for r in rows if r[1] and r[0] < 2.0]
        if c: found.append(min(c, key=lambda r: r[3]))
    print(f"{name}: lowest RMSD in any pose {np.mean(best_all):.2f} | lowest RMSD among bridging poses {np.mean(best_br):.2f} +/- {np.std(best_br):.2f}")
    print(f"   runs with bridging AND RMSD<2.0: {len(found)}/10" + (f" | rank {np.mean([r[2] for r in found]):.1f}, gap to top {np.mean([r[3] for r in found]):.2f} kcal/mol" if found else ""))
