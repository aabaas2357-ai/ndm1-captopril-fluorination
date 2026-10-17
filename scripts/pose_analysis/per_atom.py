import glob
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

EL = {"C": "C", "A": "C", "N": "N", "NA": "N", "OA": "O", "SA": "S", "S": "S", "F": "F"}
xyz = lambda l: np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
zn = [xyz(l) for l in open("5ZJ2.pdbqt") if l.startswith(("ATOM", "HETATM")) and l[77:79].strip().upper() == "ZN"]

ref, order = {}, []
for l in open("mco_crystal.pdb"):
    if l.startswith(("ATOM", "HETATM")):
        e = (l[76:78].strip() or l[12:16].strip()[:1]).upper()
        if e != "H":
            nm = l[12:16].strip()
            ref.setdefault(e, []).append((nm, xyz(l)))
            order.append(nm)

def load(f):
    models, cur = [], None
    for l in open(f):
        if l.startswith("MODEL"): cur = {"score": None, "atoms": []}
        elif l.startswith("REMARK VINA RESULT"): cur["score"] = float(l.split()[3])
        elif l.startswith(("ATOM", "HETATM")): cur["atoms"].append((l[77:79].strip(), xyz(l)))
        elif l.startswith("ENDMDL"): models.append(cur)
    return models

def devs(atoms):
    d = {}
    for t, p in atoms:
        if t in EL: d.setdefault(EL[t], []).append(p)
    out, tot, n = {}, 0.0, 0
    for e, lst in ref.items():
        names = [a for a, _ in lst]
        R = np.array([p for _, p in lst])
        c = cdist(np.array(d[e]), R)
        r, k = linear_sum_assignment(c)
        for i, j in zip(r, k):
            out[names[j]] = c[i, j]
        tot += (c[r, k] ** 2).sum(); n += len(r)
    return out, (tot / n) ** 0.5

for name in ("captopril", "F2captopril"):
    best, top = [], []
    for f in sorted(glob.glob(f"seeds/{name}_*.pdbqt")):
        ms = load(f)
        top.append(devs(ms[0]["atoms"])[0])
        cand = []
        for m in ms:
            S = [p for t, p in m["atoms"] if t in ("S", "SA")][0]
            if all(np.linalg.norm(S - z) < 3.0 for z in zn):
                dv, r = devs(m["atoms"]); cand.append((r, dv))
        if cand: best.append(min(cand, key=lambda x: x[0])[1])
    print(f"\n{name}: best-bridging pose in {len(best)}/10 seeds")
    print(f"{'atom':6s} {'best bridging (A)':>18s} {'mode 1 (A)':>12s}")
    for a in order:
        b = np.mean([x[a] for x in best]) if best else float("nan")
        t = np.mean([x[a] for x in top])
        print(f"{a:6s} {b:18.2f} {t:12.2f}")
