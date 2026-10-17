import sys
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

EL = {"C": "C", "A": "C", "N": "N", "NA": "N", "OA": "O", "SA": "S", "S": "S", "F": "F"}
rec, out, xtal = sys.argv[1:4]
xyz = lambda l: np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
zn = [xyz(l) for l in open(rec) if l.startswith(("ATOM", "HETATM")) and l[77:79].strip().upper() == "ZN"]
ref = {}
for l in open(xtal):
    if l.startswith(("ATOM", "HETATM")):
        e = (l[76:78].strip() or l[12:16].strip()[:1]).upper()
        if e != "H": ref.setdefault(e, []).append(xyz(l))
ref = {k: np.array(v) for k, v in ref.items()}

models, cur = [], None
def new(): return {"score": float("nan"), "atoms": []}
for l in open(out):
    if l.startswith("MODEL"):
        cur = new()
    elif l.startswith("REMARK VINA RESULT"):
        if cur is None: cur = new()
        cur["score"] = float(l.split()[3])
    elif l.startswith(("ATOM", "HETATM")):
        if cur is None: cur = new()
        cur["atoms"].append((l[77:79].strip(), xyz(l)))
    elif l.startswith("ENDMDL"):
        models.append(cur); cur = None
if cur is not None and cur["atoms"]:
    models.append(cur)

for i, m in enumerate(models, 1):
    S = [p for t, p in m["atoms"] if t in ("S", "SA")]
    O = [p for t, p in m["atoms"] if t == "OA"]
    s = [min(np.linalg.norm(p - z) for p in S) for z in zn]
    o = [min(np.linalg.norm(p - z) for p in O) for z in zn]
    d = {}
    for t, p in m["atoms"]:
        if t in EL: d.setdefault(EL[t], []).append(p)
    tot = n = 0
    for e in ref:
        c = cdist(np.array(d[e]), ref[e]); r, k = linear_sum_assignment(c)
        tot += (c[r, k] ** 2).sum(); n += len(r)
    print(f"mode {i}: {m['score']:7.3f} | S-Zn {s[0]:5.2f} / {s[1]:5.2f} | O-Zn {o[0]:5.2f} / {o[1]:5.2f} | RMSD vs crystal {(tot / n) ** 0.5:5.2f}")
