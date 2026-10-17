import os, re, subprocess, csv
import numpy as np
from rdkit import Chem

XTB = os.environ["XTB"]
ROOT = r".\strain"
KCAL = 627.5095
ENV = dict(os.environ, OMP_NUM_THREADS="1")
SETS = [("captopril","SH",-1), ("captopril","Sm",-2), ("F2captopril","SH",-1), ("F2captopril","Sm",-2)]
PAT = Chem.MolFromSmarts("[CX4;!R][CX3](=O)[NX3;R][CX4;R;H1]")

def dih(P, a, b, c, d):
    b0, b1, b2 = P[a]-P[b], P[c]-P[b], P[d]-P[c]
    b1 = b1/np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))

def read_xyz(path):
    L = open(path).read().splitlines()[2:]
    return np.array([[float(x) for x in l.split()[1:4]] for l in L if l.strip()])

ref = {r["set"]: r for r in csv.DictReader(open(os.path.join(ROOT, "strain_results.csv")))}
out = []
for name, form, chrg in SETS:
    tag = f"{name}_{form}"
    mol = Chem.MolFromMolFile(rf".\xtal\{name}_xtal_{form}.sdf", removeHs=False)
    sym = [a.GetSymbol() for a in mol.GetAtoms()]
    m = mol.GetSubstructMatch(PAT)
    om = (m[0], m[1], m[3], m[4])
    P0 = mol.GetConformer().GetPositions()
    d = os.path.join(ROOT, tag + "_relax")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "in.xyz"), "w") as f:
        f.write(f"{len(sym)}\n\n")
        for s, p in zip(sym, P0):
            f.write(f"{s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
    r = subprocess.run([XTB, "in.xyz", "--chrg", str(chrg), "--opt", "tight", "--alpb", "water"],
                       cwd=d, capture_output=True, text=True, encoding="utf-8", errors="replace", env=ENV)
    e = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)", r.stdout)
    if not e or "FAILED TO CONVERGE" in r.stdout:
        print(tag, "relax failed"); continue
    Er = float(e[-1]) * KCAL
    w = dih(read_xyz(os.path.join(d, "xtbopt.xyz")), *om)
    Eb, Ef = float(ref[tag]["E_bound"]), float(ref[tag]["E_free_min"])
    print(f"{tag}: omega_relaxed {w:.1f} | frozen->relaxed {Eb-Er:.2f} | relaxed->free_min {Er-Ef:.2f} | best_trans_rel {ref[tag]['best_trans_rel']}")
