import os, re, shutil, subprocess, csv
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms as T

XTB = os.environ["XTB"]
ROOT = r".\strain"
os.makedirs(ROOT, exist_ok=True)
SETS = [("captopril","SH",-1), ("captopril","Sm",-2), ("F2captopril","SH",-1), ("F2captopril","Sm",-2)]
NCONF = 60
WORKERS = max(1, (os.cpu_count() or 2) // 2)
ENV = dict(os.environ, OMP_NUM_THREADS="1")
KCAL = 627.5095

def write_xyz(path, symbols, xyz):
    with open(path, "w") as f:
        f.write(f"{len(symbols)}\n\n")
        for s, p in zip(symbols, xyz):
            f.write(f"{s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")

def read_xyz(path):
    L = open(path).read().splitlines()[2:]
    return np.array([[float(x) for x in l.split()[1:4]] for l in L if l.strip()])

def dih(P, a, b, c, d):
    b0, b1, b2 = P[a]-P[b], P[c]-P[b], P[d]-P[c]
    b1 = b1/np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))

def run_xtb(tag, symbols, xyz, chrg, fix=None):
    d = os.path.join(ROOT, tag)
    os.makedirs(d, exist_ok=True)
    write_xyz(os.path.join(d, "in.xyz"), symbols, xyz)
    cmd = [XTB, "in.xyz", "--chrg", str(chrg), "--opt", "tight", "--alpb", "water"]
    if fix:
        with open(os.path.join(d, "fix.inp"), "w") as f:
            f.write("$fix\n atoms: " + ",".join(str(i+1) for i in fix) + "\n$end\n$opt\n engine=lbfgs\n maxcycle=2000\n$end\n")
        cmd += ["--input", "fix.inp"]
    r = subprocess.run(cmd, cwd=d, capture_output=True, text=True, encoding="utf-8", errors="replace", env=ENV)
    m = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)", r.stdout)
    ok = m and os.path.exists(os.path.join(d, "xtbopt.xyz")) and "FAILED TO CONVERGE" not in r.stdout
    if not ok:
        return None, None
    return float(m[-1]) * KCAL, read_xyz(os.path.join(d, "xtbopt.xyz"))

PAT = Chem.MolFromSmarts("[CX4;!R][CX3](=O)[NX3;R][CX4;R;H1]")
rows = []
for name, form, chrg in SETS:
    sdf = rf".\xtal\{name}_xtal_{form}.sdf"
    mol = Chem.MolFromMolFile(sdf, removeHs=False)
    if mol is None:
        print("cannot read", sdf); continue
    sym = [a.GetSymbol() for a in mol.GetAtoms()]
    fix = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1 and a.GetSymbol() != "F"]
    m = mol.GetSubstructMatch(PAT)
    if not m:
        print("amide pattern not found in", sdf); continue
    om = (m[0], m[1], m[3], m[4])
    tag = f"{name}_{form}"
    print(f"\n=== {tag} (charge {chrg}) ===", flush=True)

    P0 = mol.GetConformer().GetPositions()
    print(f"crystal omega = {dih(P0, *om):.1f} deg")
    Eb, Pb = run_xtb(tag + "_bound", sym, P0, chrg, fix)
    if Eb is None:
        print("bound optimisation failed"); continue

    m2 = Chem.Mol(mol)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 7
    ps.pruneRmsThresh = 0.5
    cids = list(AllChem.EmbedMultipleConfs(m2, NCONF, ps))
    jobs = []
    for cid in cids:
        P = m2.GetConformer(cid).GetPositions()
        jobs.append((f"{tag}_f{cid}", P))
        c2 = Chem.Conformer(m2.GetConformer(cid))
        w = dih(P, *om)
        T.SetDihedralDeg(c2, *om, 0.0 if abs(w) > 90 else 180.0)
        jobs.append((f"{tag}_f{cid}x", c2.GetPositions()))
    print(f"{len(jobs)} free starting structures, {WORKERS} workers", flush=True)

    with ThreadPoolExecutor(WORKERS) as ex:
        res = list(ex.map(lambda j: run_xtb(j[0], sym, j[1], chrg), jobs))
    free = []
    for (t, _), (E, P) in zip(jobs, res):
        if E is not None:
            free.append((E, abs(dih(P, *om)) < 90, t))
    free.sort()
    if not free:
        print("no free conformer converged"); continue
    Ef = free[0][0]
    cis = [f for f in free if f[1]]
    trans = [f for f in free if not f[1]]
    print(f"converged {len(free)}/{len(jobs)}")
    print(f"E_bound {Eb:.2f} | E_free_min {Ef:.2f} ({'cis' if free[0][1] else 'trans'}) | STRAIN {Eb-Ef:.2f} kcal/mol")
    if cis and trans:
        print(f"best trans {trans[0][0]-Ef:.2f} | best cis {cis[0][0]-Ef:.2f} kcal/mol relative to free min")
    rows.append([tag, round(Eb-Ef, 2), round(Eb, 2), round(Ef, 2),
                 round(trans[0][0]-Ef, 2) if trans else "", round(cis[0][0]-Ef, 2) if cis else ""])

with open(os.path.join(ROOT, "strain_results.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["set", "strain_kcal", "E_bound", "E_free_min", "best_trans_rel", "best_cis_rel"])
    w.writerows(rows)
print("\nsaved strain\\strain_results.csv")


