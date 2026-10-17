import os, re, subprocess, sys, time
import numpy as np

BASE = r".\cluster"
XTB = os.environ["XTB"]
ENV = dict(os.environ, OMP_NUM_THREADS="4")
KCAL = 627.5095
NSITE = 88
SETS = {"captopril_Sm": -2, "F2captopril_Sm": -2}

def read_xyz(p):
    L = open(p).read().splitlines()
    n = int(L[0])
    sym, xyz = [], []
    for l in L[2:2+n]:
        t = l.split()
        sym.append(t[0]); xyz.append([float(x) for x in t[1:4]])
    return sym, np.array(xyz)

def write_xyz(p, sym, xyz):
    with open(p, "w") as f:
        f.write(f"{len(sym)}\n\n")
        for s, q in zip(sym, xyz):
            f.write(f"{s} {q[0]:.6f} {q[1]:.6f} {q[2]:.6f}\n")

def xtb(d, args, label):
    t = time.time()
    r = subprocess.run([XTB, "in.xyz"] + args, cwd=d, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", env=ENV)
    open(os.path.join(d, "out.txt"), "w", encoding="utf-8").write(r.stdout)
    m = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)", r.stdout)
    ok = bool(m)
    if "--opt" in args:
        ok = ok and "CONVERGED AFTER" in r.stdout and "FAILED TO CONVERGE" not in r.stdout
    print(f"  {label}: {'ok' if ok else 'FAILED'} ({time.time()-t:.0f} s)", flush=True)
    return (float(m[-1]) * KCAL) if ok else None

for tag, lq in SETS.items():
    print(f"\n=== {tag} ===", flush=True)
    sym, xyz = read_xyz(os.path.join(BASE, tag + "_cluster.xyz"))
    fx = open(os.path.join(BASE, tag + ".fix")).read().strip()
    d = os.path.join(BASE, tag + "_run")
    os.makedirs(d, exist_ok=True)
    write_xyz(os.path.join(d, "in.xyz"), sym, xyz)
    open(os.path.join(d, "fix.inp"), "w").write(f"$fix\n atoms: {fx}\n$end\n$opt\n engine=lbfgs\n maxcycle=2000\n$end\n")
    qc = 2 + lq
    Ec = xtb(d, ["--chrg", str(qc), "--opt", "normal", "--alpb", "water", "--input", "fix.inp"], "cluster opt")
    if Ec is None:
        continue
    sym2, xyz2 = read_xyz(os.path.join(d, "xtbopt.xyz"))
    dd = os.path.join(BASE, tag + "_frag")
    os.makedirs(dd, exist_ok=True)
    write_xyz(os.path.join(dd, "in.xyz"), sym2[:NSITE], xyz2[:NSITE])
    Es = xtb(dd, ["--chrg", "2", "--sp", "--alpb", "water"], "site sp")
    write_xyz(os.path.join(dd, "in.xyz"), sym2[NSITE:], xyz2[NSITE:])
    El = xtb(dd, ["--chrg", str(lq), "--sp", "--alpb", "water"], "ligand sp")
    if Es is None or El is None:
        continue
    print(f"  vertical interaction energy: {Ec - Es - El:.2f} kcal/mol", flush=True)
    lig = list(range(NSITE, len(sym2)))
    zn = [i for i in range(NSITE) if sym2[i] == "Zn"]
    S = [i for i in lig if sym2[i] == "S"][0]
    print("  Zn-Zn %.2f A | Zn-S: %s" % (np.linalg.norm(xyz2[zn[0]]-xyz2[zn[1]]),
          ", ".join("%.2f" % np.linalg.norm(xyz2[z]-xyz2[S]) for z in zn)))
    heavy = [i for i in range(NSITE) if sym2[i] not in ("H",)]
    for i in lig:
        if sym2[i] == "F":
            dm = sorted((float(np.linalg.norm(xyz2[i]-xyz2[j])), j) for j in heavy)[:1]
            print("  F atom", i - NSITE + 1, "nearest site heavy atom: %.2f A" % dm[0][0])
