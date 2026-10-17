import os, re, subprocess
import numpy as np

BASE = r".\cluster"
XTB = os.environ["XTB"]
ENV = dict(os.environ, OMP_NUM_THREADS="4")
KCAL = 627.5095
NSITE = 88
FIDX = [NSITE + 17, NSITE + 18]

def read_xyz(p):
    L = open(p).read().splitlines()
    n = int(L[0]); sym, xyz = [], []
    for l in L[2:2+n]:
        t = l.split(); sym.append(t[0]); xyz.append([float(x) for x in t[1:4]])
    return sym, np.array(xyz)

def write_xyz(p, sym, xyz):
    with open(p, "w") as f:
        f.write(f"{len(sym)}\n\n")
        for s, q in zip(sym, xyz):
            f.write(f"{s} {q[0]:.6f} {q[1]:.6f} {q[2]:.6f}\n")

def xtb(d, args):
    r = subprocess.run([XTB, "in.xyz"] + args, cwd=d, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", env=ENV)
    open(os.path.join(d, "out.txt"), "w", encoding="utf-8").write(r.stdout)
    m = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)", r.stdout)
    ok = bool(m) and ("--opt" not in args or ("CONVERGED AFTER" in r.stdout and "FAILED TO CONVERGE" not in r.stdout))
    return float(m[-1]) * KCAL if ok else None

def mk(d, sym, xyz, fx):
    os.makedirs(d, exist_ok=True)
    write_xyz(os.path.join(d, "in.xyz"), sym, xyz)
    open(os.path.join(d, "fix.inp"), "w").write(f"$fix\n atoms: {fx}\n$end\n$opt\n engine=lbfgs\n maxcycle=2000\n$end\n")

symC, xyzC0 = read_xyz(os.path.join(BASE, "captopril_Sm_cluster.xyz"))
symF, xyzF0 = read_xyz(os.path.join(BASE, "F2captopril_Sm_cluster.xyz"))
fx = open(os.path.join(BASE, "captopril_Sm.fix")).read().strip()
same = symC[:NSITE] == symF[:NSITE] and [symC[i] for i in FIDX] == ["H", "H"] and [symF[i] for i in FIDX] == ["F", "F"]
print("atom-order check:", "OK" if same else "MISMATCH - stop")
if not same:
    raise SystemExit
xyzC1 = read_xyz(os.path.join(BASE, "captopril_Sm_run", "xtbopt.xyz"))[1]
xyzF1 = read_xyz(os.path.join(BASE, "F2captopril_Sm_run", "xtbopt.xyz"))[1]

CASES = [
    ("captopril  from own start",      symC, xyzC0),
    ("F2         from own start",      symF, xyzF0),
    ("F2         from captopril geom", symF, xyzC1),
    ("captopril  from F2 geom",        symC, xyzF1),
]
lq = -2
res = []
for i, (label, sym, xyz) in enumerate(CASES):
    d = os.path.join(BASE, f"cross{i}")
    mk(d, sym, xyz, fx)
    Ec = xtb(d, ["--chrg", str(2 + lq), "--opt", "normal", "--alpb", "water", "--input", "fix.inp"])
    if Ec is None:
        print(label, ": cluster opt FAILED", flush=True); continue
    s2, x2 = read_xyz(os.path.join(d, "xtbopt.xyz"))
    ds = d + "_site"; mk(ds, s2[:NSITE], x2[:NSITE], fx)
    Es = xtb(ds, ["--chrg", "2", "--sp", "--alpb", "water"])
    Esr = xtb(ds, ["--chrg", "2", "--opt", "normal", "--alpb", "water", "--input", "fix.inp"])
    dl = d + "_lig"; mk(dl, s2[NSITE:], x2[NSITE:], fx)
    El = xtb(dl, ["--chrg", str(lq), "--sp", "--alpb", "water"])
    if None in (Es, Esr, El):
        print(label, ": a fragment calc FAILED", flush=True); continue
    zn = [j for j in range(NSITE) if s2[j] == "Zn"]
    zzn = np.linalg.norm(x2[zn[0]] - x2[zn[1]])
    ie = Ec - Es - El
    dfm = Es - Esr
    print(f"{label}: IE {ie:.2f} | site deformation {dfm:.2f} | IE+def {ie + dfm:.2f} | Zn-Zn {zzn:.2f}", flush=True)
