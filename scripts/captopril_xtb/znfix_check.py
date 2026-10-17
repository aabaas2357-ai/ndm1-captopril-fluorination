# NOTE: this script failed its own sanity check (frozen Zn atoms moved); its results are NOT valid. Kept for transparency.
import os, re, subprocess
import numpy as np

BASE = r".\cluster"
XTB = os.environ["XTB"]
ENV = dict(os.environ, OMP_NUM_THREADS="4")
KCAL = 627.5095
NSITE = 88

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
    m = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)", r.stdout)
    ok = bool(m) and ("--opt" not in args or ("CONVERGED AFTER" in r.stdout and "FAILED TO CONVERGE" not in r.stdout))
    return float(m[-1]) * KCAL if ok else None

fx0 = [int(x) for x in open(os.path.join(BASE, "captopril_Sm.fix")).read().strip().split(",")]
for tag in ["captopril_Sm", "F2captopril_Sm"]:
    sym, xyz = read_xyz(os.path.join(BASE, tag + "_cluster.xyz"))
    zn = [i + 1 for i in range(NSITE) if sym[i] == "Zn"]
    fx = sorted(set(fx0 + zn))
    d = os.path.join(BASE, tag + "_znfix"); os.makedirs(d, exist_ok=True)
    write_xyz(os.path.join(d, "in.xyz"), sym, xyz)
    open(os.path.join(d, "fix.inp"), "w").write("$fix\n atoms: " + ",".join(map(str, fx)) + "\n$end\n$opt\n engine=lbfgs\n maxcycle=2000\n$end\n")
    Ec = xtb(d, ["--chrg", "0", "--opt", "normal", "--alpb", "water", "--input", "fix.inp"])
    if Ec is None:
        print(tag, "cluster opt FAILED"); continue
    s2, x2 = read_xyz(os.path.join(d, "xtbopt.xyz"))
    for sub, a, b, q in [("site", 0, NSITE, 2), ("lig", NSITE, len(s2), -2)]:
        dd = os.path.join(BASE, tag + "_znfix_" + sub); os.makedirs(dd, exist_ok=True)
        write_xyz(os.path.join(dd, "in.xyz"), s2[a:b], x2[a:b])
    Es = xtb(os.path.join(BASE, tag + "_znfix_site"), ["--chrg", "2", "--sp", "--alpb", "water"])
    El = xtb(os.path.join(BASE, tag + "_znfix_lig"), ["--chrg", "-2", "--sp", "--alpb", "water"])
    zi = [i for i in range(NSITE) if s2[i] == "Zn"]
    S = [i for i in range(NSITE, len(s2)) if s2[i] == "S"][0]
    print(f"{tag}: IE {Ec - Es - El:.2f} kcal/mol | Zn-Zn {np.linalg.norm(x2[zi[0]]-x2[zi[1]]):.2f} | Zn-S {np.linalg.norm(x2[zi[0]]-x2[S]):.2f}, {np.linalg.norm(x2[zi[1]]-x2[S]):.2f}", flush=True)
