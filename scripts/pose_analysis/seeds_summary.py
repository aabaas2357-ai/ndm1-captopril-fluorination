import glob
import numpy as np
zn = [np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) for l in open("5ZJ2.pdbqt")
      if l.startswith(("ATOM", "HETATM")) and l[77:79].strip().upper() == "ZN"]
for name in ("captopril", "F2captopril"):
    tops, bridges, ranks, gaps = [], [], [], []
    for f in sorted(glob.glob(f"seeds/{name}_*.pdbqt")):
        models, cur = [], None
        for l in open(f):
            if l.startswith("MODEL"): cur = {"s": None, "S": None}
            elif l.startswith("REMARK VINA RESULT"): cur["s"] = float(l.split()[3])
            elif l.startswith(("ATOM", "HETATM")) and l[77:79].strip() in ("S", "SA"):
                cur["S"] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            elif l.startswith("ENDMDL"): models.append(cur)
        top = models[0]["s"]
        br = [(i + 1, m["s"]) for i, m in enumerate(models) if all(np.linalg.norm(m["S"] - z) < 3.0 for z in zn)]
        tops.append(top)
        if br:
            ranks.append(br[0][0]); bridges.append(br[0][1]); gaps.append(br[0][1] - top)
    n = len(tops)
    print(f"{name}: runs {n} | top energy {np.mean(tops):.2f} +/- {np.std(tops):.2f} | runs with a bridging pose {len(bridges)}/{n}")
    if bridges:
        print(f"   best bridging pose: energy {np.mean(bridges):.2f} +/- {np.std(bridges):.2f}, rank {np.mean(ranks):.1f}, gap to top {np.mean(gaps):.2f} +/- {np.std(gaps):.2f} kcal/mol")
