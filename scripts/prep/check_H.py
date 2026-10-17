import math
def heavy(path):
    d = {}
    for l in open(path):
        if l.startswith(("ATOM", "HETATM")):
            n = l[12:16].strip()
            if n[0] in "CNOSZ":
                d[(l[21], int(l[22:26]), n)] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return d
a, b = heavy("5ZJ2_clean.pdb"), heavy("5ZJ2_H.pdb")
common = a.keys() & b.keys()
dev = sorted(((math.dist(a[k], b[k]), k) for k in common), reverse=True)
print("heavy atoms: clean", len(a), "| H-added", len(b), "| common", len(common))
print("moved > 0.10 A:", sum(1 for d, k in dev if d > 0.10), "| max shift %.2f A at %s" % dev[0])
zn = [v for k, v in b.items() if k[2] == "ZN"]
for k, v in b.items():
    if k[2] in ("OD2", "SG", "NE2", "ND1") and any(math.dist(v, z) < 2.7 for z in zn):
        print("coordinating", k, "%.2f A" % min(math.dist(v, z) for z in zn))
