import numpy as np
xyz = lambda l: np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
at = []
for l in open("mco_crystal.pdb"):
    if l.startswith(("ATOM", "HETATM")):
        e = (l[76:78].strip() or l[12:16].strip()[:1]).upper()
        if e != "H": at.append((l[12:16].strip(), e, xyz(l)))
zn = [xyz(l) for l in open("5ZJ2.pdbqt") if l.startswith(("ATOM","HETATM")) and l[77:79].strip().upper()=="ZN"]
for n, e, p in at:
    nb = [m for m, f, q in at if m != n and np.linalg.norm(p - q) < 1.75]
    dz = [round(float(np.linalg.norm(p - z)), 2) for z in zn]
    print(f"{n:4s} {e:2s} bonded to: {', '.join(nb):14s} dist to Zn: {dz}")
