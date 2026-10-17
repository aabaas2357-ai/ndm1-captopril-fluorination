import os
import numpy as np
from rdkit import Chem

BASE = r"."
PDB = os.path.join(BASE, "prot_H.pdb")
OUT = os.path.join(BASE, "cluster")
os.makedirs(OUT, exist_ok=True)
PT = Chem.GetPeriodicTable()
RES = [70, 73, 120, 122, 124, 189, 208, 250]
BACKBONE = {"N","CA","C","O","OXT","H","H1","H2","H3","HA","HA2","HA3"}

atoms = []
for line in open(PDB):
    if not line.startswith(("ATOM", "HETATM")):
        continue
    name = line[12:16].strip()
    resn = line[17:20].strip()
    el = line[76:78].strip().capitalize()
    if not el:
        el = "Zn" if resn == "ZN" else [c for c in name if c.isalpha()][0].upper()
    atoms.append(dict(name=name, resn=resn, chain=line[21].strip(), resi=int(line[22:26]), el=el,
                      xyz=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))

frag, charge = [], 0
for r in RES:
    ra = [a for a in atoms if a["chain"] == "A" and a["resi"] == r]
    if not ra:
        print("residue", r, "NOT FOUND"); continue
    resn = ra[0]["resn"]
    ca = next(a for a in ra if a["name"] == "CA")
    cb = next(a for a in ra if a["name"] == "CB")
    keep = [a for a in ra if a["name"] not in BACKBONE]
    ringH = [a["name"] for a in keep if a["el"] == "H" and not a["name"].startswith("HB")]
    if resn == "CYS":
        print(f"CYS{r}: H before removal = {ringH}; thiolate modelled")
        keep = [a for a in keep if a["name"] != "HG"]
        charge -= 1
    elif resn in ("ASP", "ASH"):
        print(f"ASP{r}: extra H = {ringH}")
        if not ringH: charge -= 1
    elif resn.startswith("HI"):
        print(f"{resn}{r}: ring H = {ringH}")
        if resn == "HIP": charge += 1
    for a in keep:
        frag.append((a["el"], a["xyz"], f"{resn}{r}:{a['name']}"))
    u = cb["xyz"] - ca["xyz"]
    u = u / np.linalg.norm(u)
    frag.append(("H", cb["xyz"] + 1.09 * u, f"cap{r}"))

zn = [a for a in atoms if a["el"] == "Zn"]
print("Zn found:", len(zn))
for a in zn:
    frag.append(("Zn", a["xyz"], f"ZN{a['resi']}")); charge += 2

# N atoms near Zn: should NOT carry H
for el, p, lab in frag:
    if el == "N" and any(np.linalg.norm(p - z["xyz"]) < 2.6 for z in zn):
        nh = sum(1 for e2, p2, _ in frag if e2 == "H" and np.linalg.norm(p - p2) < 1.25)
        print(f"  {lab} near Zn, attached H = {nh}", "  <-- CHECK" if nh else "")

print("protein-site charge:", charge, "| atoms:", len(frag))
LIGS = [("captopril","SH",-1), ("captopril","Sm",-2), ("F2captopril","SH",-1), ("F2captopril","Sm",-2)]
for name, form, lq in LIGS:
    mol = Chem.MolFromMolFile(os.path.join(BASE, "xtal", f"{name}_xtal_{form}.sdf"), removeHs=False)
    P = mol.GetConformer().GetPositions()
    lig = [(a.GetSymbol(), P[a.GetIdx()], f"LIG:{a.GetIdx()}") for a in mol.GetAtoms()]
    allat = frag + lig
    q = charge + lq
    nel = sum(PT.GetAtomicNumber(e) for e, _, _ in allat) - q
    base = f"{name}_{form}"
    with open(os.path.join(OUT, base + "_cluster.xyz"), "w") as f:
        f.write(f"{len(allat)}\n{base} charge {q}\n")
        for e, p, _ in allat:
            f.write(f"{e} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
    fx = [i + 1 for i, (_, _, lab) in enumerate(allat) if lab.endswith(":CB")]
    open(os.path.join(OUT, base + ".fix"), "w").write(",".join(map(str, fx)))
    print(f"{base}: charge {q}, atoms {len(allat)}, electrons {nel} ({'even OK' if nel % 2 == 0 else 'ODD - problem'}), frozen CB {len(fx)}")
    for e, p, lab in lig:
        if e == "F":
            d = sorted((float(np.linalg.norm(p - p2)), lab2) for e2, p2, lab2 in frag if e2 != "H")[:2]
            print("   F", lab, "closest site heavy atoms:", [(round(x, 2), l) for x, l in d])
print("files in", OUT)
