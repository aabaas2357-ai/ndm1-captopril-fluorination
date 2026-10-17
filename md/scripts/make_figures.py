"""
make_figures.py - figures + results table for the README, computed straight from the trajectories.
Run in C:\\NDM1_Project\\md (needs numpy + matplotlib, and analyze_contacts.py next to it).
Writes figures/key_distances.png, figures/ligand_rmsd.png and results_summary.md
'carboxylate O' = the nearer of the two carboxylate oxygens (O1, O2), so the result does not depend on naming.
"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import analyze_contacts as ac

RUNS = [("Gasteiger charges, seed 11", "md_5ZJ2/prod_K40000.dcd", "md_5ZJ2/minimized.pdb", "md_5ZJ2/track_K40000.csv", "tab:blue"),
        ("Gasteiger charges, seed 22", "md_5ZJ2/prod_K40000_r2.dcd", "md_5ZJ2/minimized.pdb", "md_5ZJ2/track_K40000_r2.csv", "tab:purple"),
        ("AM1-BCC charges, seed 11", "md_5ZJ2_bcc/prod_K40000.dcd", "md_5ZJ2_bcc/minimized.pdb", "md_5ZJ2_bcc/track_K40000.csv", "tab:orange")]
CRYSTAL = "md_5ZJ2/solvated.pdb"
LYS, ASN = "170", "179"          # MD numbering = PDB numbering - 41 (Lys211, Asn220)
os.makedirs("figures", exist_ok=True)


def load_pdb(path):
    atoms = []
    for line in open(path, encoding="utf-8-sig", errors="replace"):
        if line.startswith(("ATOM", "HETATM")):
            atoms.append((line[12:16].strip(), line[17:20].strip(), line[22:26].strip(),
                          [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                          line[76:78].strip().capitalize()))
    return atoms


def indices(atoms):
    lig = {a[0]: i for i, a in enumerate(atoms) if a[1] == "CAP"}
    zn = [i for i, a in enumerate(atoms) if a[4] == "Zn" or a[1] == "ZN"]
    lys = [i for i, a in enumerate(atoms) if a[1] == "LYS" and a[2] == LYS and a[0] == "NZ"][0]
    asn = [i for i, a in enumerate(atoms) if a[1] == "ASN" and a[2] == ASN and a[0] == "N"][0]
    return lig, zn, lys, asn


def dist(X, i, j, box):
    d = X[:, i, :].astype(float) - X[:, j, :].astype(float)
    if box is not None:
        d = d - box * np.round(d / box)
    return np.linalg.norm(d, axis=1)


def series(X, boxes, idx):
    lig, zn, lys, asn = idx
    box = None
    if boxes and boxes[0] is not None:
        box = np.array([b if b is not None else [1e9] * 3 for b in boxes], dtype=float)
    ox = [lig["O1"], lig["O2"]]
    zdist = np.min([dist(X, o, z, box) for o in ox for z in zn], axis=0)
    kdist = np.min([dist(X, o, lys, box) for o in ox], axis=0)
    ndist = np.min([dist(X, o, asn, box) for o in ox], axis=0)
    return zdist, kdist, ndist


# crystal reference
cat = load_pdb(CRYSTAL)
cidx = indices(cat)
cX = np.array([[a[3] for a in cat]], dtype=np.float32)
cz, ck, cn = [float(v[0]) for v in series(cX, [None], cidx)]
print("crystal: carboxylate O-Zn %.2f   O-Lys211 NZ %.2f   O-Asn220 N %.2f (A)" % (cz, ck, cn))

data = []
for label, dcd, pdbf, trk, col in RUNS:
    if not (os.path.exists(dcd) and os.path.exists(pdbf)):
        print("skipping (file missing):", label)
        continue
    X, boxes = ac.read_dcd(dcd)
    idx = indices(load_pdb(pdbf))
    z, k, n = series(X, boxes, idx)
    rm = []
    if os.path.exists(trk):
        rows = [r for r in csv.DictReader(open(trk)) if r["stage"] == "prod"]
        rm = [float(r["ligand_RMSD_A"]) for r in rows]
    data.append((label, col, z, k, n, rm))
    print("%-28s frames %d" % (label, len(z)))

# ---- figure 1: key distances
fig, axs = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
panels = [("carboxylate O - nearest Zn (A)", 2, cz), ("carboxylate O - Lys211 NZ (A)", 3, ck),
          ("carboxylate O - Asn220 N (A)", 4, cn)]
for ax, (title, col_i, ref) in zip(axs, panels):
    for d in data:
        t = 0.01 * np.arange(1, len(d[col_i]) + 1)
        ax.plot(t, d[col_i], color=d[1], lw=1.2, label=d[0])
    ax.axhline(ref, color="k", ls="--", lw=1)
    ax.text(0.99, ref, " crystal %.1f" % ref, transform=ax.get_yaxis_transform(), ha="right", va="bottom", fontsize=8)
    ax.set_ylabel(title, fontsize=9)
    ax.set_ylim(0, max(7.5, ref + 1))
axs[0].legend(fontsize=8, loc="center right")
axs[-1].set_xlabel("production time (ns)")
fig.suptitle("D-captopril carboxylate vs crystal pose (5ZJ2), non-bonded Zn model", fontsize=11)
fig.tight_layout()
fig.savefig("figures/key_distances.png", dpi=200)

# ---- figure 2: ligand RMSD
fig2, ax2 = plt.subplots(figsize=(8, 3.2))
for d in data:
    if d[5]:
        ax2.plot(0.005 * np.arange(1, len(d[5]) + 1), d[5], color=d[1], lw=1.0, label=d[0])
ax2.set_xlabel("production time (ns)")
ax2.set_ylabel("ligand RMSD (A)")
ax2.set_title("whole-ligand RMSD stays near 1 A even though the carboxylate moved >3 A", fontsize=10)
ax2.legend(fontsize=8)
fig2.tight_layout()
fig2.savefig("figures/ligand_rmsd.png", dpi=200)

# ---- results table
lines = ["| run | carboxylate O-Zn mean (A) | frames <= 3.0 A | O-Lys211 NZ mean (A) | frames <= 3.5 A | O-Asn220 N mean (A) | frames <= 3.5 A | ligand RMSD mean (A) |",
         "|---|---|---|---|---|---|---|---|",
         "| crystal pose (5ZJ2) | %.2f | - | %.2f | - | %.2f | - | - |" % (cz, ck, cn)]
for label, col, z, k, n, rm in data:
    lines.append("| %s | %.2f | %.0f%% | %.2f | %.0f%% | %.2f | %.0f%% | %s |" % (
        label, z.mean(), 100 * (z <= 3.0).mean(), k.mean(), 100 * (k <= 3.5).mean(),
        n.mean(), 100 * (n <= 3.5).mean(), "%.2f" % np.mean(rm) if rm else "-"))
open("results_summary.md", "w", encoding="utf-8").write("\n".join(lines) + "\n")
print("\n" + "\n".join(lines))
print("\nsaved: figures/key_distances.png, figures/ligand_rmsd.png, results_summary.md")