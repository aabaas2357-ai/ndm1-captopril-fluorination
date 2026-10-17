"""
analyze_contacts.py - ligand (CAP) contacts in an OpenMM DCD trajectory, using only numpy + openmm.

Usage:  python analyze_contacts.py [dcd] [pdb]
        default: md_5ZJ2/prod_K40000.dcd  md_5ZJ2/minimized.pdb
Distance-only criteria: these are contacts, not proven hydrogen bonds.
"""
import sys, struct, csv, os
import numpy as np

CUT = 4.0       # residue contact cutoff (A), heavy atoms
POL = 3.5       # polar contact cutoff (A)
ZNCUT = 3.0     # ligand-Zn "bound" cutoff (A)
SOLV = {"HOH", "WAT", "NA", "CL", "Na+", "Cl-"}


def read_dcd(path):
    """minimal reader for OpenMM/CHARMM-style DCD. returns xyz (nframes, natoms, 3) in Angstrom and a
    list of per-frame box lengths (or None). Stops quietly at a truncated last frame."""
    with open(path, "rb") as f:
        def rec():
            h = f.read(4)
            if len(h) < 4:
                return None
            n = struct.unpack("<i", h)[0]
            b = f.read(n)
            t = f.read(4)
            if len(b) < n or len(t) < 4:
                return None
            return b

        r1, r2, r3 = rec(), rec(), rec()
        if r1 is None or r2 is None or r3 is None or r1[:4] != b"CORD":
            raise SystemExit("not a readable DCD file: " + path)
        natoms = struct.unpack("<i", r3)[0]
        frames, boxes = [], []
        while True:
            r = rec()
            if r is None:
                break
            box = None
            if len(r) == 48:
                c = struct.unpack("<6d", r)
                lengths = np.array([c[0], c[2], c[5]])
                angs = np.array([c[1], c[3], c[4]])
                ortho = bool(np.all(np.abs(angs) < 1e-3) or np.all(np.abs(angs - 90.0) < 1e-3))
                if lengths.max() < 20.0:      # nm instead of A
                    lengths = lengths * 10.0
                box = lengths if ortho else None
                r = rec()
                if r is None:
                    break
            ys, zs = rec(), rec()
            if ys is None or zs is None:
                break
            x = np.frombuffer(r, dtype="<f4")
            y = np.frombuffer(ys, dtype="<f4")
            z = np.frombuffer(zs, dtype="<f4")
            if not (len(x) == len(y) == len(z) == natoms):
                break
            frames.append(np.stack([x, y, z], axis=1).copy())
            boxes.append(box)
    return np.array(frames, dtype=np.float32), boxes


def mic(d, box):
    return d if box is None else d - box * np.round(d / box)


def main():
    from openmm import app
    dcd = sys.argv[1] if len(sys.argv) > 1 else "md_5ZJ2/prod_K40000.dcd"
    pdbf = sys.argv[2] if len(sys.argv) > 2 else "md_5ZJ2/minimized.pdb"
    name = os.path.splitext(os.path.basename(dcd))[0]

    xyz, boxes = read_dcd(dcd)
    pdb = app.PDBFile(pdbf)
    atoms = list(pdb.topology.atoms())
    nf, na = xyz.shape[0], xyz.shape[1]
    if na != len(atoms):
        raise SystemExit("atom count mismatch: dcd %d vs pdb %d" % (na, len(atoms)))
    print("frames read: %d   atoms: %d   periodic box used: %s"
          % (nf, na, "yes" if boxes and boxes[0] is not None else "NO (non-orthogonal or missing)"))
    if boxes and boxes[0] is not None:
        print("box lengths (A, frame 0): %.1f %.1f %.1f" % tuple(boxes[0]))

    el = np.array([a.element.symbol if a.element is not None else "X" for a in atoms])
    rn = np.array([a.residue.name for a in atoms])
    ri = np.array([a.residue.index for a in atoms])
    an = np.array([a.name for a in atoms])
    label = {a.residue.index: "%s %s" % (a.residue.name, a.residue.id) for a in atoms}

    lig = np.where(rn == "CAP")[0]
    lig_h = lig[el[lig] != "H"]
    zn = np.where(el == "Zn")[0]
    prot = np.where(~np.isin(rn, list(SOLV | {"CAP", "ZN"})) & (el != "H"))[0]
    wat_o = np.where(np.isin(rn, ["HOH", "WAT"]) & (el == "O"))[0]
    if len(lig_h) == 0:
        raise SystemExit("no ligand atoms (residue CAP) found")
    print("ligand heavy atoms:", ", ".join("%s(%s)" % (an[i], el[i]) for i in lig_h))

    lpol = lig_h[np.isin(el[lig_h], ["N", "O"])]
    lzn = lig_h[np.isin(el[lig_h], ["O", "S"])]
    ppol = prot[np.isin(el[prot], ["N", "O"])]

    rp = ri[prot]
    starts = np.r_[0, np.where(np.diff(rp) != 0)[0] + 1]
    rgroup = rp[starts]
    res_min = np.zeros((nf, len(starts)), dtype=np.float32)
    pol_cnt = np.zeros((len(ppol), len(lpol)), dtype=np.int64)
    pol_sum = np.zeros((len(ppol), len(lpol)))
    wat_n = np.zeros((nf, len(lpol)))
    zn_d = np.zeros((nf, len(zn), len(lzn)))

    for f in range(nf):
        X, box = xyz[f].astype(np.float64), boxes[f]
        d = np.linalg.norm(mic(X[prot][:, None, :] - X[lig_h][None, :, :], box), axis=2)
        res_min[f] = np.minimum.reduceat(d.min(1), starts)
        dp = np.linalg.norm(mic(X[ppol][:, None, :] - X[lpol][None, :, :], box), axis=2)
        pol_cnt += (dp <= POL)
        pol_sum += dp
        if len(wat_o):
            dw = np.linalg.norm(mic(X[wat_o][:, None, :] - X[lpol][None, :, :], box), axis=2)
            wat_n[f] = (dw <= POL).sum(0)
        if len(zn):
            zn_d[f] = np.linalg.norm(mic(X[zn][:, None, :] - X[lzn][None, :, :], box), axis=2)

    occ = (res_min <= CUT).mean(0) * 100
    mean_min = res_min.mean(0)
    order = np.argsort(-occ)
    print("\nResidues within %.1f A of the ligand (heavy atoms), top 20:" % CUT)
    print("  %-14s %9s %14s" % ("residue", "occ %", "mean min (A)"))
    for k in order[:20]:
        if occ[k] <= 0:
            break
        print("  %-14s %9.1f %14.2f" % (label[rgroup[k]], occ[k], mean_min[k]))

    print("\nPolar contacts (N/O ... N/O <= %.1f A, distance only), occupancy >= 10%%:" % POL)
    print("  %-8s %-18s %7s %12s" % ("lig", "protein atom", "occ %", "mean d (A)"))
    pairs = []
    for i in range(len(ppol)):
        for j in range(len(lpol)):
            o = 100.0 * pol_cnt[i, j] / nf
            if o >= 10:
                pairs.append((o, an[lpol[j]], "%s:%s" % (label[ri[ppol[i]]], an[ppol[i]]),
                              pol_sum[i, j] / nf))
    for o, la, pa, md in sorted(pairs, reverse=True)[:25]:
        print("  %-8s %-18s %7.1f %12.2f" % (la, pa, o, md))
    if not pairs:
        print("  (none)")

    print("\nWaters within %.1f A of each ligand polar atom (mean count per frame):" % POL)
    for j in range(len(lpol)):
        print("  %-8s %.2f" % (an[lpol[j]], wat_n[:, j].mean()))

    print("\nLigand O/S to Zn: mean distance (A) and %% of frames <= %.1f A:" % ZNCUT)
    for zi in range(len(zn)):
        for j in range(len(lzn)):
            dd = zn_d[:, zi, j]
            if dd.mean() < 6.0:
                print("  Zn%d  %-6s %6.2f  %6.1f %%" % (zn[zi], an[lzn[j]], dd.mean(), 100 * (dd <= ZNCUT).mean()))

    out = "contacts_%s.csv" % name
    with open(out, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["residue", "occupancy_pct", "mean_min_dist_A"])
        for k in order:
            if occ[k] > 0:
                wr.writerow([label[rgroup[k]], round(float(occ[k]), 1), round(float(mean_min[k]), 2)])
    print("\nsaved:", out)


if __name__ == "__main__":
    main()