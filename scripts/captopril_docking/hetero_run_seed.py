import os, re, csv, subprocess
import numpy as np
from multiprocessing import Pool
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdMolAlign, Descriptors, Crippen, rdMolDescriptors, rdCIPLabeler
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from meeko import PDBQTMolecule, RDKitMolCreate

import sys
SEED = int(sys.argv[1])
OUT = f"lib_s{SEED}"
os.makedirs(OUT, exist_ok=True)
REF = Chem.RemoveHs(Chem.MolFromMolFile("captopril_xtal_SH.sdf", removeHs=False))
BOX = ["--center_x","30.5","--center_y","40.5","--center_z","26.7",
       "--size_x","25","--size_y","25","--size_z","25"]

PH, ZN = [], []
for l in open("protein.pdbqt"):
    if l.startswith(("ATOM", "HETATM")):
        t = l.split()[-1].upper()
        xyz = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
        if t == "ZN":
            ZN.append(xyz)
        elif not t.startswith("H"):
            PH.append(xyz)
PH, ZN = np.array(PH), np.array(ZN)

LIB = [
 ("ref_Dcap",      "C[C@H](CS)C(=O)N1CCC[C@@H]1C(=O)O", "ref"),
 ("ref_CF2_44diF", "C[C@H](CS)C(=O)N1CC(F)(F)C[C@@H]1C(=O)O", "ref"),
 ("ref_CF3chain",  "FC(F)(F)[C@H](CS)C(=O)N1CCC[C@@H]1C(=O)O", "ref"),
 ("ref_CF3_44diF", "FC(F)(F)[C@H](CS)C(=O)N1CC(F)(F)C[C@@H]1C(=O)O", "ref"),
 ("h1_tetrazole",  "C[C@H](CS)C(=O)N1CCC[C@@H]1c1nnn[nH]1", "hetero"),
 ("h2_isonip_tz",  "C[C@H](CS)C(=O)N1CCC(CC1)c1nnn[nH]1", "hetero"),
 ("h3_thioPro",    "C[C@H](CS)C(=O)N1CSC[C@@H]1C(=O)O", "hetero"),
 ("h4_indoline",   "C[C@H](CS)C(=O)N1c2ccccc2C[C@@H]1C(=O)O", "hetero"),
 ("h5_thiazole",   "C[C@H](CS)C(=O)Nc1ncc(s1)C(=O)O", "hetero"),
 ("val_Lcap",      "C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O", "val"),
]
VAL_IC50 = {"val_RR": 0.3, "val_RS": 4.5, "val_SR": 6.1, "val_SS": 60.0, "val_Lcap": 7.0}

def cf3_isomers():
    base = Chem.MolFromSmiles("FC(F)(F)C(CS)C(=O)N1CCCC1C(=O)O")
    res = {}
    for iso in EnumerateStereoisomers(base):
        rdCIPLabeler.AssignCIPLabels(iso)
        ring = chain = None
        for a in iso.GetAtoms():
            if a.HasProp("_CIPCode"):
                if a.IsInRing():
                    ring = a.GetProp("_CIPCode")
                else:
                    chain = a.GetProp("_CIPCode")
        res[(ring, chain)] = Chem.MolToSmiles(iso)
    return res

def vina(lig, extra, out=None):
    cmd = ["vina", "--receptor", "protein.pdbqt", "--ligand", lig, "--cpu", "1"] + BOX + extra
    if out:
        cmd += ["--out", out]
    r = subprocess.run(cmd, capture_output=True, text=True)
    m = re.search(r"Estimated Free Energy of Binding\s*:\s*(-?\d+\.\d+)", r.stdout)
    return float(m.group(1)) if m else None

def process(task):
    name, smi, tier = task
    r = {"id": name, "tier": tier, "smiles": smi}
    try:
        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        heavy = Chem.RemoveHs(m)
        r["MW"] = round(Descriptors.MolWt(heavy), 1)
        r["cLogP"] = round(Crippen.MolLogP(heavy), 2)
        r["TPSA"] = round(rdMolDescriptors.CalcTPSA(heavy), 1)
        mcs = rdFMCS.FindMCS([REF, heavy], timeout=10, ringMatchesRingOnly=True,
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareOrder)
        patt = Chem.MolFromSmarts(mcs.smartsString)
        a, b = REF.GetSubstructMatch(patt), heavy.GetSubstructMatch(patt)
        if not a or not b:
            r["note"] = "no scaffold match"
            return r
        r["core_atoms"] = len(a)
        amap = [(b[i], a[i]) for i in range(len(a))]
        ps = AllChem.ETKDGv3()
        ps.randomSeed = SEED
        ps.numThreads = 1
        cids = list(AllChem.EmbedMultipleConfs(m, 300, ps))
        if not cids:
            r["note"] = "no conformers"
            return r
        best = None
        for cid in cids:
            rms = rdMolAlign.AlignMol(m, REF, prbCid=cid, atomMap=amap)
            if best is None or rms < best[0]:
                best = (rms, cid)
        rms, cid = best
        r["core_rmsd"] = round(rms, 2)
        props = AllChem.MMFFGetMoleculeProperties(m)
        if props is not None:
            ff = AllChem.MMFFGetMoleculeForceField(m, props, confId=cid)
            for idx in b:
                ff.MMFFAddPositionConstraint(idx, 0.2, 100.0)
        else:
            ff = AllChem.UFFGetMoleculeForceField(m, confId=cid)
            for idx in b:
                ff.UFFAddPositionConstraint(idx, 0.2, 100.0)
        ff.Minimize(maxIts=2000)
        sdf, pq = f"{OUT}/{name}.sdf", f"{OUT}/{name}.pdbqt"
        Chem.MolToMolFile(m, sdf, confId=cid)
        subprocess.run(["mk_prepare_ligand.py", "-i", sdf, "-o", pq], capture_output=True)
        if not os.path.exists(pq):
            r["note"] = "pdbqt failed"
            return r
        r["score_only"] = vina(pq, ["--score_only"])
        loc = f"{OUT}/{name}_local.pdbqt"
        r["local_only"] = vina(pq, ["--local_only"], loc)
        try:
            pm = PDBQTMolecule.from_file(loc, is_dlg=False, skip_typing=True)
            lm = RDKitMolCreate.from_pdbqt_mol(pm)[0]
            Chem.MolToMolFile(lm, f"{OUT}/{name}_local.sdf")
            pos = lm.GetConformer().GetPositions()
            sy = [x.GetSymbol() for x in lm.GetAtoms()]
            S = [i for i, s in enumerate(sy) if s == "S"]
            if S and len(ZN):
                r["S_Zn_min"] = round(min(float(np.linalg.norm(pos[i] - z)) for i in S for z in ZN), 2)
            hv = np.array([pos[i] for i, s in enumerate(sy) if s != "H"])
            d = np.linalg.norm(hv[:, None, :] - PH[None, :, :], axis=2)
            r["min_contact"] = round(float(d.min()), 2)
        except Exception as e:
            r["note"] = "local pose parse failed: " + str(e)[:60]
    except Exception as e:
        r["note"] = "error: " + str(e)[:80]
    return r

def ranks(x):
    o = np.argsort(x)
    r = np.empty(len(x))
    r[o] = np.arange(len(x))
    return r

def spearman(a, b):
    return float(np.corrcoef(ranks(np.array(a)), ranks(np.array(b)))[0, 1])

if __name__ == "__main__":
    iso = cf3_isomers()
    tasks = list(LIB)
    for nm, key in [("val_RR", ("R", "R")), ("val_RS", ("R", "S")),
                    ("val_SR", ("S", "R")), ("val_SS", ("S", "S"))]:
        if key in iso:
            tasks.append((nm, iso[key], "val"))
        else:
            print("missing CF3 isomer", key)
    nproc = max(1, (os.cpu_count() or 2) - 1)
    print(f"{len(tasks)} compounds on {nproc} workers", flush=True)
    rows = []
    with Pool(nproc) as pool:
        for r in pool.imap_unordered(process, tasks):
            print("done:", r["id"], r.get("local_only"), r.get("note", ""), flush=True)
            rows.append(r)
    by = {r["id"]: r for r in rows}
    base = by.get("ref_Dcap", {})
    for r in rows:
        for k in ("score_only", "local_only"):
            if r.get(k) is not None and base.get(k) is not None:
                r["d_" + k] = round(r[k] - base[k], 2)
    keys = ["id", "tier", "core_atoms", "core_rmsd", "score_only", "local_only", "d_score_only",
            "d_local_only", "S_Zn_min", "min_contact", "MW", "cLogP", "TPSA", "note", "smiles"]
    with open(f"library_scores_s{SEED}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["tier"], r["id"])))
    print("\nid | core_rmsd | score_only | local_only | d_local | S-Zn | min_contact | cLogP")
    for r in sorted(rows, key=lambda r: (r["tier"], r.get("local_only") or 0)):
        print(r["id"], r.get("core_rmsd"), r.get("score_only"), r.get("local_only"),
              r.get("d_local_only"), r.get("S_Zn_min"), r.get("min_contact"), r.get("cLogP"), r.get("note", ""))
    v = [(k, by[k]) for k in VAL_IC50 if k in by and by[k].get("local_only") is not None]
    if len(v) >= 4:
        ic = [VAL_IC50[k] for k, _ in v]
        for key in ("local_only", "score_only"):
            sc = [r[key] for _, r in v]
            print(f"VALIDATION {key}: Spearman rho vs IC50 = {spearman(sc, ic):.2f} (n={len(v)}; expected positive)")
        print("IC50 order (expected best to worst):", [k for k, _ in sorted(v, key=lambda t: VAL_IC50[t[0]])])
        print("local_only order (best to worst):", [k for k, _ in sorted(v, key=lambda t: t[1]["local_only"])])
