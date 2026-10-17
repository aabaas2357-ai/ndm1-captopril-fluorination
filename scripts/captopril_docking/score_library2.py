import os, re, subprocess, csv
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdMolAlign

ref = Chem.RemoveHs(Chem.MolFromMolFile("captopril_xtal_SH.sdf", removeHs=False))
BOX = ["--center_x","30.5","--center_y","40.5","--center_z","26.7",
       "--size_x","25","--size_y","25","--size_z","25"]
os.makedirs("lib2", exist_ok=True)

def vina(lig, extra, out=None):
    cmd = ["vina","--receptor","protein.pdbqt","--ligand",lig] + BOX + extra
    if out:
        cmd += ["--out", out]
    r = subprocess.run(cmd, capture_output=True, text=True)
    m = re.search(r"Estimated Free Energy of Binding\s*:\s*(-?\d+\.\d+)", r.stdout)
    return float(m.group(1)) if m else None

rows = []
for m in Chem.SDMolSupplier("library_3d.sdf", removeHs=False):
    if m is None:
        continue
    name = m.GetProp("_Name")
    heavy = Chem.RemoveHs(m)
    mcs = rdFMCS.FindMCS([ref, heavy], timeout=10, ringMatchesRingOnly=True,
                         atomCompare=rdFMCS.AtomCompare.CompareElements,
                         bondCompare=rdFMCS.BondCompare.CompareOrder)
    patt = Chem.MolFromSmarts(mcs.smartsString)
    a, b = ref.GetSubstructMatch(patt), heavy.GetSubstructMatch(patt)
    if not a or not b:
        print(name, "no scaffold match"); continue
    amap = [(b[i], a[i]) for i in range(len(a))]
    mh = Chem.Mol(m)
    mh.RemoveAllConformers()
    cids = list(AllChem.EmbedMultipleConfs(mh, 300, randomSeed=7))
    if not cids:
        print(name, "no conformers"); continue
    best = None
    for cid in cids:
        rms = rdMolAlign.AlignMol(mh, ref, prbCid=cid, atomMap=amap)
        if best is None or rms < best[0]:
            best = (rms, cid)
    rms, cid = best
    props = AllChem.MMFFGetMoleculeProperties(mh)
    ff = AllChem.MMFFGetMoleculeForceField(mh, props, confId=cid)
    for idx in b:
        ff.MMFFAddPositionConstraint(idx, 0.2, 100.0)
    ff.Minimize(maxIts=2000)
    sdf, pq = f"lib2/{name}.sdf", f"lib2/{name}.pdbqt"
    Chem.MolToMolFile(mh, sdf, confId=cid)
    subprocess.run(["mk_prepare_ligand.py","-i",sdf,"-o",pq], capture_output=True)
    if not os.path.exists(pq):
        print(name, "pdbqt failed"); continue
    s0 = vina(pq, ["--score_only"])
    s1 = vina(pq, ["--local_only"], f"lib2/{name}_local.pdbqt")
    rows.append((name, len(a), round(rms, 2), s0, s1))
    print(f"{name}: core {len(a)} atoms, core RMSD {rms:.2f} | score_only {s0} | local_only {s1}", flush=True)

with open("library_scores2.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["id", "core_atoms", "core_rmsd", "score_only", "local_only"])
    w.writerows(rows)
