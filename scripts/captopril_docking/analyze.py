import glob, re, numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from meeko import PDBQTMolecule, RDKitMolCreate

zn = []
for l in open("protein.pdbqt"):
    if l.startswith(("ATOM", "HETATM")) and l.split()[-1].upper() == "ZN":
        zn.append(np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
for lig in ["captopril", "F2captopril"]:
    ref = Chem.RemoveHs(Chem.MolFromMolFile(f"{lig}_xtal_SH.sdf", removeHs=False))
    for f in sorted(glob.glob(f"{lig}_s[0-9].pdbqt")):
        sc = [float(m) for m in re.findall(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", open(f).read())]
        mol = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule.from_file(f, is_dlg=False, skip_typing=True))[0]
        molH = Chem.RemoveHs(mol)
        rms = [rdMolAlign.CalcRMS(molH, ref, prbId=i) for i in range(molH.GetNumConformers())]
        s = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "S"][0]
        p = mol.GetConformer(0).GetAtomPosition(s); p = np.array([p.x, p.y, p.z])
        d = min(np.linalg.norm(p - z) for z in zn)
        print(f"{f}: top {sc[0]:.2f} | RMSD pose1 {rms[0]:.2f} | best RMSD in top10 {min(rms):.2f} | pose1 S-Zn {d:.2f}")
