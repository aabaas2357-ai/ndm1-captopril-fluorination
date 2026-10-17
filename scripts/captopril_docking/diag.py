import numpy as np
from rdkit import Chem
from meeko import PDBQTMolecule, RDKitMolCreate

def coords(path, key):
    out = []
    for l in open(path):
        if l.startswith(("ATOM", "HETATM")) and l.split()[-1] == key:
            out.append(np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
    return out

zn, tz = coords("protein_tz.pdbqt", "ZN"), coords("protein_tz.pdbqt", "TZ")
ref = Chem.MolFromMolFile("captopril_xtal_SH.sdf", removeHs=False)
sref = [ref.GetConformer().GetAtomPosition(a.GetIdx()) for a in ref.GetAtoms() if a.GetSymbol() == "S"][0]
sref = np.array([sref.x, sref.y, sref.z])
fmt = lambda p, L: " ".join("%.2f" % np.linalg.norm(p - q) for q in L)
print("crystal S -> Zn:", fmt(sref, zn), "| -> TZ:", fmt(sref, tz))
pm = PDBQTMolecule.from_file("captopril_ad4_out.pdbqt", is_dlg=False, skip_typing=True)
mol = RDKitMolCreate.from_pdbqt_mol(pm)[0]
sidx = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "S"][0]
for i in range(min(5, mol.GetNumConformers())):
    p = mol.GetConformer(i).GetAtomPosition(sidx)
    p = np.array([p.x, p.y, p.z])
    print(f"pose {i+1}: S -> Zn: {fmt(p, zn)} | -> TZ: {fmt(p, tz)}")
