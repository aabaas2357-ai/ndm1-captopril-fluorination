import sys
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from meeko import PDBQTMolecule, RDKitMolCreate

out, ref_sdf = sys.argv[1], sys.argv[2]
ref = Chem.RemoveHs(Chem.MolFromMolFile(ref_sdf, removeHs=False))
pm = PDBQTMolecule.from_file(out, is_dlg=False, skip_typing=True)
mol = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pm)[0])
for i in range(mol.GetNumConformers()):
    try:
        r = rdMolAlign.CalcRMS(mol, ref, prbId=i)
        print(f"pose {i+1}: RMSD to crystal = {r:.2f} A")
    except Exception as e:
        print(f"pose {i+1}: RMSD failed -> {e}")
