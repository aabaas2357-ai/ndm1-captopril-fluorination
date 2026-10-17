import glob
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

for f in sorted(glob.glob(r".\xtal\*.sdf")):
    m = Chem.MolFromMolFile(f, removeHs=False)
    Chem.AssignStereochemistryFrom3D(m)
    rdCIPLabeler.AssignCIPLabels(m)
    lab = [(a.GetIdx(), a.GetProp("_CIPCode")) for a in m.GetAtoms() if a.HasProp("_CIPCode")]
    print(f.split("\\")[-1], lab)
